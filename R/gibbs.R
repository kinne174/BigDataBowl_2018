#To do 5_9
# -Only look at plays that occur between 10-20 yard lines
# -implement full gibbs sampling using a wishart prior for V^-1
# 

library(ks)
library(readr)
library(Matrix)
library(MASS)
library(MCMCpack)

#simulate data
simulation.data = function(n_samples, n_active, n_preplay){
  
  preplays.array = matrix(NA, nrow = n_samples, ncol = n_preplay)
  actives.array = matrix(NA, nrow = n_samples, ncol = n_active)
  successes.array = numeric(n_samples)
  
  pre.prob = runif(11)
  prob.vec = pre.prob/sum(pre.prob)
  
  true.beta = matrix(rnorm(n_active*(n_preplay + 1)), nrow = n_active, ncol = n_preplay+1)
  true.W = rnorm(n_active)
  
  for(i in 1:n_samples){
    preplays.array[i,] = c(c(rmultinom(n = 3, size = 1, prob = prob.vec))[-c(1, 12, 23)], rnorm(n=12))
    actives.array[i,] = c(true.beta %*% cbind(c(1, preplays.array[i,]))) + rnorm(n_active)
    successes.array[i] = 1*((1 + exp(sum(true.W * actives.array[i,])))**-1 > 0.5)
  }
  
  return(list(preplay = preplays.array, active = actives.array, success = successes.array, true.beta = true.beta))
}

set.seed(635)
simulate.out = simulation.data(n_samples = 5000, n_active = 16, n_preplay = 30+12)
mean(simulate.out$success)
a = lm(simulate.out$active ~ simulate.out$preplay)

real.data = function(){
  setwd('C:/Users/Mitch/Documents/UofM/Fall 2018/NFL/Data/my_data')
  
  success = c(as.matrix(read_csv('success.csv', col_names = F)))
  active = as.matrix(read_csv('active.csv', col_names = F))
  preplay = as.matrix(read_csv('preplay.csv', col_names = F))
  
  preplay = preplay[,c(-14, -25, -36)]
  
  return(list(preplay = preplay, active = active, success = success))
}

real.out = real.data()
mean(real.out$success)

#braindead
braindead = glm(out$success ~ out$preplay, family='binomial')

pred = 1*(predict(braindead, newdata=as.data.frame(out$preplay)) > 0.5)
sum((out$success == pred))/length(out$success)


full.gibbs = function(chain.length, data.list, quiet=FALSE){
  preplay = add_interaction_to_preplay(cbind(1, data.list$preplay))
  active_response = data.list$active
  active_explanatory = cbind(1, data.list$active)
  success = data.list$success
  
  #get estimate of beta using multivariate linear regression
  init.beta = t(as.matrix(lm(active_response ~ preplay + 0)$coefficients))
  
  #get estimate of V using covariance matrix of actives, unclear how to get prior for degrees of freedom
  # init.phi = cov(active_response)
  # init.V = cov(active_response)
  init.phi = diag(ncol(active_response))
  init.V = diag(ncol(active_response))
  init.df = ncol(active_response)*2
  
  init.z = rep(0, length(success))
  
  #initialize chains for beta, lambda, and V, might be easiest to vectorize beta and V then transform back when needed
  beta.chain = matrix(NA, nrow = ncol(preplay)*ncol(active_response), ncol = chain.length)
  V.chain = matrix(NA, nrow = (ncol(active_response)*(ncol(active_response) + 1))/2, ncol = chain.length)
  lambda.chain = numeric(chain.length)
  
  z.chain = matrix(NA, nrow = length(success), ncol = chain.length)
  W.chain = matrix(NA, nrow = ncol(active_explanatory), ncol = chain.length)
  
  #Start chains for beta (with prior estimate) and V (with prior estimate)
  # start gibbs sampling with lambda draw
  beta.chain[,1] = c(init.beta)
  V.chain[,1] = vech(init.V)
  
  z.chain[,1] = init.z
  
  #define preplay_swirl
  preplay_swirl = matrix(0, nrow = nrow(preplay)*ncol(active_response), ncol = ncol(preplay)*ncol(active_response))
  for(i in 1:nrow(preplay)){
    temp_mat = matrix(0, nrow = ncol(active_response), ncol = ncol(preplay))
    temp_mat[1,] = preplay[i,]
    temp_vec = c(temp_mat)
    for(j in 1:ncol(active_response)){
      preplay_swirl[(i-1)*ncol(active_response) + j,] = temp_vec
      temp_vec = c(temp_vec[length(temp_vec)], temp_vec[-length(temp_vec)])
    }
  }
  
  preplay_swirl = Matrix(preplay_swirl, sparse = T)
  
  #functions to calculate V* and M*
  V.star = function(preplay_swirl, V.inv){
    #should return a ncol(preplay)*ncol(active) x ncol(preplay)*ncol(active) matrix
    V.inv.kronk = kronecker(Matrix(diag(nrow(preplay_swirl)/ncol(V.inv)), sparse=T), V.inv)
    
    temp = t(preplay_swirl) %*% V.inv.kronk %*% preplay_swirl
    
    return(diag(ncol(preplay_swirl)) + temp)
  }
  
  M.star = function(preplay_swirl, active_response, V.inv, init.beta){
    #should return a ncol(preplay)*ncol(active) dimension vector (not column vector)
    V.inv.kronk = kronecker(Matrix(diag(nrow(preplay_swirl)/ncol(V.inv)), sparse=T), V.inv)
    
    temp = t(preplay_swirl) %*% V.inv.kronk %*% cbind(c(t(active_response)))
    
    return(c(init.beta) + temp)
  }
  
  #functions to calculate shape, scale parameters of lambda posterior
  shape.lambda = function(active_response, preplay){
    return((nrow(active_response)*ncol(active_response) + ncol(preplay)*ncol(active_response))*.5)
  }
  
  rate.lambda = function(beta.vec, preplay, active_response, V.inv, beta.init.vec){
    beta.mat = matrix(beta.vec, nrow = ncol(active_response), ncol = ncol(preplay))
    
    temp1 = beta.mat %*% t(preplay) - t(active_response)
    temp.sum1 = sum(diag(t(temp1) %*% V.inv %*% temp1))
    temp2 = beta.vec - beta.init.vec
    temp.sum2 = sum(temp2 * temp2)
    
    return(.5*(temp.sum1 + temp.sum2))
  }
  
  #functions to calculate scale matrix and df of V posterior
  df.V = function(active_response){
    return(ncol(active_response) + nrow(active_response))
  }
  
  scale.V = function(lambda, beta.vec, preplay, active_response, init.phi){
    beta.mat = matrix(beta.vec, nrow = ncol(active_response), ncol = ncol(preplay))
    
    beta.crossprod = tcrossprod(beta.mat %*% t(preplay) - t(active_response))
    
    out = beta.crossprod * (lambda)^-1 + init.phi
    
    return(out)
  }
  
  #functions to calculate mean and variance of W posterior
  Mu.Sigma.W.z = function(active_explanatory, z){
    inv.temp = solve(diag(ncol(active_explanatory)) + crossprod(active_explanatory))
    big.B = inv.temp %*% crossprod(active_explanatory, cbind(z))
    
    hat.matrix = active_explanatory %*% tcrossprod(inv.temp, active_explanatory)
    
    mini.w = diag(hat.matrix)/(1 - diag(hat.matrix))
    
    mini.m = active_explanatory %*% big.B - mini.w*(z - active_explanatory %*% big.B)
    mini.v = mini.w + 1
    
    
    return(list(Mu = big.B, Sigma = inv.temp, mini.m = mini.m, mini.v = mini.v))
  }
  
  #start chain
  for(i in 2:chain.length){
    #calculate V.inv
    V.inv = solve(invvech(V.chain[,i-1]))
    
    #draw lambda
    # ra.lambda = rate.lambda(beta.vec = beta.chain[,i-1], preplay = preplay, active_response = active_response,
    #                          V.inv = V.inv, beta.init.vec = beta.chain[,1])[1]
    # lambda.chain[i] = 1/rgamma(1, shape = shape.lambda(active_response, preplay), rate = ra.lambda)
    lambda.chain[i] = 1
    
    #draw beta
    out.v.star = V.star(preplay_swirl = preplay_swirl, V.inv = V.inv)
    if(min(eigen(out.v.star, symmetric = T, only.values = T)$values) <= 0 ){
      cat('Something went wrong\n The eigen values of out.v.star are', eigen(out.V, symmetric = T, only.values = T)$values)
      break
    }
    v.star.inv = solve(out.v.star)
    m.star = M.star(preplay_swirl = preplay_swirl, active_response = active_response, 
                    V.inv = V.inv, init.beta = beta.chain[,1])
    beta.chain[,i] = c(as.matrix(mvrnorm(1, mu = v.star.inv %*% m.star, Sigma = lambda.chain[i] * v.star.inv)))
    
    #draw V
    sc.V = scale.V(lambda = lambda.chain[i], beta.vec = beta.chain[,i], preplay = preplay, 
                   active_response = active_response, init.phi = init.phi)
    out.V = riwish(v = df.V(active_response = active_response), S = sc.V)
    if(min(eigen(out.V, symmetric = T, only.values = T)$values) <= 0 ){
      cat('Something went wrong\n The eigen values of out.V are', eigen(out.V, symmetric = T, only.values = T)$values)
      break
    }
    V.chain[,i] = vech(out.V)
    
    #temp list of needed variables
    temp.actives = Mu.Sigma.W.z(active_explanatory = active_explanatory, z = z.chain[,i-1])
    
    #draw z
    m.z = temp.actives$mini.m
    sig.z = temp.actives$mini.v
    # z.temp = mvrnorm(n = 1, mu = m.z, Sigma = diag(sig.z))
    z.temp = rnorm(n = length(m.z), mean = m.z, sd = sqrt(sig.z))
    z.chain[,i] = z.temp * 1*(z.temp *(success - .5) > 0)
    
    #draw W
    m.W = temp.actives$Mu
    sig.W = temp.actives$Sigma
    W.chain[,i] = mvrnorm(n = 1, mu = m.W, Sigma = sig.W)
    
    if(!quiet){
      cat('Iteration', i, 'complete.\n')
    }
  }
  
  return(list(beta.chain = beta.chain, V.chain = V.chain, W.chain = W.chain))
}

set.seed(635)
gibb.chains.interaction = full.gibbs(chain.length = 100, data.list = real.out)
saveRDS(gibb.chains, 'gibb.chains.rds')
gibb.chains = readRDS('gibb.chains.rds')
small.chains = list(beta.chain = gibb.chains$beta.chain, W.chain = gibb.chains$W.chain)
saveRDS(small.chains, 'small.chains.rds')
small.chains = readRDS('small.chains.rds')

set.seed(635)
ind.test = sample(length(real.out$success), size = ceiling(.2*length(real.out$success)))
preplay.test = real.out$preplay[ind.test,]
preplay.train = real.out$preplay[-ind.test,]

preplay.explanatory = add_interaction_to_preplay(cbind(1, preplay.test))

data.list = real.out
preplay.explanatory = add_interaction_to_preplay(cbind(1, data.list$preplay))
active.response = data.list$active
active.explanatory = cbind(1, data.list$active)

plot(1000:10000, gibb.chains$beta.chain[1,1000:10000])

beta.vec.out = rowMeans(gibb.chains$beta.chain[,1000:10000])
beta.mat.out = matrix(beta.vec.out, ncol=ncol(preplay.explanatory), nrow=ncol(active.response))
f.norm.response = mean((beta.mat.out %*% t(preplay.explanatory) - t(active.response))^2)
correlation.preplay.active = cor(c(beta.mat.out %*% t(preplay.explanatory)), c(t(active.response)))
plot(beta.mat.out %*% t(preplay.explanatory), t(active_response),main = 'Correlation of True Active Data versus Predicted Active Data', xlab = 'Predicted Active Data', ylab = 'True Active Data', pch=16)
abline(b = 1, a = 0)

quants = apply(gibb.chains$beta.chain[(10*45):(10*45+44),1000:10000], 1, quantile, c(0.025, 0.50, 0.975))
good_ones = which(apply(quants, 2, function(x){return(x[1]*x[3] < 0)}) == F)
quants[,good_ones]


W.vec.out = rowMeans(gibb.chains$W.chain[,1000:10000])
f.norm.success = mean(1*(active.explanatory %*% W.vec.out > 0) == success)

lambda.chain = gibb.chains$lambda.chain
lambda.chain[10:20]
summary(lambda.chain[-1])

posterior.predictive = function(beta.mat, W.vec, data.list){
  preplay = cbind(1, data.list$preplay)
  
  active.response = beta.mat %*% t(preplay)
  
  active.explanatory = cbind(1, t(active.response))
  
  pred_success = 1*(active.explanatory %*%  W.vec > 0)
  
  return(pred_success)
}

pred_success = posterior.predictive(beta.mat.out, W.vec.out, real.out)
mean(pred_success == success)

add_interaction_to_preplay = function(preplay){
  #assign column names
  colnames(preplay) = c('Intercept', 'intended_height', 'intended_top_speed', 'intended_weight', 'close1_height', 'close1_top_speed', 'close1_weight', 'close2_height', 'close2_top_speed', 'close2_weight', 'in_redzone', 'shotgun', 'wr_bunch', 'num_DBs', 'intended_route2','intended_route3','intended_route4','intended_route5','intended_route6','intended_route7','intended_route8','intended_route9','intended_route10','intended_route11', 'close1_route2','close1_route3','close1_route4','close1_route5','close1_route6','close1_route7','close1_route8','close1_route9','close1_route10','close1_route11', 'close2_route2','close2_route3','close2_route4','close2_route5','close2_route6','close2_route7','close2_route8','close2_route9','close2_route10','close2_route11')
  
  preplay = as.data.frame(preplay)
  
  #interaction between top speed of intended receiver and route run by intended receiver
  preplay$'intended_top_speed:intended_route2' = preplay$intended_top_speed * preplay$intended_route2
  preplay$'intended_top_speed:intended_route3' = preplay$intended_top_speed * preplay$intended_route3
  preplay$'intended_top_speed:intended_route4' = preplay$intended_top_speed * preplay$intended_route4
  preplay$'intended_top_speed:intended_route5' = preplay$intended_top_speed * preplay$intended_route5
  preplay$'intended_top_speed:intended_route6' = preplay$intended_top_speed * preplay$intended_route6
  preplay$'intended_top_speed:intended_route7' = preplay$intended_top_speed * preplay$intended_route7
  preplay$'intended_top_speed:intended_route8' = preplay$intended_top_speed * preplay$intended_route8
  preplay$'intended_top_speed:intended_route9' = preplay$intended_top_speed * preplay$intended_route9
  preplay$'intended_top_speed:intended_route10' = preplay$intended_top_speed * preplay$intended_route10
  preplay$'intended_top_speed:intended_route11' = preplay$intended_top_speed * preplay$intended_route11
  
  #interaction between height of intended receiver and route run by intended receiver
  preplay$'intended_height:intended_route2' = preplay$intended_height * preplay$intended_route2
  preplay$'intended_height:intended_route3' = preplay$intended_height * preplay$intended_route3
  preplay$'intended_height:intended_route4' = preplay$intended_height * preplay$intended_route4
  preplay$'intended_height:intended_route5' = preplay$intended_height * preplay$intended_route5
  preplay$'intended_height:intended_route6' = preplay$intended_height * preplay$intended_route6
  preplay$'intended_height:intended_route7' = preplay$intended_height * preplay$intended_route7
  preplay$'intended_height:intended_route8' = preplay$intended_height * preplay$intended_route8
  preplay$'intended_height:intended_route9' = preplay$intended_height * preplay$intended_route9
  preplay$'intended_height:intended_route10' = preplay$intended_height * preplay$intended_route10
  preplay$'intended_height:intended_route11' = preplay$intended_height * preplay$intended_route11
  
  #interaction between weight of intended receiver and route run by intended receiver
  preplay$'intended_weight:intended_route2' = preplay$intended_weight * preplay$intended_route2
  preplay$'intended_weight:intended_route3' = preplay$intended_weight * preplay$intended_route3
  preplay$'intended_weight:intended_route4' = preplay$intended_weight * preplay$intended_route4
  preplay$'intended_weight:intended_route5' = preplay$intended_weight * preplay$intended_route5
  preplay$'intended_weight:intended_route6' = preplay$intended_weight * preplay$intended_route6
  preplay$'intended_weight:intended_route7' = preplay$intended_weight * preplay$intended_route7
  preplay$'intended_weight:intended_route8' = preplay$intended_weight * preplay$intended_route8
  preplay$'intended_weight:intended_route9' = preplay$intended_weight * preplay$intended_route9
  preplay$'intended_weight:intended_route10' = preplay$intended_weight * preplay$intended_route10
  preplay$'intended_weight:intended_route11' = preplay$intended_weight * preplay$intended_route11
  
  # #interaction between intended_routes and close1 routes
  # preplay$'close1_route2:intended_route2' = preplay$close1_route2 * preplay$intended_route2
  # preplay$'close1_route2:intended_route3' = preplay$close1_route2 * preplay$intended_route3
  # preplay$'close1_route2:intended_route4' = preplay$close1_route2 * preplay$intended_route4
  # preplay$'close1_route2:intended_route5' = preplay$close1_route2 * preplay$intended_route5
  # preplay$'close1_route2:intended_route6' = preplay$close1_route2 * preplay$intended_route6
  # preplay$'close1_route2:intended_route7' = preplay$close1_route2 * preplay$intended_route7
  # preplay$'close1_route2:intended_route8' = preplay$close1_route2 * preplay$intended_route8
  # preplay$'close1_route2:intended_route9' = preplay$close1_route2 * preplay$intended_route9
  # preplay$'close1_route2:intended_route10' = preplay$close1_route2 * preplay$intended_route10
  # preplay$'close1_route2:intended_route11' = preplay$close1_route2 * preplay$intended_route11
  # preplay$'close1_route3:intended_route2' = preplay$close1_route3 * preplay$intended_route2
  # preplay$'close1_route3:intended_route3' = preplay$close1_route3 * preplay$intended_route3
  # preplay$'close1_route3:intended_route4' = preplay$close1_route3 * preplay$intended_route4
  # preplay$'close1_route3:intended_route5' = preplay$close1_route3 * preplay$intended_route5
  # preplay$'close1_route3:intended_route6' = preplay$close1_route3 * preplay$intended_route6
  # preplay$'close1_route3:intended_route7' = preplay$close1_route3 * preplay$intended_route7
  # preplay$'close1_route3:intended_route8' = preplay$close1_route3 * preplay$intended_route8
  # preplay$'close1_route3:intended_route9' = preplay$close1_route3 * preplay$intended_route9
  # preplay$'close1_route3:intended_route10' = preplay$close1_route3 * preplay$intended_route10
  # preplay$'close1_route3:intended_route11' = preplay$close1_route3 * preplay$intended_route11
  # preplay$'close1_route4:intended_route2' = preplay$close1_route4 * preplay$intended_route2
  # preplay$'close1_route4:intended_route3' = preplay$close1_route4 * preplay$intended_route3
  # preplay$'close1_route4:intended_route4' = preplay$close1_route4 * preplay$intended_route4
  # preplay$'close1_route4:intended_route5' = preplay$close1_route4 * preplay$intended_route5
  # preplay$'close1_route4:intended_route6' = preplay$close1_route4 * preplay$intended_route6
  # preplay$'close1_route4:intended_route7' = preplay$close1_route4 * preplay$intended_route7
  # preplay$'close1_route4:intended_route8' = preplay$close1_route4 * preplay$intended_route8
  # preplay$'close1_route4:intended_route9' = preplay$close1_route4 * preplay$intended_route9
  # preplay$'close1_route4:intended_route10' = preplay$close1_route4 * preplay$intended_route10
  # preplay$'close1_route4:intended_route11' = preplay$close1_route4 * preplay$intended_route11
  # preplay$'close1_route5:intended_route2' = preplay$close1_route5 * preplay$intended_route2
  # preplay$'close1_route5:intended_route3' = preplay$close1_route5 * preplay$intended_route3
  # preplay$'close1_route5:intended_route4' = preplay$close1_route5 * preplay$intended_route4
  # preplay$'close1_route5:intended_route5' = preplay$close1_route5 * preplay$intended_route5
  # preplay$'close1_route5:intended_route6' = preplay$close1_route5 * preplay$intended_route6
  # preplay$'close1_route5:intended_route7' = preplay$close1_route5 * preplay$intended_route7
  # preplay$'close1_route5:intended_route8' = preplay$close1_route5 * preplay$intended_route8
  # preplay$'close1_route5:intended_route9' = preplay$close1_route5 * preplay$intended_route9
  # preplay$'close1_route5:intended_route10' = preplay$close1_route5 * preplay$intended_route10
  # preplay$'close1_route5:intended_route11' = preplay$close1_route5 * preplay$intended_route11
  # preplay$'close1_route6:intended_route2' = preplay$close1_route6 * preplay$intended_route2
  # preplay$'close1_route6:intended_route3' = preplay$close1_route6 * preplay$intended_route3
  # preplay$'close1_route6:intended_route4' = preplay$close1_route6 * preplay$intended_route4
  # preplay$'close1_route6:intended_route5' = preplay$close1_route6 * preplay$intended_route5
  # preplay$'close1_route6:intended_route6' = preplay$close1_route6 * preplay$intended_route6
  # preplay$'close1_route6:intended_route7' = preplay$close1_route6 * preplay$intended_route7
  # preplay$'close1_route6:intended_route8' = preplay$close1_route6 * preplay$intended_route8
  # preplay$'close1_route6:intended_route9' = preplay$close1_route6 * preplay$intended_route9
  # preplay$'close1_route6:intended_route10' = preplay$close1_route6 * preplay$intended_route10
  # preplay$'close1_route6:intended_route11' = preplay$close1_route6 * preplay$intended_route11
  # preplay$'close1_route7:intended_route2' = preplay$close1_route7 * preplay$intended_route2
  # preplay$'close1_route7:intended_route3' = preplay$close1_route7 * preplay$intended_route3
  # preplay$'close1_route7:intended_route4' = preplay$close1_route7 * preplay$intended_route4
  # preplay$'close1_route7:intended_route5' = preplay$close1_route7 * preplay$intended_route5
  # preplay$'close1_route7:intended_route6' = preplay$close1_route7 * preplay$intended_route6
  # preplay$'close1_route7:intended_route7' = preplay$close1_route7 * preplay$intended_route7
  # preplay$'close1_route7:intended_route8' = preplay$close1_route7 * preplay$intended_route8
  # preplay$'close1_route7:intended_route9' = preplay$close1_route7 * preplay$intended_route9
  # preplay$'close1_route7:intended_route10' = preplay$close1_route7 * preplay$intended_route10
  # preplay$'close1_route8:intended_route2' = preplay$close1_route8 * preplay$intended_route2
  # preplay$'close1_route8:intended_route3' = preplay$close1_route8 * preplay$intended_route3
  # preplay$'close1_route8:intended_route4' = preplay$close1_route8 * preplay$intended_route4
  # preplay$'close1_route8:intended_route5' = preplay$close1_route8 * preplay$intended_route5
  # preplay$'close1_route8:intended_route6' = preplay$close1_route8 * preplay$intended_route6
  # preplay$'close1_route8:intended_route7' = preplay$close1_route8 * preplay$intended_route7
  # preplay$'close1_route8:intended_route8' = preplay$close1_route8 * preplay$intended_route8
  # preplay$'close1_route8:intended_route9' = preplay$close1_route8 * preplay$intended_route9
  # preplay$'close1_route8:intended_route10' = preplay$close1_route8 * preplay$intended_route10
  # preplay$'close1_route8:intended_route11' = preplay$close1_route8 * preplay$intended_route11
  # preplay$'close1_route9:intended_route2' = preplay$close1_route9 * preplay$intended_route2
  # preplay$'close1_route9:intended_route3' = preplay$close1_route9 * preplay$intended_route3
  # preplay$'close1_route9:intended_route4' = preplay$close1_route9 * preplay$intended_route4
  # preplay$'close1_route9:intended_route5' = preplay$close1_route9 * preplay$intended_route5
  # preplay$'close1_route9:intended_route6' = preplay$close1_route9 * preplay$intended_route6
  # preplay$'close1_route9:intended_route7' = preplay$close1_route9 * preplay$intended_route7
  # preplay$'close1_route9:intended_route8' = preplay$close1_route9 * preplay$intended_route8
  # preplay$'close1_route9:intended_route9' = preplay$close1_route9 * preplay$intended_route9
  # preplay$'close1_route9:intended_route10' = preplay$close1_route9 * preplay$intended_route10
  # preplay$'close1_route9:intended_route11' = preplay$close1_route9 * preplay$intended_route11
  # preplay$'close1_route10:intended_route2' = preplay$close1_route10 * preplay$intended_route2
  # preplay$'close1_route10:intended_route3' = preplay$close1_route10 * preplay$intended_route3
  # preplay$'close1_route10:intended_route4' = preplay$close1_route10 * preplay$intended_route4
  # preplay$'close1_route10:intended_route5' = preplay$close1_route10 * preplay$intended_route5
  # preplay$'close1_route10:intended_route6' = preplay$close1_route10 * preplay$intended_route6
  # preplay$'close1_route10:intended_route7' = preplay$close1_route10 * preplay$intended_route7
  # preplay$'close1_route10:intended_route8' = preplay$close1_route10 * preplay$intended_route8
  # preplay$'close1_route10:intended_route9' = preplay$close1_route10 * preplay$intended_route9
  # preplay$'close1_route10:intended_route10' = preplay$close1_route10 * preplay$intended_route10
  # preplay$'close1_route10:intended_route11' = preplay$close1_route10 * preplay$intended_route11
  # preplay$'close1_route11:intended_route2' = preplay$close1_route11 * preplay$intended_route2
  # preplay$'close1_route11:intended_route3' = preplay$close1_route11 * preplay$intended_route3
  # preplay$'close1_route11:intended_route4' = preplay$close1_route11 * preplay$intended_route4
  # preplay$'close1_route11:intended_route5' = preplay$close1_route11 * preplay$intended_route5
  # preplay$'close1_route11:intended_route6' = preplay$close1_route11 * preplay$intended_route6
  # preplay$'close1_route11:intended_route7' = preplay$close1_route11 * preplay$intended_route7
  # preplay$'close1_route11:intended_route8' = preplay$close1_route11 * preplay$intended_route8
  # preplay$'close1_route11:intended_route9' = preplay$close1_route11 * preplay$intended_route9
  # preplay$'close1_route11:intended_route10' = preplay$close1_route11 * preplay$intended_route10
  # preplay$'close1_route11:intended_route11' = preplay$close1_route11 * preplay$intended_route11
  # 
  # #interaction between intended_routes and close2 routes
  # preplay$'close2_route2:intended_route2' = preplay$close2_route2 * preplay$intended_route2
  # preplay$'close2_route2:intended_route3' = preplay$close2_route2 * preplay$intended_route3
  # preplay$'close2_route2:intended_route4' = preplay$close2_route2 * preplay$intended_route4
  # preplay$'close2_route2:intended_route5' = preplay$close2_route2 * preplay$intended_route5
  # preplay$'close2_route2:intended_route6' = preplay$close2_route2 * preplay$intended_route6
  # preplay$'close2_route2:intended_route7' = preplay$close2_route2 * preplay$intended_route7
  # preplay$'close2_route2:intended_route8' = preplay$close2_route2 * preplay$intended_route8
  # preplay$'close2_route2:intended_route9' = preplay$close2_route2 * preplay$intended_route9
  # preplay$'close2_route2:intended_route10' = preplay$close2_route2 * preplay$intended_route10
  # preplay$'close2_route2:intended_route11' = preplay$close2_route2 * preplay$intended_route11
  # preplay$'close2_route3:intended_route2' = preplay$close2_route3 * preplay$intended_route2
  # preplay$'close2_route3:intended_route3' = preplay$close2_route3 * preplay$intended_route3
  # preplay$'close2_route3:intended_route4' = preplay$close2_route3 * preplay$intended_route4
  # preplay$'close2_route3:intended_route5' = preplay$close2_route3 * preplay$intended_route5
  # preplay$'close2_route3:intended_route6' = preplay$close2_route3 * preplay$intended_route6
  # preplay$'close2_route3:intended_route7' = preplay$close2_route3 * preplay$intended_route7
  # preplay$'close2_route3:intended_route8' = preplay$close2_route3 * preplay$intended_route8
  # preplay$'close2_route3:intended_route9' = preplay$close2_route3 * preplay$intended_route9
  # preplay$'close2_route3:intended_route10' = preplay$close2_route3 * preplay$intended_route10
  # preplay$'close2_route3:intended_route11' = preplay$close2_route3 * preplay$intended_route11
  # preplay$'close2_route4:intended_route2' = preplay$close2_route4 * preplay$intended_route2
  # preplay$'close2_route4:intended_route3' = preplay$close2_route4 * preplay$intended_route3
  # preplay$'close2_route4:intended_route4' = preplay$close2_route4 * preplay$intended_route4
  # preplay$'close2_route4:intended_route5' = preplay$close2_route4 * preplay$intended_route5
  # preplay$'close2_route4:intended_route6' = preplay$close2_route4 * preplay$intended_route6
  # preplay$'close2_route4:intended_route7' = preplay$close2_route4 * preplay$intended_route7
  # preplay$'close2_route4:intended_route8' = preplay$close2_route4 * preplay$intended_route8
  # preplay$'close2_route4:intended_route9' = preplay$close2_route4 * preplay$intended_route9
  # preplay$'close2_route4:intended_route10' = preplay$close2_route4 * preplay$intended_route10
  # preplay$'close2_route4:intended_route11' = preplay$close2_route4 * preplay$intended_route11
  # preplay$'close2_route5:intended_route2' = preplay$close2_route5 * preplay$intended_route2
  # preplay$'close2_route5:intended_route3' = preplay$close2_route5 * preplay$intended_route3
  # preplay$'close2_route5:intended_route4' = preplay$close2_route5 * preplay$intended_route4
  # preplay$'close2_route5:intended_route5' = preplay$close2_route5 * preplay$intended_route5
  # preplay$'close2_route5:intended_route6' = preplay$close2_route5 * preplay$intended_route6
  # preplay$'close2_route5:intended_route7' = preplay$close2_route5 * preplay$intended_route7
  # preplay$'close2_route5:intended_route8' = preplay$close2_route5 * preplay$intended_route8
  # preplay$'close2_route5:intended_route9' = preplay$close2_route5 * preplay$intended_route9
  # preplay$'close2_route5:intended_route10' = preplay$close2_route5 * preplay$intended_route10
  # preplay$'close2_route5:intended_route11' = preplay$close2_route5 * preplay$intended_route11
  # preplay$'close2_route6:intended_route2' = preplay$close2_route6 * preplay$intended_route2
  # preplay$'close2_route6:intended_route3' = preplay$close2_route6 * preplay$intended_route3
  # preplay$'close2_route6:intended_route4' = preplay$close2_route6 * preplay$intended_route4
  # preplay$'close2_route6:intended_route5' = preplay$close2_route6 * preplay$intended_route5
  # preplay$'close2_route6:intended_route6' = preplay$close2_route6 * preplay$intended_route6
  # preplay$'close2_route6:intended_route7' = preplay$close2_route6 * preplay$intended_route7
  # preplay$'close2_route6:intended_route8' = preplay$close2_route6 * preplay$intended_route8
  # preplay$'close2_route6:intended_route9' = preplay$close2_route6 * preplay$intended_route9
  # preplay$'close2_route6:intended_route10' = preplay$close2_route6 * preplay$intended_route10
  # preplay$'close2_route6:intended_route11' = preplay$close2_route6 * preplay$intended_route11
  # preplay$'close2_route7:intended_route2' = preplay$close2_route7 * preplay$intended_route2
  # preplay$'close2_route7:intended_route3' = preplay$close2_route7 * preplay$intended_route3
  # preplay$'close2_route7:intended_route4' = preplay$close2_route7 * preplay$intended_route4
  # preplay$'close2_route7:intended_route5' = preplay$close2_route7 * preplay$intended_route5
  # preplay$'close2_route7:intended_route6' = preplay$close2_route7 * preplay$intended_route6
  # preplay$'close2_route7:intended_route7' = preplay$close2_route7 * preplay$intended_route7
  # preplay$'close2_route7:intended_route8' = preplay$close2_route7 * preplay$intended_route8
  # preplay$'close2_route7:intended_route9' = preplay$close2_route7 * preplay$intended_route9
  # preplay$'close2_route7:intended_route10' = preplay$close2_route7 * preplay$intended_route10
  # preplay$'close2_route8:intended_route2' = preplay$close2_route8 * preplay$intended_route2
  # preplay$'close2_route8:intended_route3' = preplay$close2_route8 * preplay$intended_route3
  # preplay$'close2_route8:intended_route4' = preplay$close2_route8 * preplay$intended_route4
  # preplay$'close2_route8:intended_route5' = preplay$close2_route8 * preplay$intended_route5
  # preplay$'close2_route8:intended_route6' = preplay$close2_route8 * preplay$intended_route6
  # preplay$'close2_route8:intended_route7' = preplay$close2_route8 * preplay$intended_route7
  # preplay$'close2_route8:intended_route8' = preplay$close2_route8 * preplay$intended_route8
  # preplay$'close2_route8:intended_route9' = preplay$close2_route8 * preplay$intended_route9
  # preplay$'close2_route8:intended_route10' = preplay$close2_route8 * preplay$intended_route10
  # preplay$'close2_route8:intended_route11' = preplay$close2_route8 * preplay$intended_route11
  # preplay$'close2_route9:intended_route2' = preplay$close2_route9 * preplay$intended_route2
  # preplay$'close2_route9:intended_route3' = preplay$close2_route9 * preplay$intended_route3
  # preplay$'close2_route9:intended_route4' = preplay$close2_route9 * preplay$intended_route4
  # preplay$'close2_route9:intended_route5' = preplay$close2_route9 * preplay$intended_route5
  # preplay$'close2_route9:intended_route6' = preplay$close2_route9 * preplay$intended_route6
  # preplay$'close2_route9:intended_route7' = preplay$close2_route9 * preplay$intended_route7
  # preplay$'close2_route9:intended_route8' = preplay$close2_route9 * preplay$intended_route8
  # preplay$'close2_route9:intended_route9' = preplay$close2_route9 * preplay$intended_route9
  # preplay$'close2_route9:intended_route10' = preplay$close2_route9 * preplay$intended_route10
  # preplay$'close2_route9:intended_route11' = preplay$close2_route9 * preplay$intended_route11
  # preplay$'close2_route10:intended_route2' = preplay$close2_route10 * preplay$intended_route2
  # preplay$'close2_route10:intended_route3' = preplay$close2_route10 * preplay$intended_route3
  # preplay$'close2_route10:intended_route4' = preplay$close2_route10 * preplay$intended_route4
  # preplay$'close2_route10:intended_route5' = preplay$close2_route10 * preplay$intended_route5
  # preplay$'close2_route10:intended_route6' = preplay$close2_route10 * preplay$intended_route6
  # preplay$'close2_route10:intended_route7' = preplay$close2_route10 * preplay$intended_route7
  # preplay$'close2_route10:intended_route8' = preplay$close2_route10 * preplay$intended_route8
  # preplay$'close2_route10:intended_route9' = preplay$close2_route10 * preplay$intended_route9
  # preplay$'close2_route10:intended_route10' = preplay$close2_route10 * preplay$intended_route10
  # preplay$'close2_route10:intended_route11' = preplay$close2_route10 * preplay$intended_route11
  # preplay$'close2_route11:intended_route2' = preplay$close2_route11 * preplay$intended_route2
  # preplay$'close2_route11:intended_route3' = preplay$close2_route11 * preplay$intended_route3
  # preplay$'close2_route11:intended_route4' = preplay$close2_route11 * preplay$intended_route4
  # preplay$'close2_route11:intended_route5' = preplay$close2_route11 * preplay$intended_route5
  # preplay$'close2_route11:intended_route6' = preplay$close2_route11 * preplay$intended_route6
  # preplay$'close2_route11:intended_route7' = preplay$close2_route11 * preplay$intended_route7
  # preplay$'close2_route11:intended_route8' = preplay$close2_route11 * preplay$intended_route8
  # preplay$'close2_route11:intended_route9' = preplay$close2_route11 * preplay$intended_route9
  # preplay$'close2_route11:intended_route10' = preplay$close2_route11 * preplay$intended_route10
  # preplay$'close2_route11:intended_route11' = preplay$close2_route11 * preplay$intended_route11
  
  
  preplay = as.matrix(preplay)
  
  return(preplay)
}

gibbs_no_active = function(chain.length, preplay, success, quiet=FALSE){
  preplay = add_interaction_to_preplay(cbind(1, data.list$preplay))
  success = data.list$success
  
  init.z = rep(0, length(success))
  
  z.chain = matrix(NA, nrow = length(success), ncol = chain.length)
  W.chain = matrix(NA, nrow = ncol(preplay), ncol = chain.length)
  
  z.chain[,1] = init.z
  
  Mu.Sigma.W.z = function(preplay, z){
    inv.temp = solve(diag(ncol(preplay)) + crossprod(preplay))
    big.B = inv.temp %*% crossprod(preplay, cbind(z))
    
    hat.matrix = preplay %*% tcrossprod(inv.temp, preplay)
    
    mini.w = diag(hat.matrix)/(1 - diag(hat.matrix))
    
    mini.m = preplay %*% big.B - mini.w*(z - preplay %*% big.B)
    mini.v = mini.w + 1
    
    return(list(Mu = big.B, Sigma = inv.temp, mini.m = mini.m, mini.v = mini.v))
  }
  
  #start chain
  for(i in 2:chain.length){
    #temp list of needed variables
    temp.actives = Mu.Sigma.W.z(preplay = preplay, z = z.chain[,i-1])
    
    #draw z
    m.z = temp.actives$mini.m
    sig.z = temp.actives$mini.v
    # z.temp = mvrnorm(n = 1, mu = m.z, Sigma = diag(sig.z))
    z.temp = rnorm(n = length(m.z), mean = m.z, sd = sqrt(sig.z))
    z.chain[,i] = z.temp * 1*(z.temp *(success - .5) > 0)
    
    #draw W
    m.W = temp.actives$Mu
    sig.W = temp.actives$Sigma
    W.chain[,i] = mvrnorm(n = 1, mu = m.W, Sigma = sig.W)
    
    if(!quiet){
      cat('Iteration', i, 'complete.\n')
    }
  }
  
  return(list(W.chain = W.chain, z.chain = z.chain))
}

set.seed(635)
ind.test = sample(length(real.out$success), size = ceiling(.2*length(real.out$success)))
preplay.test = real.out$preplay[ind.test,]
preplay.train = real.out$preplay[-ind.test,]
success.test = real.out$success[ind.test]
success.train = real.out$success[-ind.test]

gibbs.without.active = gibbs_no_active(chain.length = 10000, preplay = preplay.train, success = success.train)

preplay.explanatory = add_interaction_to_preplay(cbind(1, preplay.test))
success = success.test

W = rowMeans(gibbs.without.active$W.chain[,-1])
(f.norm.success = mean(1*(preplay.explanatory %*% W > 0) == success))

