library(readr)
library(e1071)

wd = "C:/Users/Mitch/Documents/UofM/Fall 2018/NFL/Data/my_data"
setwd(wd)

#build logistic regression model

model_df = read_csv('model_df.csv')

#get a training and testing model
set.seed(1001)
test_ind = sample(nrow(model_df), ceiling(nrow(model_df)/5))
test_df = model_df[test_ind,]
train_df = model_df[-test_ind,]


simple_train_df = model_df[,c(2:3, 5:9, 16:21, 23, 25, 36)]

lr_simple = glm(pass_success~.^2, data = simple_train_df, family = 'binomial')
lr_full = glm(pass_success~. + cDef1_outcome_direction*cDef1_outcome_position + cDef1_outcome_distance*cDef1_outcome_position + cDef1_outcome_speed*cDef1_outcome_position + cDef2_outcome_direction*cDef2_outcome_position + cDef2_outcome_distance*cDef2_outcome_position + cDef2_outcome_speed*cDef2_outcome_position + cDef1_pass_direction*cDef1_pass_position + cDef1_pass_distance*cDef1_pass_position + cDef1_pass_speed*cDef1_pass_position + cDef2_pass_direction*cDef2_pass_position + cDef2_pass_distance*cDef2_pass_position + cDef2_pass_speed*cDef2_pass_position + playerL3_depth*playerL3_position + playerL2_depth*playerL2_position + playerL1_depth*playerL1_position + playerM2_depth*playerM2_position + playerM1_depth*playerM1_position, data = train_df, family='binomial')
summary(lr)

step_model_simple = step(lr_simple, direction='backward', trace = 0)
summary(step_model_simple)

#using all the data without interaction or squared effects perform a backwards step process using AIC
lr = glm(pass_success~., data = train_df, family = 'binomial')
step_model = step(lr, direction = 'backward', trace = 0)
step_model = update(step_model, ~.-playerL3_position - playerL3_depth - playerL1_depth)
summary(step_model)

#get ~73% correct prediction rate
sum(test_df$pass_success == (predict(step_model_simple, newdata = test_df[,-25], type = 'response') >= 0.5))/nrow(test_df)

#########################################################################################################################################

#see what routes are successful by careful tabulation

routes_df = read_csv('~/UofM/Fall 2018/NFL/Data/my_data/routes_df.csv')

#L3 is the inteded receiver
table(routes_df[routes_df$playerL3_intended,]$pass_success, routes_df[routes_df$playerL3_intended,]$playerL3_route)[2,]/(table(routes_df[routes_df$playerL3_intended,]$pass_success, routes_df[routes_df$playerL3_intended,]$playerL3_route)[2,] + table(routes_df[routes_df$playerL3_intended,]$pass_success, routes_df[routes_df$playerL3_intended,]$playerL3_route)[1,]) 
table(routes_df[routes_df$playerL3_intended,]$pass_success, routes_df[routes_df$playerL3_intended,]$playerL3_route) #best in order are dig, block, out, post, flat, comeback

#using dig as intended route for L3 what is good L2/L1 routes
table(routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'dig',]$pass_success, routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'dig',]$playerL2_route) #flat or out
table(routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'dig',]$pass_success, routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'dig',]$playerL1_route) #slant or block
#get indices of plays using successful route combinations
L3_dig_L1_slant_ind = which(routes_df$playerL3_intended & routes_df$playerL3_route == 'dig' & routes_df$playerL1_route == 'slant')


#using blcok as intended route for L3 what is good L2/L1 routes
table(routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'bubble.block',]$pass_success, routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'bubble.block',]$playerL2_route) #block
table(routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'bubble.block',]$pass_success, routes_df[routes_df$playerL3_intended & routes_df$playerL3_route == 'bubble.block',]$playerL1_route) #flat


#L2 is the inteded receiver
table(routes_df[routes_df$playerL2_intended,]$pass_success, routes_df[routes_df$playerL2_intended,]$playerL2_route)[2,]/(table(routes_df[routes_df$playerL2_intended,]$pass_success, routes_df[routes_df$playerL2_intended,]$playerL2_route)[2,] + table(routes_df[routes_df$playerL2_intended,]$pass_success, routes_df[routes_df$playerL2_intended,]$playerL2_route)[1,]) 
table(routes_df[routes_df$playerL2_intended,]$pass_success, routes_df[routes_df$playerL2_intended,]$playerL2_route) #best in order are post, dig, corner, flat

#Using post as intended L2 route what are good L3/L1 routes
table(routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'post',]$pass_success, routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'post',]$playerL3_route) #out
table(routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'post',]$pass_success, routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'post',]$playerL1_route) #flat
#get indices of plays using successful route combinations
L2_post_L1_flat_ind = which(routes_df$playerL2_intended & routes_df$playerL2_route == 'post' & routes_df$playerL1_route == 'flat')

#Using dig as intended L2 route what are good L3/L1 routes
table(routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'dig',]$pass_success, routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'dig',]$playerL3_route) #post or wheel
table(routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'dig',]$pass_success, routes_df[routes_df$playerL2_intended & routes_df$playerL2_route == 'dig',]$playerL1_route) #post or flat


#Using L1 as intended receiever
table(routes_df[routes_df$playerL1_intended,]$pass_success, routes_df[routes_df$playerL1_intended,]$playerL1_route)[2,]/(table(routes_df[routes_df$playerL1_intended,]$pass_success, routes_df[routes_df$playerL1_intended,]$playerL1_route)[2,] + table(routes_df[routes_df$playerL1_intended,]$pass_success, routes_df[routes_df$playerL1_intended,]$playerL1_route)[1,]) 
table(routes_df[routes_df$playerL1_intended,]$pass_success, routes_df[routes_df$playerL1_intended,]$playerL1_route) #best in order are corner, curl, out

#Using corner as intended L1 route what are good L2/M1 routes
table(routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'corner',]$pass_success, routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'corner',]$playerL2_route) #flat
#get indices of plays using successful route combinations
L1_corner_L2_flat_ind = which(routes_df$playerL1_route == 'corner' & routes_df$playerL1_intended & routes_df$playerL2_route == 'flat')
table(routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'corner',]$pass_success, routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'corner',]$playerM1_route) #flat

#Using curl as intended L1 route what are good L2/M1 routes
table(routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'curl',]$pass_success, routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'curl',]$playerL2_route) #nothing
table(routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'curl',]$pass_success, routes_df[routes_df$playerL1_intended & routes_df$playerL1_route == 'curl',]$playerM1_route) #nothing

#Using M1 as intended receiver
table(routes_df[routes_df$playerM1_intended,]$pass_success, routes_df[routes_df$playerM1_intended,]$playerM1_route)[2,]/(table(routes_df[routes_df$playerM1_intended,]$pass_success, routes_df[routes_df$playerM1_intended,]$playerM1_route)[2,] + table(routes_df[routes_df$playerM1_intended,]$pass_success, routes_df[routes_df$playerM1_intended,]$playerM1_route)[1,]) 
table(routes_df[routes_df$playerM1_intended,]$pass_success, routes_df[routes_df$playerM1_intended,]$playerM1_route) #best in order are out, dig

#Using out as intended M1 route what are good L1/M2 routes
table(routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'out',]$pass_success, routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'out',]$playerL1_route) #flat
#get indices of plays using successful route combinations
M1_out_L1_flat_indices = which(routes_df$playerM1_route == 'out' & routes_df$playerL1_route == 'flat' & routes_df$playerM1_intended)
table(routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'out',]$pass_success, routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'out',]$playerM2_route) #nothing

#Using dig as intended M1 route what are good L1/M2 routes
table(routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'dig',]$pass_success, routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'dig',]$playerL1_route) #flat or slant
table(routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'dig',]$pass_success, routes_df[routes_df$playerM1_intended & routes_df$playerM1_route == 'dig',]$playerM2_route) #dig
#get indices of plays using successful route combinations
M1_dig_M2_dig_inds = which(routes_df$playerM1_intended & routes_df$playerM1_route == 'dig' & routes_df$playerM2_route == 'dig')

#Using M2 as intedend receiver
table(routes_df[routes_df$playerM2_intended,]$pass_success, routes_df[routes_df$playerM2_intended,]$playerM2_route)[2,]/(table(routes_df[routes_df$playerM2_intended,]$pass_success, routes_df[routes_df$playerM2_intended,]$playerM2_route)[2,] + table(routes_df[routes_df$playerM2_intended,]$pass_success, routes_df[routes_df$playerM2_intended,]$playerM2_route)[1,]) 
table(routes_df[routes_df$playerM2_intended,]$pass_success, routes_df[routes_df$playerM2_intended,]$playerM2_route) #best in order are post, dig, flat

#Using post as intended M2 route what are good L1/M1 routes
table(routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'post',]$pass_success, routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'post',]$playerL1_route) #block
table(routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'post',]$pass_success, routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'post',]$playerM1_route) #slant
#get indices of plays using successful route combinations
M2_post_M1_slant_inds = which(routes_df$playerM2_intended & routes_df$playerM2_route == 'post' & routes_df$playerM1_route == 'slant')

#Using dig as intended M2 route what are good L1/M1 routes
table(routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'dig',]$pass_success, routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'dig',]$playerL1_route) #block or flat
table(routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'dig',]$pass_success, routes_df[routes_df$playerM2_intended & routes_df$playerM2_route == 'dig',]$playerM1_route) #nothing

#get model attributes of the routes_df, see if they are predicted to do well too

model_routes_df = read_csv('~/UofM/Fall 2018/NFL/Data/my_data/model_routes_df.csv')
sum((predict(step_model, newdata = model_routes_df[M1_out_L1_flat_indices, -25], type = 'response') >= 0.5) == model_routes_df[M1_out_L1_flat_indices, 25])/nrow(model_routes_df[M1_out_L1_flat_indices,])

sum((predict(step_model, newdata = model_routes_df[M1_out_L1_flat_indices, -25], type = 'response') >= 0.5))
sum(model_routes_df[M1_out_L1_flat_indices, 25, drop = T])

#use draw play to see if there is any consistent pattern on how the routes were run
routes_df[M1_out_L1_flat_indices, c('gameId', 'playId')] #good fail: 1705; 1255, good successes: 091100;2504

sum((predict(step_model, newdata = model_routes_df[L2_post_L1_flat_ind, -25], type = 'response') >= 0.5) == model_routes_df[L2_post_L1_flat_ind, 25])/nrow(model_routes_df[L2_post_L1_flat_ind,])

#use draw play to see if there is any consistent pattern on how the routes were run
routes_df[L2_post_L1_flat_ind, c('gameId', 'playId')] #good success: 092406;1802

sum((predict(step_model, newdata = model_routes_df[M2_post_M1_slant_inds, -25], type = 'response') >= 0.5) == model_routes_df[M2_post_M1_slant_inds, 25])/nrow(model_routes_df[M2_post_M1_slant_inds,])

#use draw play to see if there is any consistent pattern on how the routes were run
routes_df[M2_post_M1_slant_inds, c('gameId', 'playId')] #good: 092409;2600/2622, 101507;2800

sum((predict(step_model, newdata = model_routes_df[M1_dig_M2_dig_inds, -25], type = 'response') >= 0.5) == model_routes_df[M1_dig_M2_dig_inds, 25])/nrow(model_routes_df[M1_dig_M2_dig_inds,])

#use draw play to see if there is any consistent pattern on how the routes were run
routes_df[M1_dig_M2_dig_inds, c('gameId', 'playId', 'pass_success')] #good: 1, 4

predict(step_model, new_data=model_routes_df[,-25], type = 'response')
which(routes_df$gameId == 2017092409 & routes_df$playId == 2600)

#started working on seeing if could find plays that had good qualities using quantiles within the columns themselves, did not proceed far with this idea though

# percentile = function(vec){
#   out_vec = c()
#   perc = ecdf(vec)
#   for(i in 1:length(vec)){
#     out_vec[i] = perc(vec[i])
#   }
#   return(out_vec)
# }
# 
# 
# outcome_direction_perc = percentile(model_routes_df$cDef1_outcome_direction)
# outcome_distance_perc = percentile(model_routes_df$cDef1_outcome_distance)
# outcome_speed_perc = percentile(model_routes_df$cDef1_outcome_speed)
# pass_direction_perc = percentile(model_routes_df$cDef1_pass_direction)
# pass_distance_perc = percentile(model_routes_df$cDef1_pass_distance)
# closest_sideline_perc = percentile(model_routes_df$closest_sideline)
# intended_distance_perc = percentile(model_routes_df$intended_distance)
# intended_direction_perc = percentile(model_routes_df$intended_direction)
# percentile_df = data.frame(outcome_direction_perc, outcome_distance_perc, outcome_speed_perc, pass_direction_perc, pass_distance_perc, closest_sideline_perc, intended_distance_perc)
# 
# good_pred_inds = which(percentile_df$outcome_distance_perc > .5 & percentile_df$outcome_speed_perc > .5 & percentile_df$pass_distance_perc < .5 & percentile_df$closest_sideline_perc > .5 & percentile_df$intended_distance_perc < .5 & percentile_df$pass_direction_perc > .5)
# 
# for(i in 1:length(good_pred_inds)){
#   print(routes_df[good_pred_inds[i], c('playerL3_route', 'playerL2_route', 'playerL1_route', 'playerM1_route', 'playerM2_route'), drop = T])
# }

