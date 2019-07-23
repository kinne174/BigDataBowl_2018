library(readr)
library(dplyr)

wd = "C:/Users/Mitch/Documents/UofM/Fall 2018/NFL/Data/all_plays"
setwd(wd)
all_files = list.files(path = wd)

all_df = list()
for(i in 1:length(all_files)){
  all_df[[paste(i)]] = z <- read_csv(paste("~/UofM/Fall 2018/NFL/Data/all_plays/",all_files[i], sep=""))
}

plays = bind_rows(all_df)

#I decided to squish the player positions together rather than leave columns with a lot of missing

#reorder columns
plays = plays[,c(1:28, 39:43, 34:38, 29:33, 44:76)]

#remove Nones in positions
remove_None = function(vec){
  for(i in 2:length(vec)){
    if(is.na(vec[i - 1]) & !is.na(vec[i])){
      j = i
      while(is.na(vec[j - 1])){
        #switch and move non-None to the left until its next to another non-None
        vec[j - 1] = vec[j]
        vec[j] = NA
        j = j - 1
        if(j == 1){
          break
        }
      }
    }
  }
  return(vec)
}

route_cols = colnames(plays)[29:68]

#R adopted the tibble which doesn't work the same as a dataframe otherwise I would have just used apply
for(i in 1:nrow(plays)){
  plays[i,route_cols] = remove_None(plays[i, route_cols])
}

#remove columns with all None that aren't being used anymore, as well as any rows that contain None - there were a few plays (84 rows) where the player position player was bad but not horrendous enough to be flagged earlier in the process
plays = plays[,-1*c(54:68)]
plays = plays[complete.cases(plays),]

#took awhile so saving my progress
write_csv(x=plays, path='~/UofM/Fall 2018/NFL/Data/my_data/plays1.csv')

plays = read_csv('~/UofM/Fall 2018/NFL/Data/my_data/plays1.csv')

#filter down position names
pos_filter = function(vec){
  if(vec == "FB"){
    vec = "RB"
  }
  if(vec == 'block.bubble'){
    vec = 'bubble.block'
  }
  if(vec %in% c('CB', 'DB', 'FS', 'SS')){
    vec = 'DB'
  }
  if(vec %in% c('DE', 'DT', 'NT')){
    vec = 'DL'
  }
  if(vec %in% c('ILB', 'LB', 'MLB', 'OLB')){
    vec = 'LB'
  }
  return(vec)
}

for(i in 1:nrow(plays)){
  for(j in 1:length(plays[i,])){
    plays[i,j] = pos_filter(plays[i,j]) 
  }
}

#change RB to DL in defender positions, some oddities got through (5 rows)
plays$cDef1_outcome_position = ifelse(plays$cDef1_outcome_position == 'RB', rep('DL', nrow(plays)), plays$cDef1_outcome_position)
plays$cDef1_pass_position = ifelse(plays$cDef1_pass_position == 'RB', rep('DL', nrow(plays)), plays$cDef1_pass_position)
plays$cDef2_outcome_position = ifelse(plays$cDef2_outcome_position == 'RB', rep('DL', nrow(plays)), plays$cDef2_outcome_position)
plays$cDef2_pass_position = ifelse(plays$cDef2_pass_position == 'RB', rep('DL', nrow(plays)), plays$cDef2_pass_position)

#took awhile so saving my progress
write_csv(x=plays, path='~/UofM/Fall 2018/NFL/Data/my_data/plays1.csv')

plays = read_csv('~/UofM/Fall 2018/NFL/Data/my_data/plays1.csv')

#finding which columns correspond with route information and all other rows will go into logistic regression
routes_intended_cols = colnames(plays)[c(19, 27:28, 31, 33, 36, 38, 41, 43, 46, 48, 51, 53)]

all_else_cols = colnames(plays)[c(1:18, 20, 22:27, 29, 32, 34, 37, 39, 42, 44, 47, 49, 52, 54, 57:58)]

#set my seed for reproducability
set.seed(1001)
model_rows = sample(nrow(plays), ceiling(nrow(plays)/3)) #1:2 split

routes_df = plays[-model_rows, routes_intended_cols]
model_routes_df = plays[-model_rows, all_else_cols] #also need to save the columns for routes to test the info in the logistic regression
model_df = plays[model_rows, all_else_cols]

write_csv(routes_df, '~/UofM/Fall 2018/NFL/Data/my_data/routes_df.csv')
write_csv(model_routes_df, '~/UofM/Fall 2018/NFL/Data/my_data/model_routes_df.csv')
write_csv(model_df, '~/UofM/Fall 2018/NFL/Data/my_data/model_df.csv')








