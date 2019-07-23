#logistic regression on success

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

na2none = function(vec){
  for(i in 1:length(vec)){
    if(is.na(vec[i])){
      vec[i] = 'None'
    }
  }
  return(vec)
}

na2zero = function(vec){
  for(i in 1:length(vec)){
    if(is.na(vec[i])){
      vec[i] = 0
    }
  }
  return(vec)
}

none_cols = c()
zero_cols = c()
for(c_name in colnames(plays)){
  if(nchar(c_name) < 14){
    next
  }
  if(substr(c_name, 10, 14)  == 'route'){
    none_cols = append(none_cols, c_name)
  }
  if(substr(c_name, 10, 17) == 'position'){
    none_cols = append(none_cols, c_name)
  }
  
  if(substr(c_name, 10, 14) == 'depth'){
    zero_cols = append(zero_cols, c_name)
  }
  if(substr(c_name, 10, 17) == 'distance'){
    zero_cols = append(zero_cols, c_name)
  }
}

zero_cols = zero_cols[-1]

plays[,zero_cols] = apply(plays[,zero_cols], 2, na2zero)
plays[,none_cols] = apply(plays[,none_cols], 2, na2none)

plays = plays[complete.cases(plays),]

#change block.bubble
route_columns = c("playerL3_route", "playerL2_route", "playerL1_route", "playerM2_route", "playerM1_route", "playerR1_route", "playerR2_route", "playerR3_route")

block2bubble = function(vec){
  for(i in 1:length(vec)){
    if(vec[i] == 'block.bubble'){
      vec[i] = 'bubble.block'
    }
  }
  return(vec)
}

plays[,route_columns] = apply(plays[,route_columns], 2, block2bubble)

#classify defense
classify_def_positions = function(vec){
  for(i in 1:length(vec)){
    if(vec[i] %in% c('CB', 'DB', 'FS', 'SS')){
      vec[i] = 'DB'
    }
    if(vec[i] %in% c('DE', 'DT', 'NT')){
      vec[i] = 'DL'
    }
    if(vec[i] %in% c('ILB', 'LB', 'MLB', 'OLB')){
      vec[i] = 'LB'
    }
  }
  return(vec)
}
plays[,c("cDef1_outcome_position", "cDef2_outcome_position", "cDef1_pass_position", "cDef2_pass_position")] = apply(plays[,c("cDef1_outcome_position", "cDef2_outcome_position", "cDef1_pass_position", "cDef2_pass_position")], 2, classify_def_positions)

cols_to_use = c("blitz", "cDef1_outcome_distance", "cDef1_outcome_position", "cDef1_outcome_speed", "cDef1_pass_distance", "cDef1_pass_position", "cDef1_pass_speed", "cDef2_outcome_distance", "cDef2_outcome_position", "cDef2_outcome_speed", "cDef2_pass_distance", "cDef2_pass_position", "cDef2_pass_speed", "closest_sideline", "in_pocket", "intended_direction", "intended_distance", "intended_speed", "num_DBs", "num_closest_defenders_to_qb", "playerL1_depth", "playerL1_position", "playerL1_route", "playerL2_depth", "playerL2_position", "playerL2_route", "playerL3_depth", "playerL3_position", "playerL3_route", "playerM1_depth", "playerM1_position", "playerM1_route", "playerM2_depth", "playerM2_position", "playerM2_route", "playerR1_depth", "playerR1_position", "playerR1_route", "playerR2_depth", "playerR2_position", "playerR2_route", "playerR3_depth", "playerR3_position", "playerR3_route", "shotgun", "pass_success", "time_between_pass_and_outcome", "time_between_snap_and_pass")

#changing reference level to None
plays = within(plays, {
  playerL3_route <- relevel(as.factor(playerL3_route), ref = 'None')
  playerL2_route <- relevel(as.factor(playerL2_route), ref = 'None')
  playerL1_route <- relevel(as.factor(playerL1_route), ref = 'None')
  playerM1_route <- relevel(as.factor(playerM1_route), ref = 'None')
  playerM2_route <- relevel(as.factor(playerM2_route), ref = 'None')
  playerR1_route <- relevel(as.factor(playerR1_route), ref = 'None')
  playerR2_route <- relevel(as.factor(playerR2_route), ref = 'None')
  playerR3_route <- relevel(as.factor(playerR3_route), ref = 'None')
  playerL3_position <- relevel(as.factor(playerL3_position), ref = 'None')
  playerL2_position <- relevel(as.factor(playerL2_position), ref = 'None')
  playerL1_position <- relevel(as.factor(playerL1_position), ref = 'None')
  playerM1_position <- relevel(as.factor(playerM1_position), ref = 'None')
  playerM2_position <- relevel(as.factor(playerM2_position), ref = 'None')
  playerR1_position <- relevel(as.factor(playerR1_position), ref = 'None')
  playerR2_position <- relevel(as.factor(playerR2_position), ref = 'None')
  playerR3_position <- relevel(as.factor(playerR3_position), ref = 'None')
  
  })

options(max.print = 1500)
lr = glm(pass_success ~ ., data = plays[,cols_to_use], family = 'binomial')
summary(lr)

table(plays$success, plays$playerM1_route)

plays_copy = plays[,cols_to_use]
plays_copy$playerL1_route = ifelse((plays$playerL1_route == 'bubble.block') & (plays$playerL1_position == 'FB'), rep('None', nrow(plays)), as.character(plays$playerL1_route))
plays_copy$playerL2_route = ifelse((plays$playerL2_route == 'swing') | (plays$playerL2_route == 'v.in') | (plays$playerL2_route == 'v.out'), rep('None', nrow(plays)), as.character(plays$playerL2_route))
plays_copy$playerL3_route = ifelse((plays$playerL3_route == 'swing') | (plays$playerL3_route == 'v.in') | (plays$playerL3_route == 'v.out') | (plays$playerL3_route == 'stick.out') | (plays$playerL3_route == 'flat'), rep('None', nrow(plays)), as.character(plays$playerL3_route))
plays_copy$playerR3_route = ifelse((plays$playerR3_route == 'swing') | (plays$playerR3_route == 'v.in') | (plays$playerR3_route == 'v.out') | (plays$playerR3_route == 'stick.out') | (plays$playerR3_route == 'flat'), rep('None', nrow(plays)), as.character(plays$playerR3_route))
plays_copy$playerR2_route = ifelse((plays$playerR2_route == 'swing') | (plays$playerR2_route == 'v.in') | (plays$playerR2_route == 'v.out'), rep('None', nrow(plays)), as.character(plays$playerR2_route))
plays_copy$playerR1_route = ifelse((plays$playerR1_route == 'v.in'), rep('None', nrow(plays)), as.character(plays$playerR1_route))
plays_copy$playerM2_route = ifelse((plays$playerM2_route == 'bubble.block') & (plays$playerM2_position == 'FB'), rep('None', nrow(plays)), as.character(plays$playerM2_route))
plays_copy$playerM1_route = ifelse((plays$playerM1_route == 'bubble.block') & (plays$playerM1_position == 'FB'), rep('None', nrow(plays)), as.character(plays$playerM1_route))

plays_copy = within(plays_copy, {
  playerL3_route <- relevel(as.factor(playerL3_route), ref = 'None')
  playerL2_route <- relevel(as.factor(playerL2_route), ref = 'None')
  playerL1_route <- relevel(as.factor(playerL1_route), ref = 'None')
  playerM1_route <- relevel(as.factor(playerM1_route), ref = 'None')
  playerM2_route <- relevel(as.factor(playerM2_route), ref = 'None')
  playerR1_route <- relevel(as.factor(playerR1_route), ref = 'None')
  playerR2_route <- relevel(as.factor(playerR2_route), ref = 'None')
  playerR3_route <- relevel(as.factor(playerR3_route), ref = 'None')
})

lr_copy = glm(success ~., data = plays_copy, family = 'binomial')
summary(lr_copy)

cols_without_routes = c("blitz", "cDef1_outcome_distance", "cDef1_outcome_position", "cDef1_outcome_speed", "cDef1_pass_distance", "cDef1_pass_position", "cDef1_pass_speed", "cDef2_outcome_distance", "cDef2_outcome_position", "cDef2_outcome_speed", "cDef2_pass_distance", "cDef2_pass_position", "cDef2_pass_speed", "closest_sideline", "in_pocket", "intended_direction", "intended_distance", "intended_speed", "num_DBs", "num_closest_defenders_to_qb", "shotgun", "success", "time_between_pass_and_outcome", "time_between_snap_and_pass")

col_routes_only = c("playerL1_depth", "playerL1_position", "playerL1_route", "playerL2_depth", "playerL2_position", "playerL2_route", "playerL3_depth", "playerL3_position", "playerL3_route", "playerM1_depth", "playerM1_position", "playerM1_route", "playerM2_depth", "playerM2_position", "playerM2_route", "playerR1_depth", "playerR1_position", "playerR1_route", "playerR2_depth", "playerR2_position", "playerR2_route", "playerR3_depth", "playerR3_position", "playerR3_route", "success")

lr = glm(success ~., data = plays_copy, family = 'binomial')
summary(lr)
