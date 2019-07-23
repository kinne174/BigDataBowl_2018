#basic analysis of stuff to get anything of basic interest
#still not clear how to do regression with this stuff though

#To get an idea of what kind of data I was working with and what I could expect I ran an initial testing on how success was related to each parameter

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

colnames(plays)
# [1] "blitz"                         "cDef1_outcome_direction"       "cDef1_outcome_distance"        "cDef1_outcome_position"       
# [5] "cDef1_outcome_speed"           "cDef1_pass_direction"          "cDef1_pass_distance"           "cDef1_pass_position"          
# [9] "cDef1_pass_speed"              "cDef2_outcome_direction"       "cDef2_outcome_distance"        "cDef2_outcome_position"       
# [13] "cDef2_outcome_speed"           "cDef2_pass_direction"          "cDef2_pass_distance"           "cDef2_pass_position"          
# [17] "cDef2_pass_speed"              "closest_sideline"              "gameId"                        "in_pocket"                    
# [21] "in_redzone"                    "intended_direction"            "intended_distance"             "intended_speed"               
# [25] "num_DBs"                       "num_closest_defenders_to_qb"   "pass_success"                  "playId"                       
# [29] "playerL1_depth"                "playerL1_distance"             "playerL1_position"             "playerL1_route"               
# [33] "playerL2_depth"                "playerL2_distance"             "playerL2_position"             "playerL2_route"               
# [37] "playerL3_depth"                "playerL3_distance"             "playerL3_position"             "playerL3_route"               
# [41] "playerM1_depth"                "playerM1_distance"             "playerM1_position"             "playerM1_route"               
# [45] "playerM2_depth"                "playerM2_distance"             "playerM2_position"             "playerM2_route"               
# [49] "playerR1_depth"                "playerR1_distance"             "playerR1_position"             "playerR1_route"               
# [53] "playerR2_depth"                "playerR2_distance"             "playerR2_position"             "playerR2_route"               
# [57] "playerR3_depth"                "playerR3_distance"             "playerR3_position"             "playerR3_route"               
# [61] "shotgun"                       "success"                       "thrown_away"                   "time_between_pass_and_outcome"
# [65] "time_between_snap_and_pass"    "total_yards"                   "wr_bunch"                      "yards_after_catch"

############################################################################################################################################

#blitz
addmargins(table(plays$blitz, plays$success, dnn=list('blitz', 'success')))
#         success
# blitz   FALSE TRUE  Sum
# FALSE  1489 2665 4154
# TRUE    733  955 1688
# Sum    2222 3620 5842

#looks like there are much more non-blitz plays, but probability of sucess when blitzed (0.57) is less than when not blitzed (0.64)

chisq.test(table(plays$blitz, plays$success, dnn=list('blitz', 'success')))
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  table(plays$blitz, plays$success, dnn = list("blitz", "success"))
# X-squared = 28.934, df = 1, p-value = 7.489e-08

#clearly there is a relationship (see if there still is after a bernoulli correction), blitzing is associated with lower rate of success

############################################################################################################################################

#cDef1_outcome_distance
boxplot(plays$cDef1_outcome_distance ~ plays$success)
#looks pretty similar
summary(plays[plays$success,]$cDef1_outcome_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.050   1.504   2.845   3.393   4.725  24.177      13 
summary(plays[!plays$success,]$cDef1_outcome_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.041   1.564   2.981   3.782   5.191  35.055       8

#pretty similar still not a whole lot of difference here, indicates that defender closeness does not affect sucess

t.test(plays[plays$success,]$cDef1_outcome_distance, plays[!plays$success,]$cDef1_outcome_distance)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef1_outcome_distance and plays[!plays$success, ]$cDef1_outcome_distance
# t = -4.9738, df = 3916.7, p-value = 6.848e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.5428142 -0.2358702
# sample estimates:
#   mean of x mean of y 
# 3.392622  3.781964 

#t test says there is a difference, indicates defenders are closer on successful plays than unsuccessful

############################################################################################################################################

#cDef1_outcome_speed
boxplot(plays$cDef1_outcome_speed ~ plays$success)
#pretty similar, slower defenders for unsucessful plays, seems backwards
summary(plays[plays$success,]$cDef1_outcome_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   3.110   4.540   4.578   5.990  10.140      42
summary(plays[!plays$success,]$cDef1_outcome_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.080   2.820   4.270   4.349   5.780   9.980      21
#agrees with boxplot, slower defenders associated with unscuessful plays

t.test(plays[plays$success,]$cDef1_outcome_speed, plays[!plays$success,]$cDef1_outcome_speed)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef1_outcome_speed and plays[!plays$success, ]$cDef1_outcome_speed
# t = 4.2115, df = 4548.2, p-value = 2.586e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1224595 0.3357682
# sample estimates:
#   mean of x mean of y 
# 4.577705  4.348592 

#higher speed on sucessful plays than unsuccessful, seems backwards

############################################################################################################################################

#cDef1_pass_distance
boxplot(plays$cDef1_pass_distance ~ plays$success)
#successful plays have shorter distance, still seems backwards
summary(plays[plays$success,]$cDef1_pass_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.148   4.702   6.561   7.287   8.845  33.714      13
summary(plays[!plays$success,]$cDef1_pass_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.240   4.973   8.418  10.274  14.148  37.342       8
#still seems backwards

t.test(plays[plays$success,]$cDef1_pass_distance, plays[!plays$success,]$cDef1_pass_distance)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef1_pass_distance and plays[!plays$success, ]$cDef1_pass_distance
# t = -18.49, df = 3162.7, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3.303865 -2.670345
# sample estimates:
#   mean of x mean of y 
# 7.287176 10.274280

#clear that distances are shorter for successful plays than unsucessful plays, this has got to be backwards

############################################################################################################################################

#cDef1_pass_speed
boxplot(plays$cDef1_pass_speed ~ plays$success)
#higher speed for unsuccessful plays, this makes sense
summary(plays[plays$success,]$cDef1_pass_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   2.040   3.300   3.438   4.610   9.310      43 
summary(plays[!plays$success,]$cDef1_pass_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.130   2.460   3.940   4.002   5.450   9.340      2

t.test(plays[plays$success,]$cDef1_pass_speed, plays[!plays$success,]$cDef1_pass_speed)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef1_pass_speed and plays[!plays$success, ]$cDef1_pass_speed
# t = -11.037, df = 4348.3, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.6648769 -0.4642981
# sample estimates:
#   mean of x mean of y 
# 3.437571  4.002158 

#clear difference between speeds for successful and unsucessful plays

############################################################################################################################################

#cDef2_outcome_distance
boxplot(plays$cDef2_outcome_distance ~ plays$success)
#looks pretty similar
summary(plays[plays$success,]$cDef2_outcome_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.535   4.864   7.043   7.444   9.373  30.333      13  
summary(plays[!plays$success,]$cDef2_outcome_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.500   4.628   7.415   8.388  11.108  42.146       8

#still pretty similar, unsuccesful plays have a little more distance

t.test(plays[plays$success,]$cDef2_outcome_distance, plays[!plays$success,]$cDef2_outcome_distance)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef2_outcome_distance and plays[!plays$success, ]$cDef2_outcome_distance
# t = -7.6656, df = 3639.8, p-value = 2.269e-14
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.1846324 -0.7020718
# sample estimates:
#   mean of x mean of y 
# 7.444204  8.387556

#confirmed by t test, backwards from what is expected

############################################################################################################################################

#cDef1_outcome_speed
boxplot(plays$cDef2_outcome_speed ~ plays$success)
#pretty similar
summary(plays[plays$success,]$cDef2_outcome_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   2.740   4.100   4.128   5.430   9.670      43
summary(plays[!plays$success,]$cDef2_outcome_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   2.680   4.130   4.137   5.500   9.470      21
#almost identical

t.test(plays[plays$success,]$cDef2_outcome_speed, plays[!plays$success,]$cDef2_outcome_speed)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef2_outcome_speed and plays[!plays$success, ]$cDef2_outcome_speed
# t = -0.17348, df = 4431.2, p-value = 0.8623
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.10881050  0.09111959
# sample estimates:
#   mean of x mean of y 
# 4.128342  4.137188 

#shows no difference

############################################################################################################################################

#cDef1_pass_distance
boxplot(plays$cDef2_pass_distance ~ plays$success)
#successful plays have shorter distance, still seems backwards
summary(plays[plays$success,]$cDef2_pass_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.485   7.187   9.311  10.402  12.404  41.714      13 
summary(plays[!plays$success,]$cDef2_pass_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.101   8.213  12.832  14.257  19.158  42.981       8
#still seems backwards

t.test(plays[plays$success,]$cDef2_pass_distance, plays[!plays$success,]$cDef2_pass_distance)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef2_pass_distance and plays[!plays$success, ]$cDef2_pass_distance
# t = -21.717, df = 3294.8, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -4.202786 -3.506745
# sample estimates:
#   mean of x mean of y 
# 10.40188  14.25665 

#clear that distances are shorter for successful plays than unsucessful plays, this has got to be backwards
#might make sense for distances since unsuccessful passes will possibly miss receivers which implies defenders will be farther away too
#this doesn't make sense for the speed discrepancies tho

############################################################################################################################################

#cDef1_pass_speed
boxplot(plays$cDef2_pass_speed ~ plays$success)
#higher speed for unsuccessful plays, this makes sense
summary(plays[plays$success,]$cDef2_pass_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   2.120   3.350   3.454   4.570   9.330      43  
summary(plays[!plays$success,]$cDef2_pass_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.010   2.630   3.990   4.006   5.310   8.960      21

t.test(plays[plays$success,]$cDef2_pass_speed, plays[!plays$success,]$cDef2_pass_speed)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$cDef2_pass_speed and plays[!plays$success, ]$cDef2_pass_speed
# t = -11.351, df = 4423.7, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.6472032 -0.4565715
# sample estimates:
#   mean of x mean of y 
# 3.454196  4.006084 

#clear difference between speeds for successful and unsucessful plays

############################################################################################################################################

#closest sideline

boxplot(plays$closest_sideline ~ plays$success)
#pretty clear that unsuccessful plays happen near the sideline more often
summary(plays[plays$success,]$closest_sideline)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   6.795  13.710  13.156  19.060  26.620      13 
summary(plays[!plays$success,]$closest_sideline)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   2.390   7.585   9.701  16.650  26.640       8

t.test(plays[plays$success,]$closest_sideline, plays[!plays$success,]$closest_sideline)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$closest_sideline and plays[!plays$success, ]$closest_sideline
# t = 16.536, df = 4391.7, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   3.045335 3.864572
# sample estimates:
#   mean of x mean of y 
# 13.155523  9.700569 

#pretty clear that the distance to the sideline is impactful of a successful catch

############################################################################################################################################

#in_pocket

addmargins(table(plays$in_pocket, plays$success, dnn=list('in_pocket', 'success')))

#         success
#in_pocketFALSE TRUE  Sum
# FALSE   407  393  800
# TRUE   1815 3227 5042
# Sum    2222 3620 5842

#most passes came from within the pocket but when the QB throws outside of the pocket success is sharply diminished

chisq.test(table(plays$in_pocket, plays$success, dnn=list('in_pocket', 'success')))
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  table(plays$in_pocket, plays$success, dnn = list("in_pocket",     "success"))
# X-squared = 64.212, df = 1, p-value = 1.117e-15

#there is some relationship between passing in the pocket and success

############################################################################################################################################

#num_closest_defenders_to_qb

boxplot(plays$num_closest_defenders_to_qb ~ plays$success)
#boxplot shows that a lower number of defenders leads to sucess

addmargins(table(plays$num_closest_defenders_to_qb, plays$success, dnn=list('num_defenders', 'success')))
#         success
# num_defenders FALSE TRUE  Sum
# 0     781 1848 2629
# 1     892 1123 2015
# 2     411  496  907
# 3     114  127  241
# 4      21   26   47
# 5       2    0    2
# 6       1    0    1
# Sum  2222 3620 5842

#the degredation is clear in the table, the closer the defenders are the better chance at sucess there is for a QB

############################################################################################################################################

#shotgun

addmargins(table(plays$shotgun, plays$success, dnn=list('shotgun', 'success')))
#         success
# shotgun FALSE TRUE  Sum
# FALSE   474  813 1287
# TRUE   1748 2807 4555
# Sum    2222 3620 5842

#A lot of shotgun, surprising but I am only looking at passing plays, seems pretty equal though

chisq.test(table(plays$shotgun, plays$success, dnn=list('shotgun', 'success')))

#yeah no relationship here

############################################################################################################################################

#time_between_pass_and_outcome

boxplot(plays$time_between_pass_and_outcome ~ plays$success)
#shorter pass times correspond with more success
summary(plays[plays$success,]$time_between_pass_and_outcome)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0     8.0    11.0    11.6    14.0    38.0 
summary(plays[!plays$success,]$time_between_pass_and_outcome)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   12.00   15.00   16.44   20.00   53.00 


t.test(plays[plays$success,]$time_between_pass_and_outcome, plays[!plays$success,]$time_between_pass_and_outcome)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$time_between_pass_and_outcome and plays[!plays$success, ]$time_between_pass_and_outcome
# t = -28.881, df = 3455.8, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -5.178545 -4.520136
# sample estimates:
#   mean of x mean of y 
# 11.59530  16.44464 

#t test confirms that shorter passes are better for success

############################################################################################################################################

#time_between_snap_and_pass

boxplot(plays$time_between_snap_and_pass ~ plays$success)
#shorter times qb holding ball correspond with more success
summary(plays[plays$success,]$time_between_snap_and_pass)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3      19      24      25      29      89 
summary(plays[!plays$success,]$time_between_snap_and_pass)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00   22.00   27.00   28.79   33.00   83.00  


t.test(plays[plays$success,]$time_between_snap_and_pass, plays[!plays$success,]$time_between_snap_and_pass)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$time_between_snap_and_pass and plays[!plays$success, ]$time_between_snap_and_pass
# t = -13.849, df = 3918.6, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -4.323101 -3.250889
# sample estimates:
#   mean of x mean of y 
# 25.00193  28.78893 

############################################################################################################################################

#thrown_away

addmargins(table(plays$thrown_away, plays$success, dnn=list('thrown_away', 'success')))
# success
# thrown_away FALSE TRUE  Sum
# FALSE  1862 3588 5450
# TRUE    360   32  392
# Sum    2222 3620 5842

chisq.test(table(plays$thrown_away, plays$success, dnn=list('thrown_away', 'success')))
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  table(plays$thrown_away, plays$success, dnn = list("thrown_away",     "success"))
# X-squared = 513.63, df = 1, p-value < 2.2e-16

############################################################################################################################################

#intended_speed

boxplot(plays$intended_speed ~ plays$success)
#faster speeds correspond with unsuccessful plays
summary(plays[plays$success,]$intended_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   3.900   5.530   5.332   6.850  10.170      17
summary(plays[!plays$success,]$intended_speed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   3.290   4.940   4.758   6.200   9.810      37


t.test(plays[plays$success,]$intended_speed, plays[!plays$success,]$intended_speed)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$intended_speed and plays[!plays$success, ]$intended_speed
# t = -10.335, df = 4445, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.6825222 -0.4648709
# sample estimates:
#   mean of x mean of y 
# 4.757809  5.331506 

#successful plays have slower intended receivers

############################################################################################################################################

#intended_distance

boxplot(plays$intended_distance ~ plays$success)
#successful plays are closer to outcome
summary(plays[plays$success,]$intended_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0100  0.3065  0.4900  0.7687  0.7525 28.9230      13
summary(plays[!plays$success,]$intended_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.050   1.422   2.890   4.089   5.439  46.508       8


t.test(plays[plays$success,]$intended_distance, plays[!plays$success,]$intended_distance)
# Welch Two Sample t-test
# 
# data:  plays[plays$success, ]$intended_distance and plays[!plays$success, ]$intended_distance
# t = -38.907, df = 2582.5, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3.487779 -3.153086
# sample estimates:
#   mean of x mean of y 
# 0.7686873 4.0891197 

#successful plays have intended recievers that are closer

#possible interactions to explore
#routes of the different positions
#blitz/num close defenders and time between snap and pass
#get rid of wr bunch




















