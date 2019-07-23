#success of third downs
library(readr)

wd = "C:/Users/Mitch/Documents/UofM/Fall 2018/NFL/Data"
setwd(wd)

all_plays = read_csv(file = 'plays.csv')

#find rows where a third down occurred
third_down_rows = all_plays$down == 3
#find rows where the third down was converted
converted_rows = all_plays$PlayResult >= all_plays$yardsToGo & third_down_rows
#find rows where the third down play was not converted (not just the opposite since I didn't subset)
not_converted_rows = all_plays$PlayResult < all_plays$yardsToGo & third_down_rows

#output visuals
hist(all_plays[not_converted_rows,]$yardsToGo, breaks = 1:max(all_plays[not_converted_rows,]$yardsToGo))

#blue is converted, red not converted
hist(all_plays[c(converted_rows | not_converted_rows),]$yardsToGo, breaks = 0:max(all_plays[c(converted_rows | not_converted_rows),]$yardsToGo), col = 'Red', main = 'Blue is converted 3rd down, Red not converted', xlab = 'yard to gain')
hist(all_plays[converted_rows,]$yardsToGo, breaks = 0:max(all_plays[c(converted_rows | not_converted_rows),]$yardsToGo), col = 'Blue', add = T)
table(all_plays[converted_rows,]$yardsToGo)/(table(all_plays[not_converted_rows & all_plays$yardsToGo <= 21,]$yardsToGo) + table(all_plays[converted_rows,]$yardsToGo))
#odds ratio of 0.75 at 5 so 3/(3+4) chance approximately of making it, makes 4 a good cutoff

