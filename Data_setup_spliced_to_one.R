# this is spliced all together to start at the same timepoint
library(data.table)
library(plyr)
library(ggplot2)
# Data setup, splicing seasons together

PHIRST_NP <- data.table(read.csv("~/Documents/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Nasopharyngeal specimens 2016-2018 Flu & RSV 2022-04-19.csv",
                                 na.strings = ""))
# format correctly
main_db <- PHIRST_NP[,c("indid", "npsdatecol", "flu", "rsv")]
main_db$indid <- as.factor(main_db$indid)
main_db$npsdatecol <- as.Date(main_db$npsdatecol, "%d-%m-%Y")
# remove the rows that have NA obsevations for flu/rsv (i.e. no data recoreded)
main_db <- na.omit(main_db)
# earliest date recprded
start_date <- min(main_db$npsdatecol)
# convert to numeric time
main_db[,timestep := (npsdatecol - start_date)]
main_db[,timestep := as.numeric(timestep)]

##### Investigate the numbers that get infected multiple times ####
num_days_neg <- 14
flu_subset <- main_db[flu == "Pos"]
num_flu_infected <- length(unique(flu_subset$indid))
flu_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =indid]
flu_tab <- table(table(flu_subset[time_diff > num_days_neg]$indid))
print(sum(flu_tab[-1])/num_flu_infected*100)

### how many double positives for flu, assuming X days negative in between)
rsv_subset <- main_db[rsv == "Pos"]
num_rsv_infected <- length(unique(rsv_subset$indid))
rsv_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =indid]
rsv_tab <- table(table(rsv_subset[time_diff > num_days_neg]$indid))
print(sum(rsv_tab[-1])/num_rsv_infected*100)

### respecify as positive those that are within 14 days of previous, 
### and remove the ones with subsequent re-infections
repeat_rsvs <- rsv_subset[time_diff > num_days_neg]
repeat_flus <- flu_subset[time_diff > num_days_neg]
main_db$remove_repeat <- F
# mark the flus to remove
for(i in 1:nrow(repeat_flus)){
  main_db[indid == unlist(repeat_flus[i,"indid"]) &
            timestep >= unlist(repeat_flus[i,"timestep"]) , remove_repeat := T]
}
# mark the rsvs to remove
for(i in 1:nrow(repeat_rsvs)){
  main_db[indid == unlist(repeat_rsvs[i,"indid"]) &
            timestep >= unlist(repeat_rsvs[i,"timestep"]) , remove_repeat := T]
  
}

main_db <- main_db[remove_repeat == F ]
main_db[, factor_id := as.numeric(factor(indid))]

#recode to same format
main_db[flu=="Neg", flu := 0]
main_db[flu=="Pos", flu := 1]
main_db[rsv=="Neg", rsv := 0]
main_db[rsv=="Pos", rsv := 1]

# convert them to all be the same year
main_db[, combi_step := substring(npsdatecol,6)]
main_db[, combi_step := paste0("2015-",combi_step)]
min_artificial_date <- min(main_db$combi_step)
main_db[,timestep := (as.Date(combi_step) - as.Date(min_artificial_date))]
main_db[,timestep := as.numeric(timestep)]

# remove unneccesary columns
main_db <- main_db[,c("factor_id","rsv", "flu", "timestep")]
main_db$state <- "tbd"



for(i in unique(main_db$factor_id)){
  # subset for each person
  temp <- main_db[factor_id == i,]
  # calculate timings of infection
  flu_time <- as.numeric(unlist(temp[flu == 1, "timestep"]))
  rsv_time <- as.numeric(unlist(temp[rsv == 1, "timestep"]))
  if(is.na(rsv_time[1])){rsv_time = Inf}
  if(is.na(flu_time[1])){flu_time = Inf}
  # specify the state at each timepoint
  temp[timestep < flu_time[1] & timestep < rsv_time[1], state := "SS"]
  temp[timestep > flu_time[1] & timestep < rsv_time[1], state := "SPR"]
  temp[timestep < flu_time[1] & timestep >rsv_time[1], state := "PRS"]
  temp[timestep > flu_time[1] & timestep > rsv_time[1], state := "RR"]
  temp[timestep == flu_time[1] & timestep< rsv_time[1], state := "SI"]
  temp[timestep < flu_time[1] & timestep == rsv_time[1], state := "IS"]
  temp[timestep > flu_time[1] & timestep == rsv_time[1], state := "IPR"]
  temp[timestep == flu_time[1] & timestep > rsv_time[1], state := "PRI"]
  temp[timestep == flu_time[1] & timestep == rsv_time[1], state := "II"]
  temp[timestep > tail(flu_time,1) & timestep > tail(rsv_time,1), state := "RR"]
  
  main_db[factor_id == i,] <- temp
}
# recode to numeric
main_db[state == "SS", num_state :=1]
main_db[state == "IS", num_state :=2]
main_db[state == "PRS", num_state :=99]
main_db[state == "SI", num_state :=5]
main_db[state == "II", num_state :=6]
main_db[state == "PRI", num_state :=7]
main_db[state == "SPR", num_state :=101]
main_db[state == "IPR", num_state :=9]
main_db[state == "RR", num_state :=11]

#save(main_db, file="main_db_together_april.Rdata")
library(plyr)
main_db[, rounded := round_any(timestep, 10, f = ceiling)]
tester <- main_db[,.N, by = c("rounded", "state")]

ggplot(tester, aes(x = rounded, y = N, colour = state)) + geom_line()
