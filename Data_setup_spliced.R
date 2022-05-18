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

# remove unneccesary columns
main_db <- main_db[,c("npsdatecol","factor_id","rsv", "flu", "timestep")]
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

#### splice the timeslots together

# times when there are none: 182-258, 545 - 622, 912+
# check they are empty
main_db[timestep  %in% c(182:258)]
main_db[timestep  %in% c(545:622)]
main_db[timestep  %in% c(913:1200)]
main_db$track_colour <- "first"
main_db[timestep > 258,track_colour := "second"]
main_db[timestep > 258, timestep := timestep-77]
main_db[timestep >469, track_colour := "third"]
main_db[timestep > 469, timestep := timestep-78]

#save(main_db, file="main_db_seasons_april.Rdata")

# total numbers of each by timestep
main_db$flu <- as.numeric(as.character(main_db$flu))
main_db$rsv <- as.numeric(as.character(main_db$rsv))
main_db$timestep <- as.numeric(main_db$timestep)
to_plot <- aggregate(main_db[, c("flu", "rsv")], by=list(main_db$timestep), sum)
# use this if want to look at daily reports
to_plot_m1 <- melt(to_plot, id.vars = "Group.1")
ggplot(to_plot_m1, aes(x = Group.1, y = value, colour = variable)) +
  geom_line()
# calulate which month of the year it is (i.e. seasonality)
to_plot$month <- month(as.Date(to_plot$Group.1, origin = start_date))
to_plot <- aggregate(to_plot[, c("flu", "rsv")], by = list(to_plot$month), sum)
to_plot_m2 <- melt(to_plot, id.vars = "Group.1")
# plot the month of year vs cases (i.e. seasonality)
ggplot(to_plot_m2, aes(x = Group.1, y = value, colour = variable)) +
  geom_line() + 
  labs(x = "MOnth of the year", y = "Total cases")

main_db[, rounded := round_any(timestep, 10, f = ceiling)]
tester <- main_db[,.N, by = c("rounded", "state")]

ggplot(tester, aes(x = rounded, y = N, colour = state)) + geom_line()
