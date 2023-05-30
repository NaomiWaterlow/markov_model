# crude analysis
library(plyr)
library(zoo)
library(ggplot2)
library(gridExtra)
library(binom)
colour1 <- "#fc8d62"
colour2 <- "#66c2a5"
colour3 <- "#8da0cb"

# load in the data
library(data.table)
PHIRST_NP <- data.table(read.csv("~/Documents/OneDrive - London School of Hygiene and Tropical Medicine/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Nasopharyngeal specimens 2016-2018 Flu & RSV 2022-04-19.csv",
                                 na.strings = ""))
PHIRST_BACKGROUND <- data.table(read.csv("~/Documents/OneDrive - London School of Hygiene and Tropical Medicine/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Master 2016-2018 Flu & RSV 2022-04-19.csv"))
# format correctly
main_db <- PHIRST_NP[,c("indid", "npsdatecol", "flu", "rsv", "npsinfa", "npsinfb")]
main_db[,flutype := "unknown"]
main_db[npsinfa == "1",flutype := "A"]
main_db[npsinfb == "1",flutype := "B"]
main_db[, npsinfa := NULL]
main_db[, npsinfb := NULL]
colnames(PHIRST_BACKGROUND)[1] <- "indid"
main_db[PHIRST_BACKGROUND, on = c("indid"),age_grp := i.age_cat_at_consent]

main_db$indid <- as.factor(main_db$indid)
main_db[, factor_id := as.numeric(factor(indid))]
main_db$npsdatecol <- as.Date(main_db$npsdatecol, "%d-%m-%Y")

#main_db <- main_db[!is.na(flu) | !is.na(rsv)]

# earliest date recprded
start_date <- min(main_db$npsdatecol)
# change from date to timestp
main_db[,timestep := (npsdatecol - start_date)]
main_db[,timestep := as.numeric(timestep)]

#recode to same format
main_db[flu=="Neg", flu := 0]
main_db[flu=="Pos", flu := 1]
main_db[rsv=="Neg", rsv := 0]
main_db[rsv=="Pos", rsv := 1]

### respecify as positive those that are within 14 days of previous, 

# Attack rate
main_db[, season := 2]
main_db[npsdatecol < as.Date("2016-11-01"), season := 1]
main_db[npsdatecol > as.Date("2018-01-01"), season := 3]

##### Investigate the numbers that get infected multiple times ####
num_days_neg <- 14
flu_subset <- main_db[flu == 1]
num_flu_infected <- length(unique(flu_subset$factor_id))
flu_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =factor_id]
flu_tab <- table(table(flu_subset[time_diff > num_days_neg]$factor_id))

print(paste0(round(sum(flu_tab)/num_flu_infected*100), "% of infected individuals with Flu, get Flu at least once more"))
# 88 flu reinfections

### how many double positives for flu, assuming X days negative in between)
rsv_subset <- main_db[rsv == 1]
# number of people who get an infection 
num_rsv_infected <- length(unique(rsv_subset$factor_id))
rsv_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =factor_id]
rsv_tab <- table(table(rsv_subset[time_diff > num_days_neg]$factor_id))
# number of peopple who have more than one infection / number of people who get infected at least once
print(paste0(round(sum(rsv_tab)/num_rsv_infected*100), "% of infected individuals with RSV, get RSV at least once more"))
# 43 RSV reinfections

### mark the subsequent re-infections - htis is REINFECTIONS marked, NOT REPEAT SWABS FROM SAME INFECTION
repeat_rsvs <- rsv_subset[time_diff > num_days_neg]
repeat_flus <- flu_subset[time_diff > num_days_neg]
main_db$not_primary_inf_flu <- F
main_db$not_primary_inf_rsv <- F
# mark the flus susbequent infections
for(i in 1:nrow(repeat_flus)){
  main_db[indid == unlist(repeat_flus[i,"indid"]) &
            timestep >= unlist(repeat_flus[i,"timestep"]) , not_primary_inf_flu := T]
}
# mark the rsvs subsequent infections
for(i in 1:nrow(repeat_rsvs)){
  main_db[indid == unlist(repeat_rsvs[i,"indid"]) &
            timestep >= unlist(repeat_rsvs[i,"timestep"]) , not_primary_inf_rsv := T]
  
}

# attack rate <- number of individual infections / population at risk
# at the same time, mark whichis the first swab for each infection
for(i in 1:3){
  rsv_subset[, first_swab := T]
  rsv_subset[time_diff<=num_days_neg, first_swab := F]
  rsv_first_swabs <- rsv_subset[season == i,sum(first_swab)]
  pop_at_risk <- nrow(unique(main_db[season==i,"factor_id"]))
  print(paste0("RSV attack rate for season ", i, " is ", round((rsv_first_swabs/pop_at_risk)*100), "%"))
  print(paste0(i, "season: ", rsv_first_swabs, "/", pop_at_risk ))
}
for(i in 1:3){
  flu_subset[, first_swab := T]
  flu_subset[time_diff<=num_days_neg, first_swab := F]
  flu_first_swabs <- flu_subset[season == i,sum(first_swab)]
  pop_at_risk <- nrow(unique(main_db[season==i,"factor_id"]))
  print(paste0("Influenza attack rate for season ", i, " is ", round((flu_first_swabs/pop_at_risk)*100), "%"))
  print(paste0(i, "season: ", rsv_first_swabs, "/", pop_at_risk ))
}
print(paste("Note, the first season is only a partial season and that this analysis includes repeat infections as seperate infections if more than ", num_days_neg, " days apart"))

#add the first swab details back into the main database
main_db[, first_swab_flu := F]
main_db[, first_swab_rsv := F]
main_db[, flu_first_swab_and_infection := F]
main_db[, rsv_first_swab_and_infection := F]
main_db[flu_subset, on = c("indid", "npsdatecol"), first_swab_flu := i.first_swab]
main_db[rsv_subset, on = c("indid", "npsdatecol"), first_swab_rsv := i.first_swab]
main_db[not_primary_inf_flu == F & first_swab_flu ==T, flu_first_swab_and_infection := T]
main_db[not_primary_inf_rsv == F & first_swab_rsv ==T, rsv_first_swab_and_infection := T]

#### CALCULATE OVERLAPPING EPISODES #### 

# this should be false for overlapping episodes, but true for extended overlapping episodes (additional 7 days after)
add_7_days <- F

 main_db_duos <- copy(main_db)
 main_db_duos[,considered_flu := F]
 main_db_duos[,considered_rsv := F]
 
 if(add_7_days ==T){
   add_amount <- 7
 } else {add_amount <- 0}
 
 # for each person
 for(i in unique(main_db_duos$indid)){
   # find start and end of the episodes
   temp <- main_db_duos[indid == i]
   
  #### RSV
   
   #find the positives
   temp_pos <- temp[rsv ==1]
   # need to run for each episode!
   temp_first_episodes <- c(unlist(temp_pos[first_swab_rsv == T, "timestep"]), Inf)
   
   for(j in 1:(length(temp_first_episodes)-1)){
     
     temp_times <- temp_pos[timestep < temp_first_episodes[j+1] &
                              timestep >= temp_first_episodes[j]]$timestep
     
     start_ep <- head(temp_times,1)
     end_ep <- tail(temp_times,1) +7
     
     if(!length(start_ep)==0 & !length(end_ep)==0){
       # if have both values do all between
       main_db_duos[indid ==i & timestep %in% c(start_ep:end_ep), 
                    considered_rsv := T]
     } 
   }
   
   
   ##### FLU
   
   #find the positives
   temp_pos <- temp[flu ==1]
   # need to run for each episode!
   temp_first_episodes <- c(unlist(temp_pos[first_swab_flu == T, "timestep"]), Inf)
   
   for(j in 1:(length(temp_first_episodes)-1)){
     
     temp_times <- temp_pos[timestep < temp_first_episodes[j+1] &
                              timestep >= temp_first_episodes[j]]$timestep
     
     start_ep <- head(temp_times,1)
     end_ep <- tail(temp_times,1) +7
     
     if(!length(start_ep)==0 & !length(end_ep)==0){
       # if have both values do all between
       main_db_duos[indid ==i & timestep %in% c(start_ep:end_ep), 
                    considered_flu := T]
     } 
   }
 }
 
#Total overlapping infections!
length(unique(main_db_duos[considered_flu == T & considered_rsv == T]$indid))
 
# 
# for(i in 1:nrow(flu_subset[first_swab==T])){
#   t_indid <- flu_subset[first_swab==T][i,"indid"]
#   t_start <- flu_subset[first_swab==T][i,"timestep"]
#   main_db_duos[indid == unlist(t_indid) & timestep %in% c(unlist(t_start):(unlist(t_start)+14)), 
#                considered_flu := T]
# }
# 
# for(i in 1:nrow(rsv_subset[first_swab==T])){
#   t_indid <- rsv_subset[first_swab==T][i,"indid"]
#   t_start <- rsv_subset[first_swab==T][i,"timestep"]
#   main_db_duos[indid == unlist(t_indid) & timestep %in% c(unlist(t_start):(unlist(t_start)+14)), 
#                considered_rsv := T]
# }
# 
# # total number of individuals with a potentially overlapping infection.   
# nrow(unique(main_db_duos[considered_flu == T & considered_rsv == T, "indid"]))

# make to weekly (over 7 days)
main_db[, rounded := round_any(timestep, 7, f = ceiling)]
# #calculate the number per week.
 weekly <- data.table( main_db[,sum(first_swab_rsv), by=c("rounded", "factor_id", "season")])
weekly_flu<- data.table(main_db[,sum(first_swab_flu), by=c("rounded", "factor_id", "season")])
weekly[weekly_flu, flu := i.V1, on=c("rounded", "factor_id", "season")]
 colnames(weekly)[4] <- "rsv"

 # age by season
 age_sumamry <- main_db[first_swab_flu==T,.N, by = c("season", "age_grp")]
 age_sumamry_rsv <- main_db[first_swab_rsv==T,.N, by = c("season", "age_grp")]
age_sumamry[age_sumamry_rsv, on = c("season", "age_grp"), RSV := i.N] 
colnames(age_sumamry)[3] <- "Flu"
age_sumamry_m <- melt.data.table(age_sumamry, id.vars = c("season", "age_grp"))
age_sumamry_m$age_grp <- factor(age_sumamry_m$age_grp, levels = c(
  "< 5 years","5 - 18 years", "19 - 65 years", ">65 years"
))
age_sumamry_m[season ==1, season1 := "Season 1"]
age_sumamry_m[season ==2, season1 := "Season 2"]
age_sumamry_m[season ==3, season1 := "Season 3"]



AGE_PLOT <- ggplot(age_sumamry_m, aes(x = age_grp, y = value, fill = variable)) + 
  facet_grid(season1~.) + 
  geom_bar(stat="identity", position = "dodge") + 
  theme_linedraw() +
  scale_fill_manual(values=c("#fc8d62", "#66c2a5")) +
  labs(x = "Age group", y = "Number of episodes", fill = "Infection", 
       title = "B") + 
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=15),
        strip.background = element_rect(fill = "azure4"), 
        legend.position = "none")

 
rsv_swab_daily <- data.table("current_date" = NA, 
                          "rsv_current" = NA, 
                          "flu_target" = NA, 
                          "target_date"= NA)

season_limits <- data.frame( c(1,2,3), 
                             c(1, 266, 623), 
                             c(182, 546 ,917))
colnames(season_limits) = c("season", "start", "end")


days_of_interest <- c(-21:21)
n_samples <- 20
set.seed(63)

main_db_firsts <- main_db[rsv_first_swab_and_infection == T | 
                            flu_first_swab_and_infection ==T]
firsts_rsv <- main_db_firsts[rsv == 1]
firsts_flu <- main_db_firsts[flu == 1]
no_controls <- c()

#####rsv calculations ######
#for each positive first swab
for(pos_sample in 1:nrow(firsts_rsv)){
  
  target_date <- unlist(firsts_rsv[pos_sample,"timestep"])
  target_age <- unlist(firsts_rsv[pos_sample, "age_grp"])
  positive_id <- unlist(firsts_rsv[pos_sample, "factor_id"])
  season_look <- unlist(firsts_rsv[pos_sample, "season"])
  season_dates <- season_limits[season_look, 2:3]
  
  #find 5 controls on same date and same age but negative
  controls_all <-  main_db[timestep == target_date &
                             age_grp == target_age &
                             rsv == 0,]
  no_controls <- c(no_controls, length(unique(controls_all$factor_id)))

  if(no_controls[pos_sample] >= n_samples){
  controls_sub <- unique(controls_all[sample(x = c(1:nrow(controls_all)),
                                             size = n_samples),])$factor_id
  } else { controls_sub <- unique(controls_all$factor_id)}
  ids_to_run <-  c(positive_id, controls_sub)
  
  
  # for each week of interest method.
  for(j in 1:length(c(days_of_interest))){
    search_day <- days_of_interest[j]
    
    # need to remove the weeks that don't exist#
    # e.g. can't look -3 backs if youre cases in the first week of data
    if(!(target_date + search_day) <= season_dates[1]& !(target_date + search_day) >= season_dates[2] ){
      # for each individual
      for(k in 1:length(ids_to_run)){
        # subset only that individuals data
        temp <- main_db[factor_id ==ids_to_run[k],]
        this_day_flu <- NA
        target_day_rsv <- NA
        # In the current day, what is rsv?
        if(temp[timestep ==target_date,"first_swab_rsv"] == 0){this_day_rsv <- 0
        } else if(temp[timestep ==target_date,"first_swab_rsv"] == 1){this_day_rsv <- 1
        }
        # in the target week what is flu
        if((target_date+search_day) %in% c(temp$timestep)){
        if(temp[timestep == (target_date + search_day), "first_swab_flu"] == 0){target_day_flu <- 0
        } else if(temp[timestep == (target_date + search_day),
                       "first_swab_flu"] == 1){target_day_flu <- 1} } else{
                         target_day_flu <- NA
                       }
        temp2 <- data.frame(target_date, this_day_rsv, target_day_flu, search_day)
        colnames(temp2) <- colnames(rsv_swab_daily)
        
        rsv_swab_daily <- rbind(rsv_swab_daily,temp2 )
      }}
  }
  
}

rsv_swab_daily <- na.omit(rsv_swab_daily)
# when there is no rsv, what proportion of the target day have flu
rsv0 <- rsv_swab_daily[rsv_current == 0, (sum(flu_target)/.N), by = "target_date"]
# when there is rsv, what proportion of the targt day have flu
rsv1 <- rsv_swab_daily[rsv_current == 1, (sum(flu_target)/.N), by = "target_date"]
colnames(rsv0) = colnames(rsv1) = c("target_date", "prop_in_target")
rsv0$infected <- F
rsv1$infected <- T

rsv0$N <- rsv_swab_daily[rsv_current == 0, .N, by = "target_date"]$N
rsv1$N <- rsv_swab_daily[rsv_current == 1, .N, by = "target_date"]$N

rsv_together <- rbind(rsv0,rsv1)
rsv_together$base_virus <- "rsv"

rsv_swab_daily[target_date %in% c(0:6) , target_week := 0  ]
rsv_swab_daily[target_date %in% c(7:13) , target_week := 1  ]
rsv_swab_daily[target_date %in% c(14:20) , target_week := 2  ]
#rsv_swab_daily[target_date %in% c(18:24) , target_week := 3  ]
rsv_swab_daily[target_date %in% c(-1:-7) , target_week := -1  ]
rsv_swab_daily[target_date %in% c(-8:-14) , target_week := -2 ]
rsv_swab_daily[target_date %in% c(-15:-21) , target_week := -3 ]

rsv_swab_daily <- na.omit(rsv_swab_daily)
# function for calculating risk ratio
calc_rr_rsv_exposure <- function(rsv_swab_daily, no_week){

tester <- rsv_swab_daily[target_week == no_week,]
tester_tab <-data.table(table(tester[,c("rsv_current", "flu_target")]))

a = tester_tab[which(tester_tab$rsv_current==1 & tester_tab$flu_target ==1),]$N
b = tester_tab[which(tester_tab$rsv_current==1 & tester_tab$flu_target ==0),]$N
c = tester_tab[which(tester_tab$rsv_current==0 & tester_tab$flu_target ==1),]$N
d = tester_tab[which(tester_tab$rsv_current==0 & tester_tab$flu_target ==0),]$N

if(length(a) != 1){ a = 0}
if(length(b) != 1){ b = 0}
if(length(c) != 1){ c = 0}
if(length(d) != 1){ d = 0}

rr <- (a/(a+b)) / (c/(c+d))
lower_ln <- exp(log(rr) + 1.96* sqrt( ((b/a)/(b+a)) + ((d/c)/(d+c))))
upper_ln <- exp(log(rr) - 1.96* sqrt(  ((b/a)/(b+a)) + ((d/c)/(d+c))))

if(!is.finite(rr)){rr <- NA}
if(length(lower_ln)<1){lower_ln <- NA}
if(length(upper_ln)<1){upper_ln <- NA}

if(length(c(rr, lower_ln, upper_ln))<3){combo <- c(rr,NA,NA,no_week)} else{
  combo = c(rr, min(lower_ln, upper_ln), max(lower_ln, upper_ln), no_week)}

return(combo)}

# setup for risk ratios
rr_store <- data.frame(matrix(ncol = 4, nrow = length(unique(rsv_swab_daily$target_week))))
colnames(rr_store) <- c("risk ratio", "lower", "upper", "target_week")

# run risk ratios
for(i in 1:length(unique(rsv_swab_daily$target_week))){
  no_week <- unique(rsv_swab_daily$target_week)[i]
  temp <- calc_rr_rsv_exposure(rsv_swab_daily, no_week = no_week)
  rr_store[i,] <- temp
}

###### flu calculations #####


flu_swab_daily <- data.table("current_date" = NA, 
                             "flu_current" = NA, 
                             "rsv_target" = NA, 
                             "target_date"= NA)
no_controls <- c()
#for each positive first swab
for(pos_sample in 1:nrow(firsts_flu)){
  
  target_date <- unlist(firsts_flu[pos_sample,"timestep"])
  target_age <- unlist(firsts_flu[pos_sample, "age_grp"])
  positive_id <- unlist(firsts_flu[pos_sample, "factor_id"])
  season_look <- unlist(firsts_flu[pos_sample, "season"])
  season_dates <- season_limits[season_look, 2:3]
  
  #find 5 controls on same date and same age but negative
  controls_all <-  main_db[timestep == target_date &
                             age_grp == target_age &
                             flu == 0,]
  no_controls <- c(no_controls, length(unique(controls_all$factor_id)))
  
  if(no_controls[pos_sample] >= n_samples){
    controls_sub <- unique(controls_all[sample(x = c(1:nrow(controls_all)),
                                               size = n_samples),])$factor_id
  } else { controls_sub <- unique(controls_all$factor_id)}
  ids_to_run <-  c(positive_id, controls_sub)
  
  # for each week of interest method.
  for(j in 1:length(c(days_of_interest))){
    search_day <- days_of_interest[j]
    
    # need to remove the weeks that don't exist#
    # e.g. can't look -3 backs if youre cases in the first week of data
    if(!(target_date + search_day) <= season_dates[1]& !(target_date + search_day) >= season_dates[2] ){
      # for each individual
      for(k in 1:length(ids_to_run)){
        # subset only that individuals data
        temp <- main_db[factor_id ==ids_to_run[k],]
        this_day_rsv <- NA
        target_day_flu <- NA
        # In the current dte, what is rsv?
        if(temp[timestep ==target_date,"first_swab_flu"] == 0){this_day_flu <- 0
        } else if(temp[timestep ==target_date,"first_swab_flu"] == 1){this_day_flu <- 1}
        # in the target week what is flu
        if((target_date+search_day) %in% c(temp$timestep)){
          if(temp[timestep == (target_date + search_day), "first_swab_rsv"] == 0){target_day_rsv <- 0
          } else if(temp[timestep == (target_date + search_day),
                         "first_swab_rsv"] == 1){target_day_rsv<- 1} } else{
                           target_day_rsv <- NA
                         }
        
        temp2 <- data.frame(target_date, this_day_flu, target_day_rsv, search_day)
        colnames(temp2) <- colnames(flu_swab_daily)
        flu_swab_daily <- rbind(flu_swab_daily,temp2 )
      }}
  }
  
}

flu_swab_daily <- na.omit(flu_swab_daily)
flu0 <- flu_swab_daily[flu_current == 0, (sum(rsv_target)/.N), by = "target_date"]
flu1 <- flu_swab_daily[flu_current == 1, (sum(rsv_target)/.N), by = "target_date"]

colnames(flu0) = colnames(flu1) = c("target_date", "prop_in_target")
flu0$infected <- F
flu1$infected <- T

flu0$N <- flu_swab_daily[flu_current == 0, .N, by = "target_date"]$N
flu1$N <- flu_swab_daily[flu_current == 1, .N, by = "target_date"]$N
flu_together <- rbind(flu0,flu1)
flu_together$base_virus <- "flu"


# flu_swab_daily[target_date %in% c(-3:3) , target_week := 0  ]
# flu_swab_daily[target_date %in% c(4:10) , target_week := 1  ]
# flu_swab_daily[target_date %in% c(11:17) , target_week := 2  ]
# flu_swab_daily[target_date %in% c(18:24) , target_week := 3  ]
# flu_swab_daily[target_date %in% c(-4:-10) , target_week := -1  ]
# flu_swab_daily[target_date %in% c(-11:-17) , target_week := -2 ]
# flu_swab_daily[target_date %in% c(-18:-24) , target_week := -3 ]

flu_swab_daily[target_date %in% c(0:6) , target_week := 0  ]
flu_swab_daily[target_date %in% c(7:13) , target_week := 1  ]
flu_swab_daily[target_date %in% c(14:20) , target_week := 2  ]
#flu_swab_daily[target_date %in% c(18:24) , target_week := 3  ]
flu_swab_daily[target_date %in% c(-1:-7) , target_week := -1  ]
flu_swab_daily[target_date %in% c(-8:-14) , target_week := -2 ]
flu_swab_daily[target_date %in% c(-15:-21) , target_week := -3 ]

flu_swab_daily <- na.omit(flu_swab_daily)
# function for calculating risk ratio
calc_rr_flu_exposure <- function(flu_swab_daily, no_week){
  
  tester <- flu_swab_daily[target_week == no_week,]
  
  tester_tab <-data.table(table(tester[,c("flu_current", "rsv_target")]))
  a = tester_tab[which(tester_tab$flu_current==1 & tester_tab$rsv_target ==1),]$N
  b = tester_tab[which(tester_tab$flu_current==1 & tester_tab$rsv_target ==0),]$N
   c = tester_tab[which(tester_tab$flu_current==0 & tester_tab$rsv_target ==1),]$N
  d = tester_tab[which(tester_tab$flu_current==0 & tester_tab$rsv_target ==0),]$N

  if(length(a) != 1){ a = 0}
  if(length(b) != 1){ b = 0}
  if(length(c) != 1){ c = 0}
  if(length(d) != 1){ d = 0}

  rr <- (a/(a+b)) / (c/(c+d))
  lower_ln <- exp(log(rr) + 1.96* sqrt( ((b/a)/(b+a)) + ((d/c)/(d+c))))
  upper_ln <- exp(log(rr) - 1.96* sqrt(  ((b/a)/(b+a)) + ((d/c)/(d+c))))

  if(!is.finite(rr)){rr <- NA}
  if(length(lower_ln)<1){lower_ln <- NA}
  if(length(upper_ln)<1){upper_ln <- NA}

  if(length(c(rr, lower_ln, upper_ln))<3){combo <- c(rr,NA,NA,no_week)} else{
   combo = c(rr, min(lower_ln, upper_ln), max(lower_ln, upper_ln), no_week)}
  
  return(combo)}

# setup for risk ratios
rr_store_flu <- data.frame(matrix(ncol = 4, nrow = length(unique(flu_swab_daily$target_week))))
colnames(rr_store_flu) <- c("risk ratio", "lower", "upper", "target_week")

# run risk ratios
for(i in 1:length(unique(flu_swab_daily$target_week))){
  no_week <- unique(flu_swab_daily$target_week)[i]
  
  temp <- calc_rr_flu_exposure(flu_swab_daily, no_week = no_week)
  rr_store_flu[i,] <- temp
}

rr_store_flu$exposure_virus <- "influenza"
rr_store$exposure_virus <- "rsv"

rr_together <- rbind(rr_store, rr_store_flu)

# need to do smoothing?
colnames(rr_together)[1] <- "rr"
#rr_together <- na.omit(rr_together)
rr_together$exposure_virus <- factor(rr_together$exposure_virus, levels = c("rsv", "influenza"))
rr_together <- data.table(rr_together)
rr_together <- rr_together[order(exposure_virus,target_week)]

saveRDS(rr_together, file="rr_together.rds")


rr_together<- na.omit(rr_together)

rr_together[target_week == 0, label_x := "0:6"]
rr_together[target_week == 1, label_x := "7:13"]
rr_together[target_week == 2, label_x := "14:20"]
rr_together[target_week == -1, label_x := "-7:-1"]
rr_together[target_week == -2, label_x := "-14:-8"]
rr_together[target_week == -3, label_x := "-21:-15"]

rr_together[exposure_virus == "rsv", exposure_virus := "RSV"]
rr_together[exposure_virus == "influenza", exposure_virus := "Influenza"]

rr_together$label_x <- factor(rr_together$label_x, levels = 
                                c("-21:-15", "-14:-8", "-7:-1", 
                                  "0:6", "7:13", "14:20"))

# initials fo a plot
RISKRATIO <- ggplot(rr_together, aes(x = label_x, y = rr, colour = exposure_virus, group = exposure_virus)) + 
  geom_line()+
  geom_errorbar(aes(ymin = lower, ymax = upper, group = exposure_virus), 
                width =0.3, alpha = 0.8) +
  theme_linedraw() + 
  geom_point() +
  labs( x = "Days relative to swab of virus", y = "Risk ratio of being infected with other virus", 
        title = "A", colour = "Exposure virus", fill = "Exposure virus") + 
  scale_colour_manual(values = c(colour2, colour1)) + 
  scale_fill_manual(values = c(colour2, colour1)) +
 facet_grid(.~exposure_virus)+
  theme(axis.text = element_text(size = 15), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size = 15)) +
   geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(trans='log10') +
  geom_hline(yintercept = 1, linetype = "dashed") 

RISKRATIO
# 
# 
# # below this is for the analysis without risk ratios. No longer used 
# together <- rbind(flu_together,rsv_together)
# together[, target := N*prop_in_target ]
# together$mean <- 1000
# together$upper <- 1000
# together$lower <- 1000
# 
# for(i in 1:nrow(together)){
#   cis <- binom.confint(x = as.numeric(together[i,"target"]), 
#                        n = as.numeric(together[i, "N"]), 
#                        methods = "exact")
#   together[i,"mean"] <- cis$mean
#   together[i,"lower"] <- cis$lower
#   together[i,"upper"] <- cis$upper
#   
#   
# }
# 
# 
# together[base_virus =="flu", base_virus := "Influenza"]
# together[base_virus =="rsv", base_virus := "RSV"]
# 
# together$base_virus <- factor(together$base_virus, levels = c("RSV", "Influenza"))
# 
# together <- together[order(base_virus,infected,target_date)]
# 
# no_to_roll <- 3
# together$prop_rolled <- c(rep(NA,no_to_roll-2),rollmean(together$mean, k=no_to_roll),rep(NA,no_to_roll-2))
# together$lower_rolled <- c(rep(NA,no_to_roll-2),rollmean(together$lower, k=no_to_roll),rep(NA,no_to_roll-2))
# together$upper_rolled <- c(rep(NA,no_to_roll-2),rollmean(together$upper, k=no_to_roll),rep(NA,no_to_roll-2))
# 
# 
# 
# 
# 
# PERCENT <- ggplot(together, aes(x = target_date, y = mean ,fill = infected)) +
#   geom_line(aes(colour = infected)) +geom_point(aes(colour = infected)) + 
#   geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.3) +
#   #geom_bar(stat="identity", position = "dodge") +
#   facet_grid(.~base_virus ) + 
#   theme_linedraw() + 
#   labs( x = "Day relative to swab of virus", y = "Proportion infected with other virus", 
#         title = "A", colour = "Virus infection", fill = "Virus infection") + 
#   scale_colour_manual(values = c("gray20", colour3)) + 
#   scale_fill_manual(values = c("gray20", colour3)) +
#   theme(axis.text = element_text(size = 15), 
#         axis.title = element_text(size=15), 
#         title = element_text(size=17), 
#         legend.text = element_text(size=15), 
#         strip.text = element_text(size = 15)) 
# 
# PERCENT
# # for each individual in the data set
# for(j in 1:length(c(swab_weekly$week))){
#   for(i in 1:length(unique(weekly$factor_id))){
#   # subset the relevant steps
#   search_week <- swab_weekly$week[j]
#   temp <- weekly[factor_id == i]
#   rsv_weeks <- which(temp$rsv >0)
#   flu_weeks <- which(temp$flu >0)
#   last_step <- nrow(temp)
#   # exclude any infections that are too close to the start or end to be included
#   if(search_week <0){
#    exclude <-  c(0:-(search_week))
#   } else {
#     exclude <- c(last_step:(last_step-(search_week)))
#   }
#   # for each rsv positive swab
#     for(k in rsv_weeks){
#       # if the relevant rsv week is 
#       if(!(k %in% exclude)) {  
#       # add any positive swabs to the target week total
#         swab_weekly[j,"number_flu"] <-  swab_weekly[j,"number_flu"] + temp[k+search_week,]$flu
#         # add the number of swabs this week to the total looked at
#         swab_weekly[j,"total_rsv"] <-  swab_weekly[j,"total_rsv"] + temp[k,]$rsv
# 
#       } 
#     }
#   
#   # for each flu positive swab
#   for(k in flu_weeks){
#     # if the relevant rsv week is 
#     if(!(k %in% exclude)) {  
#       # add any positive swabs to the target week total
#       swab_weekly[j,"number_rsv"] <-  swab_weekly[j,"number_rsv"] + temp[k+search_week,]$rsv
#       # add the number of swabs this week to the total looked at
#       swab_weekly[j,"total_flu"] <-  swab_weekly[j,"total_flu"] + temp[k,]$flu
#       
#     } 
#   }
# }
# }
# 
# swab_weekly[,RSV := (number_flu/total_rsv)*100]
# swab_weekly[,influenza := (number_rsv/total_flu)*100]
# 
# swab_weekly_m <- melt.data.table(swab_weekly, id.vars = c("week"),
#                                  measure.vars=c("RSV", "influenza"))
# swab_weekly_m$week <- as.factor(swab_weekly_m$week)
# 
#  RELATIVE_WEEKS <- ggplot(swab_weekly_m, aes(x = variable, y = value, group = week)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#  # geom_line(size = 1.5) +  
#   theme_linedraw() +
#   labs(x = "Week relative to positive swab of base virus", 
#        y = "Percentage infected with other virus in relevant week", 
#        colour = "Base virus", title = "A") +
#   geom_vline(xintercept = 0, linetype = "dotted") + 
#    scale_colour_manual(values = c(colour2, colour1)) +
#    theme(axis.text = element_text(size = 15), 
#          axis.title = element_text(size=15), 
#          title = element_text(size=17), 
#          legend.text = element_text(size=15))
# 
#  RELATIVE_WEEKS

######## Incidence of dual positive swabs ######

# - need incidence of RSV by week
# - need inciddence of flu by  week
# - if independent, chance of having both is the multiple?
# - so can compare incidence of dual swab against incidence of indep*incidence of dual
# - however this will be weekly/daily... can sum total incidence of each meauser of the whole time period? 

# try with daily for now. 
main_db$rsv <- as.numeric(main_db$rsv)
main_db$flu <- as.numeric(main_db$flu)

# by timestep, number of positive swabs by virus
incidence_comparator <- main_db[, sum(rsv, na.rm = T), by =c("npsdatecol", "season")]
colnames(incidence_comparator)[3] <- "RSV"
incidence_comparator$flu <- main_db[, sum(flu, na.rm = T), by =c("npsdatecol", "season")]$V1
# who has dual infections?
temp <- main_db[rsv==1 & flu==1]
temp <- temp[,sum(rsv), by = "npsdatecol"]
#mark how many of the incidences were actually dual incidences
incidence_comparator[, Observed := 0]
incidence_comparator[temp, Observed := i.V1, on=c("npsdatecol") ]
N_swabs <- main_db[, N_swabs :=  .N, by = "npsdatecol"]
incidence_comparator[N_swabs, on = c("npsdatecol"), N_swabs := i.N_swabs]

# Need the total number of swabs that day

incidence_comparator[,incidence_RSV := (RSV/N_swabs)]
incidence_comparator[,incidence_FLU := flu/N_swabs]
incidence_comparator[,Observed_incidence := Observed/N_swabs]
incidence_comparator[, Expected_incidence := incidence_FLU*incidence_RSV]
incidence_comparator[,Expected := Expected_incidence*N_swabs]


N_observed_dual <- sum(incidence_comparator$Observed)
N_expected_dual <-sum(incidence_comparator$Expected)
N_swabs <-   sum(incidence_comparator$N_swabs)

test_out <- prop.test(x = c(N_observed_dual, N_expected_dual), n = c(N_swabs, N_swabs))
test_out
(N_expected_dual + (test_out$conf.int[1]*N_swabs)) /N_expected_dual
(N_expected_dual + (test_out$conf.int[2]*N_swabs)) /N_expected_dual

  # incidence_comparator[, upper := (prop_in_target + z*sqrt(prop_in_target * (1-prop_in_target)/N))*100]


incidence_comparator_m <- melt(incidence_comparator, 
                               id.vars = "npsdatecol", 
                               measure.vars = c("Observed", "Expected"))
library(ISOweek)
incidence_comparator_m[, rounded := ISOweek(npsdatecol)]
plotter <- incidence_comparator_m[, sum(value), by = c("variable", "rounded")]
plotter[, cumulative := cumsum(V1), by = "variable"]
plotter[,date := ISOweek2date(paste0(rounded, "-1"))]


EXPECTED <- ggplot(plotter, aes(x=date, y = V1, colour = variable)) + 
  scale_y_continuous(
     name = "Weekly number of dual positive swabs (solid)", 
     sec.axis = sec_axis(trans = ~.*3, name = "Cumulative dual positive swabs (dashed)")
  ) +
  geom_line(size = 1.2) +
  geom_line(aes(y = cumulative/3), size = 1, linetype =2)+
  labs(y = "Weekly number of positive swabs for both influenza and RSV", 
       x = "Date", 
       colour = "Type", title = "B") + 
  theme_linedraw() + 
  scale_colour_manual(values = c("#8da0cb","gray20" )) +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15))
EXPECTED  

tiff("Crude_analysis.tiff", width=1200, height=800)
grid.arrange(RISKRATIO, EXPECTED,
             layout_matrix = rbind(c(1,1),
                                   c(2,2)))
dev.off()



flu1 <- unique(main_db[season==1 &flu == 1]$indid)
flu1_val <- round(PHIRST_BACKGROUND[indid %in% flu1, median(age_at_consent)])
flu1_min <- round(PHIRST_BACKGROUND[indid %in% flu1, min(age_at_consent)])
flu1_max <- round(PHIRST_BACKGROUND[indid %in% flu1, max(age_at_consent)])

flu2 <- unique(main_db[season==2 &flu == 1]$indid)
flu2_val <- round(PHIRST_BACKGROUND[indid %in% flu2, median(age_at_consent)])
flu2_min <- round(PHIRST_BACKGROUND[indid %in% flu2, min(age_at_consent)])
flu2_max <- round(PHIRST_BACKGROUND[indid %in% flu2, max(age_at_consent)])

flu3 <- unique(main_db[season==3 &flu == 1]$indid)
flu3_val <- round(PHIRST_BACKGROUND[indid %in% flu3, median(age_at_consent)])
flu3_min <- round(PHIRST_BACKGROUND[indid %in% flu3, min(age_at_consent)])
flu3_max <- round(PHIRST_BACKGROUND[indid %in% flu3, max(age_at_consent)])

rsv1 <- unique(main_db[season==1 &rsv == 1]$indid)
rsv1_val <- round(PHIRST_BACKGROUND[indid %in% rsv1, median(age_at_consent)])
rsv1_min <- round(PHIRST_BACKGROUND[indid %in% rsv1, min(age_at_consent)])
rsv1_max <- round(PHIRST_BACKGROUND[indid %in% rsv1, max(age_at_consent)])

rsv2 <- unique(main_db[season==2 &rsv == 1]$indid)
rsv2_val <- round(PHIRST_BACKGROUND[indid %in% rsv2, median(age_at_consent)])
rsv2_min <- round(PHIRST_BACKGROUND[indid %in% rsv2, min(age_at_consent)])
rsv2_max <- round(PHIRST_BACKGROUND[indid %in% rsv2, max(age_at_consent)])

rsv3 <- unique(main_db[season==3 &rsv == 1]$indid)
rsv3_val <- round(PHIRST_BACKGROUND[indid %in% rsv3, median(age_at_consent)])
rsv3_min <- round(PHIRST_BACKGROUND[indid %in% rsv2, min(age_at_consent)])
rsv3_max <- round(PHIRST_BACKGROUND[indid %in% rsv2, max(age_at_consent)])

print(paste0("season 1: flu median age is ", flu1_val, " and  rsv mean age is ", rsv1_val))
print(paste0("range flu is ", flu1_min, " to ", flu1_max))
print(paste0("range rsv is ", rsv1_min, " to ", rsv1_max))
print(paste0("season 2: flu median age is ", flu2_val, " and  rsv mean age is ", rsv2_val))
print(paste0("range flu is ", flu2_min, " to ", flu2_max))
print(paste0("range rsv is ", rsv2_min, " to ", rsv2_max))
print(paste0("season 3: flu median age is ", flu3_val, " and  rsv mean age is ", rsv3_val))
print(paste0("range flu is ", flu3_min, " to ", flu3_max))
print(paste0("range rsv is ", rsv3_min, " to ", rsv3_max))






library(ISOweek)
# plot states over time
timer_plot <- main_db[,.N, by= c("state", "timestep")]
timer_plot[, rounded := ISOweek(npsdatecol)]
timer_plot[, N := sum(N), by = c("state", "rounded")]
timer_plot[, total := sum(N), by = "rounded"]
timer_plot[, prop := N/total]
timer_plot[, rounded := paste0(rounded, "-1")]
timer_plot[,Date := ISOweek2date(rounded)]


ggplot(timer_plot, aes(x = timestep, y = N, colour = state)) + 
  geom_point()





tmp <- rle(flu_subset$first_swab)
flu_subset$Seq_flu <- rep(tmp$lengths >= 1,times = tmp$lengths)
# add the sequence number to each
seq_to_add <- c()
start_seq <- 1
for(i in 1:length(tmp$lengths)){
  
  length_run <- tmp$lengths[i]
  tester <- flu_subset[sum(tmp$lengths[1:i]),"first_swab"]
  if(tester == F){ 
    seq_to_add <- c(seq_to_add, rep(0,length_run))
  } else {
    seq_to_add <- c(seq_to_add, rep(start_seq,length_run))
    start_seq <- start_seq + 1
  }
}

flu_subset$seq_flu <- seq_to_add

for(i in 1:nrow(flu_subset)){
  if(flu_subset[i, "seq_flu"]==0){
  flu_subset[i, "seq_flu"] <-flu_subset[i-1, "seq_flu"] }
}

# if the time difference is bigger than specified value then it's acutally reset
flu_subset[time_diff >= num_days_neg, time_diff := 0]

quantile(flu_subset[, sum(time_diff, na.rm = T), by = "seq_flu"]$V1, probs=c(0.025, 0.5, 0.975))




tmp <- rle(rsv_subset$first_swab)
rsv_subset$Seq_rsv <- rep(tmp$lengths >= 1,times = tmp$lengths)
# add the sequence number to each
seq_to_add <- c()
start_seq <- 1
for(i in 1:length(tmp$lengths)){
  
  length_run <- tmp$lengths[i]
  tester <- rsv_subset[sum(tmp$lengths[1:i]),"first_swab"]
  if(tester == F){ 
    seq_to_add <- c(seq_to_add, rep(0,length_run))
  } else {
    seq_to_add <- c(seq_to_add, rep(start_seq,length_run))
    start_seq <- start_seq + 1
  }
}

rsv_subset$seq_rsv <- seq_to_add

for(i in 1:nrow(rsv_subset)){
  if(rsv_subset[i, "seq_rsv"]==0){
    rsv_subset[i, "seq_rsv"] <-rsv_subset[i-1, "seq_rsv"] }
}

# if the time difference is bigger than specified value then it's acutally reset
rsv_subset[time_diff >= num_days_neg, time_diff := 0]

quantile(rsv_subset[, sum(time_diff, na.rm = T), by = "seq_rsv"]$V1, probs=c(0.025, 0.5, 0.975))


