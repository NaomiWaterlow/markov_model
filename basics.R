# crude analysis

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

### how many double positives for flu, assuming X days negative in between)
rsv_subset <- main_db[rsv == 1]
# number of people who get an infection 
num_rsv_infected <- length(unique(rsv_subset$factor_id))
rsv_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =factor_id]
rsv_tab <- table(table(rsv_subset[time_diff > num_days_neg]$factor_id))
# number of peopple who have more than one infection / number of people who get infected at least once
print(paste0(round(sum(rsv_tab)/num_rsv_infected*100), "% of infected individuals with RSV, get RSV at least once more"))

# probably want an attack rate by season, given only have half a season at the start? 

# attack rate <- number of individual infections / population at risk
for(i in 1:3){
  rsv_subset[, first_swab := T]
  rsv_subset[time_diff<num_days_neg, first_swab := F]
  rsv_first_swabs <- rsv_subset[season == i,sum(first_swab)]
  pop_at_risk <- nrow(unique(main_db[season==i,"factor_id"]))
  print(paste0("RSV attack rate for season ", i, " is ", round((rsv_first_swabs/pop_at_risk)*100), "%"))
}
for(i in 1:3){
  flu_subset[, first_swab := T]
  flu_subset[time_diff<num_days_neg, first_swab := F]
  flu_first_swabs <- flu_subset[season == i,sum(first_swab)]
  pop_at_risk <- nrow(unique(main_db[season==i,"factor_id"]))
  print(paste0("Influenza attack rate for season ", i, " is ", round((flu_first_swabs/pop_at_risk)*100), "%"))
}
print("Note, the first season is only a partial season and that this includes repeat infections")



