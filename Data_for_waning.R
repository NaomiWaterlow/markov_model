## Data setup, all in one with waning
library(data.table)
library(ggplot2)
library(msm)
library(readr)

PHIRST_NP <- data.table(read_csv("~/Documents/OneDrive - London School of Hygiene and Tropical Medicine/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Nasopharyngeal specimens 2016-2018 Flu & RSV 2022-04-19.csv"))
MASTER_PHIRST <- data.table(read_csv("~/Documents/OneDrive - London School of Hygiene and Tropical Medicine/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Master 2016-2018 Flu & RSV 2022-04-19.csv"))
# format correctly
main_db <- PHIRST_NP[,c("indid", "npsdatecol", "flu", "rsv")]
colnames(MASTER_PHIRST)[1] <- "indid"
main_db[MASTER_PHIRST, on="indid", age_cat := i.age_cat_at_consent]

main_db$indid <- as.factor(main_db$indid)
main_db$age_cat <- as.factor(main_db$age_cat)
main_db$npsdatecol <- as.Date(main_db$npsdatecol, "%d-%m-%Y")
# remove the rows that have NA obsevations for flu/rsv (i.e. no data recoreded)
main_db <- na.omit(main_db)
# earliest date recprded
start_date <- min(main_db$npsdatecol)
# convert to numeric time
main_db[,timestep := (npsdatecol - start_date)]
main_db[,timestep := as.numeric(timestep)]

#recode to same format
main_db[flu=="Neg", flu := 0]
main_db[flu=="Pos", flu := 1]
main_db[rsv=="Neg", rsv := 0]
main_db[rsv=="Pos", rsv := 1]

main_db[, factor_id := as.numeric(factor(indid))]

# remove unneccesary columns
main_db <- main_db[,c("npsdatecol","factor_id","rsv", "flu", "age_cat", "timestep")]
# recode to numeric
main_db[rsv==0 & flu ==0, state := 101]
main_db[rsv==1 & flu ==0, state := 102]
main_db[rsv==0 & flu ==1, state := 103]
main_db[rsv==0 & flu ==1, state := 104]

main_db$rsv <- as.numeric(main_db$rsv)
main_db$flu <- as.numeric(main_db$flu)

# format to create a plot
to_plot <- main_db[ ,sum(flu), by = timestep]
to_plot$RSV <- main_db[ ,sum(rsv), by = timestep]$V1
colnames(to_plot) <- c("timestep", "flu", "rsv")
to_plot_m <- melt.data.table(to_plot, id.vars=c("timestep"))
ggplot(to_plot_m, aes(x = timestep, y = value, colour = variable )) + 
  geom_line()

#### create the covariate based on weekly time
main_db[, week := week(npsdatecol)]
main_db[, year := year(npsdatecol)]
main_db[, week_year := paste0(year, "_", week)]
# calculate total of each by week
background <- main_db[, mean(rsv), by=week_year]
background$flu <- main_db[,mean(flu), by=week_year]$V1
colnames(background) <- c("week_year", "rsv_week", "flu_week")

main_db[background, on ="week_year", rsv_week := i.rsv_week]
main_db[background, on ="week_year", flu_week := i.flu_week]
main_db$rsv_week <- as.numeric(main_db$rsv_week)
ggplot(main_db, aes(x = timestep, y = flu_week)) + geom_line()

# better age group format
main_db[age_cat == "< 5 years", age_cat := "under_5"]
main_db[age_cat == "5 - 18 years", age_cat := "5_to_18"]
main_db[age_cat == "19 - 65 years", age_cat := "19_to_65"]
main_db[age_cat == ">65 years", age_cat := "over_65"]
main_db$age_cat <- factor(main_db$age_cat, levels = c("under_5", "5_to_18", "19_to_65", "over_65"))
# the final data to input into the model
input_data <- main_db[,c("factor_id", "timestep", "state", "rsv_week", "flu_week", "age_cat")]


save(input_data, file= "data_waning.Rdata")



