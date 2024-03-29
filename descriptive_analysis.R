# this is spliced all together to start at the same timepoint
library(data.table)
library(plyr)
library(ggplot2)
library(gridExtra)
# Data setup, splicing seasons together
cc <- scales::seq_gradient_pal("darkmagenta",  "darkorange", "Lab")(seq(0,1,length.out=5))
colour1 <- "#fc8d62"
colour2 <- "#66c2a5"
colour3 <- "#8da0cb"

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
main_db$npsdatecol <- as.Date(main_db$npsdatecol, "%d-%m-%Y")
# remove the rows that have NA obsevations for flu/rsv (i.e. no data recoreded)
main_db <- na.omit(main_db)
# remove first season (the partial one)
#main_db <- main_db[npsdatecol >= as.Date("2017-01-16")]
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
(flu_tab/num_flu_infected)*100
### how many double positives for flu, assuming X days negative in between)
rsv_subset <- main_db[rsv == "Pos"]
num_rsv_infected <- length(unique(rsv_subset$indid))
rsv_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =indid]
rsv_tab <- table(table(rsv_subset[time_diff > num_days_neg]$indid))
print(sum(rsv_tab[-1])/num_rsv_infected*100)

summary_infections <- data.table("N_infections" = c(1,2,3),
                                 "RSV" = c(num_rsv_infected,unname(rsv_tab)[2:3]),
                                 "Flu" = c(num_flu_infected,unname(flu_tab)[2:3]))
summary_infections_m <- melt(summary_infections, id.vars = "N_infections",
                             measure.vars = c("RSV","Flu"))
summary_infections_m$N_infections <- factor(summary_infections_m$N_infections)
N_INFECTIONS <- ggplot(summary_infections_m,aes( x = variable, y = value, fill=N_infections)) +
  geom_bar(stat="identity", position ="dodge") + 
  theme_linedraw() + 
  labs(x = "Virus", fill = "No. Infections", y = "No. individuals", 
       title = "C") + 
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15)) + 
  scale_fill_manual(values=c( "#8da0cb","#b2df8a","#fb9a99"))

#recode to same format
main_db[flu=="Neg", flu := 0]
main_db[flu=="Pos", flu := 1]
main_db[rsv=="Neg", rsv := 0]
main_db[rsv=="Pos", rsv := 1]

main_db$flu <- as.numeric(main_db$flu)
main_db$rsv <- as.numeric(main_db$rsv)
summary_plot <- main_db[,sum(flu), by = npsdatecol]
summary_plot$rsv <-main_db[,sum(rsv), by = npsdatecol]$V1
colnames(summary_plot) <- c("Date", "Flu", "RSV")
summary_plot_m <- melt(summary_plot, id.vars = "Date")
INFECTIONS <- ggplot(summary_plot_m, aes(x = Date, y = value, colour = variable)) + 
  geom_line(size=1, alpha = 0.7) + 
  theme_linedraw() + 
  labs(y = "No. of positive tests by day", colour = "Virus", title = "A") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15)) + 
  scale_colour_manual(values=c(colour1,colour2))

main_db[,month := month(npsdatecol)]
main_db[,year := year(npsdatecol)]
summary_plot <- main_db[,sum(flu), by = c("month", "year")]
summary_plot$rsv <-main_db[,sum(rsv), by = c("month", "year")]$V1
colnames(summary_plot) <- c("month_num","year", "Flu", "RSV")
totals <- summary_plot[,sum(Flu), by = year]
totals$RSV <- summary_plot[,sum(RSV), by = year]$V1
summary_plot[totals, on="year", flu_tot := i.V1]
summary_plot[totals, on="year", rsv_tot := i.RSV]
summary_plot[, Flu :=( Flu/flu_tot)*100]
summary_plot[, RSV := (RSV/rsv_tot)*100]
summary_plot_m <- melt(summary_plot, id.vars = c("month_num", "year"))
summary_plot_m[,Month := month.abb[month_num]]
summary_plot_m$Month <- factor(summary_plot_m$Month, levels = month.abb)
plotter1 <- summary_plot[, mean(Flu), by = "month_num"]
plotter1$RSV <- summary_plot[, mean(RSV), by = "month_num"]$V1
colnames(plotter1) <- c("Month", "Flu", "RSV")
plotter1[,Month := month.abb[Month]]
plotter1 <- melt(plotter1, id.vars = "Month")
plotter1[Month == "Jan", value := NA]
plotter1 <- as.data.frame(plotter1)
plotter1 <- rbind(plotter1, c("Nov", "RSV", NA))
plotter1 <- rbind(plotter1, c("Nov", "Flu", NA))
plotter1 <- rbind(plotter1, c("Dec", "Flu", NA))
plotter1 <- rbind(plotter1, c("Dec", "RSV", NA))
plotter1$Month <- factor(plotter1$Month, levels = c("Jan", "Feb", "Mar", "Apr", 
                                                  "May", "Jun", "Jul", "Aug", 
                                                  "Sep", "Oct", "Nov", "Dec"))
# plotter1$Month <- factor(plotter1$Month)
plotter <- as.data.table(summary_plot_m[variable  %in% c("Flu", "RSV")])
plotter[Month == "Jan", value := NA]
plotter <- as.data.frame(plotter)
plotter <- rbind(plotter, c(11, 2018, "RSV", NA, "Nov"))
plotter <- rbind(plotter, c(12, 2018, "RSV", NA, "Dec"))
plotter$value <- as.numeric(plotter$value)
plotter$Month <- factor(plotter$Month, levels = c("Jan", "Feb", "Mar", "Apr", 
                                                  "May", "Jun", "Jul", "Aug", 
                                                  "Sep", "Oct", "Nov", "Dec"))
plotter1$variable <- as.factor(plotter1$variable)
plotter1$value <- as.numeric(plotter1$value)

# plotter$year <- factor(plotter$year)
plotter1 <- as.data.frame(plotter1)



SEASONALITY <- ggplot(plotter,
                      aes(x = Month, y = value, colour = variable, 
                          group=interaction(year,variable))) + 
  geom_rect(data = plotter, xmin = as.numeric(plotter$Month[[7]]), 
            xmax = as.numeric(plotter$Month[[10]])+0.75,
            ymin = 0, ymax = Inf, alpha = 0.01, fill = "lightgrey",color = NA) +
  geom_line(alpha = 0.5) +
   theme_linedraw() + 
   geom_line(data=plotter1, aes(x = Month, y = value, colour = variable, 
                                group = variable), size = 1.5)+
   labs(y = "Percentage of annual positive \n tests by month of the year", colour = "Virus", title = "D")+ 
   theme(axis.text = element_text(size = 15), 
         axis.title = element_text(size=15), 
         title = element_text(size=17), 
         legend.text = element_text(size=15))  + 
  scale_colour_manual(values=c(colour1, colour2)) + 
  geom_rect(data = plotter, xmin = as.numeric(plotter$Month[[52]])+0.25, 
            xmax = as.numeric(plotter$Month[[54]]),
            ymin = 0, ymax = Inf, alpha = 0.01, fill = "grey18",color = NA) + 
  geom_rect(data = plotter, xmin = as.numeric(plotter$Month[[7]]), 
            xmax = as.numeric(plotter$Month[[7]])+0.75,
            ymin = 0, ymax = Inf, alpha = 0.01, fill = "grey18",color = NA) 
  
SEASONALITY# ### respecify as positive those that are within 14 days of previous, 
# ### and remove the ones with subsequent re-infections
# repeat_rsvs <- rsv_subset[time_diff > num_days_neg]
# repeat_flus <- flu_subset[time_diff > num_days_neg]
# main_db$remove_repeat <- F
# # mark the flus to remove
# for(i in 1:nrow(repeat_flus)){
#   main_db[indid == unlist(repeat_flus[i,"indid"]) &
#             timestep >= unlist(repeat_flus[i,"timestep"]) , remove_repeat := T]
# }
# # mark the rsvs to remove
# for(i in 1:nrow(repeat_rsvs)){
#   main_db[indid == unlist(repeat_rsvs[i,"indid"]) &
#             timestep >= unlist(repeat_rsvs[i,"timestep"]) , remove_repeat := T]
#   
# }
# 
# main_db <- main_db[remove_repeat == F ]
main_db[, factor_id := as.numeric(factor(indid))]
# 
# 
# # convert them to all be the same year
# main_db[, combi_step := substring(npsdatecol,6)]
# main_db[, combi_step := paste0("2015-",combi_step)]
# min_artificial_date <- min(main_db$combi_step)
# main_db[,timestep := (as.Date(combi_step) - as.Date(min_artificial_date))]
# main_db[,timestep := as.numeric(timestep)]
# 
# # remove unneccesary columns
# main_db <- main_db[,c("factor_id","rsv", "flu", "timestep")]
# main_db$state <- "tbd"
# 
# 
# for(i in unique(main_db$factor_id)){
#   # subset for each person
#   temp <- main_db[factor_id == i,]
#   # calculate timings of infection
#   flu_time <- as.numeric(unlist(temp[flu == 1, "timestep"]))
#   rsv_time <- as.numeric(unlist(temp[rsv == 1, "timestep"]))
#   if(is.na(rsv_time[1])){rsv_time = Inf}
#   if(is.na(flu_time[1])){flu_time = Inf}
#   # specify the state at each timepoint
#   temp[timestep < flu_time[1] & timestep < rsv_time[1], state := "SS"]
#   temp[timestep > flu_time[1] & timestep < rsv_time[1], state := "SPR"]
#   temp[timestep < flu_time[1] & timestep >rsv_time[1], state := "PRS"]
#   temp[timestep > flu_time[1] & timestep > rsv_time[1], state := "RR"]
#   temp[timestep == flu_time[1] & timestep< rsv_time[1], state := "SI"]
#   temp[timestep < flu_time[1] & timestep == rsv_time[1], state := "IS"]
#   temp[timestep > flu_time[1] & timestep == rsv_time[1], state := "IPR"]
#   temp[timestep == flu_time[1] & timestep > rsv_time[1], state := "PRI"]
#   temp[timestep == flu_time[1] & timestep == rsv_time[1], state := "II"]
#   temp[timestep > tail(flu_time,1) & timestep > tail(rsv_time,1), state := "RR"]
#   
#   main_db[factor_id == i,] <- temp
# }
# # recode to numeric
# main_db[state == "SS", num_state :=1]
# main_db[state == "IS", num_state :=2]
# main_db[state == "PRS", num_state :=99]
# main_db[state == "SI", num_state :=5]
# main_db[state == "II", num_state :=6]
# main_db[state == "PRI", num_state :=7]
# main_db[state == "SPR", num_state :=101]
# main_db[state == "IPR", num_state :=9]
# main_db[state == "RR", num_state :=11]
# 
# #save(main_db, file="main_db_together_april_nopartial.Rdata")
# 
# # Number infected first with flu and then RSV
# length(unique(main_db[state %in%c ("IPR")]$factor_id))
# # Number infected first with RSV and then flu
# length(unique(main_db[state %in%c ("PRI")]$factor_id))


main_db


  
tiff("Data_basics.tiff", width=1300, height=800)
grid.arrange(INFECTIONS, N_INFECTIONS, SEASONALITY, AGE_PLOT, 
             layout_matrix = rbind(c(1,1,1,1,4,4),
                                   c(2,2,3,3,3,3)))
dev.off()

