# Evaluate the output
library(ggplot2)
library(msm)
library(data.table)
# see the overall output 
output

pci_splits <- c(30,60,90,120,150,181,
                211,241,271,301,331,361,391,421,451,467,
                497,527,557,587,617,647,677,707,737)

# what is the impact of interaction?
# strength
load("~/Documents/GitHub/markov_model/july_fixed_7_rep.Rdata")
when_inf_first <- qratio.msm(output, ind1 = c(5,6), ind2 = c(1,2))
when_rsv_first <- qratio.msm(output, ind1 = c(2,6), ind2 = c(1,5))
print(paste0("Transmission from S to I with RSV is ", round(when_inf_first[1],2)," (", round(when_inf_first[3],2),
             " - ", round(when_inf_first[4],2),
             ") more likely when already infected with flu"))
print(paste0("Transmission from S to I with flu is ", round(when_rsv_first[1],2)," (", round(when_rsv_first[3],2),
             " - ", round(when_rsv_first[4],2),
             ") more likely when already infected with rsv"))
when_inf_first_7 <- c(when_inf_first, virus_infected = "Influenza", Interaction_duration = "7 days")
when_rsv_first_7 <- c(when_rsv_first, virus_infected = "RSV", Interaction_duration = "7 days")

load("~/Documents/GitHub/markov_model/july_fixed_14_rep.Rdata")
when_inf_first <- qratio.msm(output, ind1 = c(5,6), ind2 = c(1,2))
when_rsv_first <- qratio.msm(output, ind1 = c(2,6), ind2 = c(1,5))
print(paste0("Transmission from S to I with RSV is ", round(when_inf_first[1],2)," (", round(when_inf_first[3],2),
             " - ", round(when_inf_first[4],2),
             ") more likely when already infected with flu"))
print(paste0("Transmission from S to I with flu is ", round(when_rsv_first[1],2)," (", round(when_rsv_first[3],2),
             " - ", round(when_rsv_first[4],2),
             ") more likely when already infected with rsv"))
when_inf_first_14 <- c(when_inf_first, virus_infected = "Influenza", Interaction_duration = "14 days")
when_rsv_first_14 <- c(when_rsv_first, virus_infected = "RSV", Interaction_duration = "14 days")

load("~/Documents/GitHub/markov_model/july_fixed_21_rep.Rdata")
when_inf_first <- qratio.msm(output, ind1 = c(5,6), ind2 = c(1,2))
when_rsv_first <- qratio.msm(output, ind1 = c(2,6), ind2 = c(1,5))
print(paste0("Transmission from S to I with RSV is ", round(when_inf_first[1],2)," (", round(when_inf_first[3],2),
             " - ", round(when_inf_first[4],2),
             ") more likely when already infected with flu"))
print(paste0("Transmission from S to I with flu is ", round(when_rsv_first[1],2)," (", round(when_rsv_first[3],2),
             " - ", round(when_rsv_first[4],2),
             ") more likely when already infected with rsv"))
when_inf_first_21 <- c(when_inf_first, virus_infected = "Influenza", Interaction_duration = "21 days")
when_rsv_first_21 <- c(when_rsv_first, virus_infected = "RSV", Interaction_duration = "21 days")


load("~/Documents/GitHub/markov_model/july_fixed_28_rep.Rdata")
when_inf_first <- qratio.msm(output, ind1 = c(5,6), ind2 = c(1,2))
when_rsv_first <- qratio.msm(output, ind1 = c(2,6), ind2 = c(1,5))
print(paste0("Transmission from S to I with RSV is ", round(when_inf_first[1],2)," (", round(when_inf_first[3],2),
             " - ", round(when_inf_first[4],2),
             ") more likely when already infected with flu"))
print(paste0("Transmission from S to I with flu is ", round(when_rsv_first[1],2)," (", round(when_rsv_first[3],2),
             " - ", round(when_rsv_first[4],2),
             ") more likely when already infected with rsv"))
when_inf_first_28 <- c(when_inf_first, virus_infected = "Influenza", Interaction_duration = "28 days")
when_rsv_first_28 <- c(when_rsv_first, virus_infected = "RSV", Interaction_duration = "28 days")



# add the gaps back into the pci splits
pci_splits[which(pci_splits>182)] <- pci_splits[which(pci_splits>182)]+77
pci_splits[which(pci_splits>546)] <- pci_splits[which(pci_splits>545)]+78


interaction_info <- data.table(rbind(when_inf_first_7,when_rsv_first_7, 
                                     when_inf_first_14,when_rsv_first_14,
                                     when_inf_first_21,when_rsv_first_21,
                                     when_inf_first_28,when_rsv_first_28))
interaction_info$estimate <- as.numeric(interaction_info$estimate)
interaction_info$se <- as.numeric(interaction_info$se)
interaction_info$L <- as.numeric(interaction_info$L)
interaction_info$U <- as.numeric(interaction_info$U)
interaction_info$Interaction_duration <- factor(interaction_info$Interaction_duration, 
                                                levels = rev(c("7 days", "14 days", "21 days", 
                                                           "28 days")))

INTERACTION <- ggplot(interaction_info, aes(x = virus_infected, y = estimate, colour = Interaction_duration
                            , alpha = Interaction_duration
                             ))+ 
  geom_point(position = position_dodge(width =0.4))+
  geom_pointrange(aes(ymin = L, ymax = U), position = position_dodge(width =0.4),size=1.2) + 
  theme_linedraw() + 
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15)) +
  labs(x = "Initial virus infection", y = "Interaction multiplier", colour = "Interaction duration", 
       alpha = "Interaction duration", title = "C") +
  coord_flip() + 
  scale_colour_manual(values=rev(c(colour3,"lightsalmon3", "orchid4", "deeppink2")))+ 
  scale_alpha_manual(values=rev(c(1,0.4,0.4,0.4)) )

# What is the seasonality. Can plot over time? 
# hazard.msm creates a list of the hazard rates for each pci value
# For each time there are two parameters estimated
# extract from model

load("~/Documents/GitHub/markov_model/july_fixed_7_rep2.Rdata")
timing_hazards_list <- hazard.msm(output)
# combine into data table
timing_hazards <- data.table(names_times = names(timing_hazards_list)[4:28])
timing_hazards[,names_t := gsub(pattern = "timeperiodRSVt",replacement = "",x= names_times)]
timing_hazards$names_t <- as.numeric(timing_hazards$names_t)
timing_hazards$RSV <- sapply(timing_hazards_list[4:28], function(x) x[1,1])
timing_hazards$Influenza <- sapply(timing_hazards_list[29:53],function(x) x[2,1])

timing_hazards <- timing_hazards[,.SD[order(names_t)]]
timing_hazards$from <- c(pci_splits)
timing_hazards[,Date := as.Date(from, origin = "2016-05-02")]

timing_cis <- data.table(
  names_times = names(timing_hazards_list)[4:28]
)
timing_cis[,names_t := gsub(pattern = "timeperiodRSVt",replacement = "",x= names_times)]
timing_cis$names_t <- as.numeric(timing_cis$names_t)
timing_cis$RSV <- sapply(timing_hazards_list[4:28], function(x) x[1,2])
timing_cis$Influenza <- sapply(timing_hazards_list[29:53], function(x) x[2,2])
timing_cis <- timing_cis[,.SD[order(names_t)]]
timing_cis$from <- c(pci_splits)
timing_cis[,Date := as.Date(from, origin = "2016-05-02")]

timing_cis_up <- data.table(
  names_times = names(timing_hazards_list)[4:28]
)
timing_cis_up[,names_t := gsub(pattern = "timeperiodRSVt",replacement = "",x= names_times)]
timing_cis_up$names_t <- as.numeric(timing_cis_up$names_t)
timing_cis_up$RSV <- sapply(timing_hazards_list[4:28],function(x) x[1,3])
timing_cis_up$Influenza <- sapply(timing_hazards_list[29:53], function(x) x[2,3])
timing_cis_up <- timing_cis_up[,.SD[order(names_t)]]
timing_cis_up$from <- c(pci_splits)
timing_cis_up[,Date := as.Date(from, origin = "2016-05-02")]

# make long
timing_hazards_m <- melt.data.table(timing_hazards, 
                                    id.vars = c("names_t", "from", "Date"), 
                                    measure.vars = c("Influenza", "RSV"))
timing_cis_m <- melt.data.table(timing_cis, 
                                    id.vars = c("names_t", "from", "Date"), 
                                    measure.vars = c("Influenza", "RSV"))
timing_cis_up_m <- melt.data.table(timing_cis_up, 
                                id.vars = c("names_t", "from", "Date"), 
                                measure.vars = c("Influenza", "RSV"))

timing_hazards_m[timing_cis_m, on=c("names_t", "from", "Date", "variable"), 
                 lower:= i.value]
timing_hazards_m[timing_cis_up_m, on=c("names_t", "from", "Date", "variable"), 
                 upper:= i.value]
timing_hazards_m[,year := year(Date)]
timing_hazards_m[,week := isoweek(Date)]


# plot the output 
FOIS <- ggplot(timing_hazards_m, aes(x = week, y = log(value), colour = variable)) + 
  geom_line()+
  geom_point() + 
  theme_linedraw() + 
  facet_grid(year~.) +
  lims(x=c(5,45))+
  geom_ribbon(aes(ymin = log(lower), ymax = log(upper), fill = variable), alpha = 0.3, colour =NA) + 
  labs(y= " Log Hazard ratio compared to first month of data", fill = "Force of Infection", 
       colour = "Force of Infection", x = "Week of the year", title = "A") + 
  scale_colour_manual(values = c(colour1, colour2)) +
  scale_fill_manual(values =c(colour1,colour2)) + 
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15), 
        legend.position = "bottom") 
FOIS

FOISA <- FOIS 
FOISB <- FOIS 

###### hazards for age groups

load("~/Documents/GitHub/markov_model/july_fixed_7_rep2.Rdata")
# combine into data table
age_hazards <- data.table(names_times = names(timing_hazards_list)[1:3])
age_hazards$hazard <- sapply(timing_hazards_list[1:3], function(x) x[1,1])
age_hazards$lower <- sapply(timing_hazards_list[1:3],function(x) x[1,2])
age_hazards$upper <- sapply(timing_hazards_list[1:3],function(x) x[1,3])
age_hazards$type <- "RSV"

age_hazards1 <- data.table(names_times = names(timing_hazards_list)[1:3])
age_hazards1$hazard <- sapply(timing_hazards_list[1:3], function(x) x[2,1])
age_hazards1$lower <- sapply(timing_hazards_list[1:3],function(x) x[2,2])
age_hazards1$upper <- sapply(timing_hazards_list[1:3],function(x) x[2,3])
age_hazards1$type <- "Flu"

age_hazards <- rbind(age_hazards, age_hazards1)

age_hazards[names_times == "age_grpage2", label_age := "5 - 18 years"]
age_hazards[names_times == "age_grpage3", label_age := "19 - 65 years"]
age_hazards[names_times == "age_grpage4", label_age := ">65 years"]
age_hazards <- rbind(age_hazards,data.frame(
  names_times ="age_grpage1",
  hazard =1,
  lower =1,
  upper =1,
  type = "RSV",
  label_age ="< 5 years"))
age_hazards <- rbind(age_hazards,data.frame(
  names_times ="age_grpage1",
  hazard =1,
  lower =1,
  upper =1,
  type = "Flu",
  label_age ="< 5 years"))

age_hazards$label_age <- factor(age_hazards$label_age, levels = 
                                  c("< 5 years", "5 - 18 years","19 - 65 years",">65 years"))
age_hazards$label_age <- factor(age_hazards$label_age, levels = rev(levels(age_hazards$label_age)))

AGES <- ggplot(age_hazards, aes(x = label_age, y = log(hazard), group= type, colour = type)) + 
  geom_point( size =1.5,position=position_dodge(width = 0.30))+
  geom_pointrange(aes(ymin = log(lower), ymax = log(upper)),
                  position=position_dodge(width = 0.30),
                  size =1.5) + 
  theme_linedraw() + 
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15)) +
  scale_colour_manual(values = c(colour1, colour2 ))+
  labs(x = "Age group", y = "Log Hazard ratio compared to <5 year olds", title="B") +
  coord_flip() + theme(plot.margin=margin(10,20,10,10))

# AGES_COMBI <- ggplot(age_hazards, aes(x = label_age, y = log(hazard), fill= type, colour = type)) + 
#   geom_point( size =1.5, position = position_dodge(width =0.4),)+
#   geom_pointrange(aes(ymin = log(lower), ymax = log(upper)),  size =1.5,
#                   position = position_dodge(width = 0.4)) + 
#   theme_linedraw() + 
#   theme(axis.text = element_text(size = 15), 
#         axis.title = element_text(size=15), 
#         title = element_text(size=17), 
#         legend.text = element_text(size=15)) +
#   labs(x = "Age group", y = "Log Hazard ratio compared to <5 year olds", title="C", colour = "Replicate", 
#        fill = "Replicate") +
#   coord_flip() + theme(plot.margin=margin(10,20,10,10)) 


AGES_COMBI

tiff("Supplement_21.TIFF", width=1200, height=800)
grid.arrange(FOIS, AGES, layout_matrix = rbind(c(1,2)))
dev.off()


BLANK <- ggplot() + theme_minimal() + labs(title ="A") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15)) + 
  theme(plot.margin=margin(10,10,10,10))


tiff("Model_output.TIFF", width=1200, height=800)
grid.arrange(FOIS, AGES, BLANK, INTERACTION, layout_matrix = rbind(c(3,3,3,2,2,2),
                                                                   c(3,3,3,2,2,2),
                                                                   c(3,3,3,1,1,1),
                                                                   c(4,4,4,1,1,1),
                                                                   c(4,4,4,1,1,1)))
dev.off()

tiff("Supplement_rho70.TIFF", width=1200, height=800)
grid.arrange(FOIS, AGES, layout_matrix = rbind(c(1,1,1,2,2,2), 
                                                      c(1,1,1,2,2,2)))
dev.off()



# What is the misclassifciation

# goodness of model?
plot.prevalence.msm(output)


# work out the actual FOI values over time, rather than the hazard ratios

timing_fois <- data.table(from = rep(pci_splits,2), 
                             to = rep(c(pci_splits[-1],NA),2), 
                          variable = c("tbd"))



for(i in 1:length(pci_splits)){
  start <- as.character(pci_splits[i])
  end <- as.character(pci_splits[i+1])
  if(is.na(end)){end <- "Inf"}
  variable_label <- paste0("[",start,",", end,")" )
 
  temp <- qmatrix.msm(output, covariates=list(timeperiod=variable_label))[1,2]
  timing_fois[i,"variable"] <- "RSV"
  timing_fois[i,c("value")] <- temp[1]
  timing_fois[i,c("lower")] <- temp[3]
  timing_fois[i,c("upper")] <- temp[4]
  
  temp <- qmatrix.msm(output, covariates=list(timeperiod=variable_label))[1,5]
  timing_fois[i+length(pci_splits),"variable"] <- "Influenza"
  timing_fois[i+length(pci_splits),c("value")] <- temp[1]
  timing_fois[i+length(pci_splits),c("lower")] <- temp[3]
  timing_fois[i+length(pci_splits),c("upper")] <- temp[4]
  
}

timing_fois[,Date := as.Date(from, origin = "2016-05-02")]
#timing_fois[,Date := from]
# plot the output 
TIMING_FOIS <- ggplot(timing_fois, aes(x = Date, y = value, colour = variable)) + 
  geom_line()+
  geom_point() + 
  theme_linedraw() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), alpha = 0.3, colour =NA) + 
  labs(y= "Force of Infection (daily)", fill = "Virus", 
       colour = "Virus", x = "Month", title = "A" ) +
  scale_color_manual(values = c("darkorange1", "darkolivegreen3")) +
  scale_fill_manual(values = c("darkorange2", "darkolivegreen3"))+
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15)) + 
  scale_fill_manual(values =c(colour3))
  



BLANK <- ggplot() + theme_minimal() + labs(title ="A") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15))

source("survivial plot.R")
 SURVIVAL <- plot.msm(output)

 library(gridExtra)

 
tiff("Timing_default.TIFF", width=1200, height=800)
grid.arrange(SURVIVAL, TIMING_FOIS, BLANK,layout_matrix = rbind(c(2,2),
                                                          c(3,1)))
dev.off()

tiff("Timing_linked.TIFF", width=1200, height=800)
grid.arrange(SURVIVAL, TIMING_FOIS, layout_matrix = rbind(c(2,1)
                                                                ))
dev.off()


# total parameters
# 7 main parameters, but 2 fixed = 5
# 2*lengthpci = 18 (50 for in season)
# 1 misclassifciation = 1
# 11 initial state parameters = 11

n_par <- 67
aic = output[[2]] + 2*n_par 
aic

####
duration_out <- qmatrix.msm(output)[3,4]
paste0("duration is ", round(1/duration_out["estimate"]), " (", round(1/duration_out["upper"]), " - ", round(1/duration_out["lower"]), ") days")


pci_splits <- c(30,60,90,120,150,181,
                211,241,271,301,331,361,391,421,451,467,
                497,527,557,587,617,647,677,707,737)


# don't like the plot prevalence. Want to make it more flexible and be able to 
# see the small percentages and the confidence intervals etc. 

#therefore make a table for each one.

prev3 <- prevalence.msm(x = output, 
                     #   timezero = 468,
                        times = seq(468, 757,7),
                        piecewise.times = pci_splits[which(pci_splits>497)], 
                        piecewise.covariates = list(timeperiodRSV = c("t16","t17", "t18", 
                                                                      "t19", "t20", "t21", "t22", 
                                                                      "t23", "t24", "t25"),
                                                    timeperiodFlu = c("t16","t17", "t18", 
                                                                      "t19", "t20", "t21", "t22", 
                                                                      "t23", "t24", "t25")))

prev2 <- prevalence.msm(x = output, 
                     #   timezero = 182,
                        times=seq(182, 467, 7),
                        piecewise.times = c(211,241,271,301,331,361,391,421,451), 
                        piecewise.covariates = list(timeperiodRSV = c("t6","t7", "t8", 
                                                                      "t9", "t10", "t11", "t12", 
                                                                      "t13", "t14", "t15"),
                                                    timeperiodFlu = c("t6","t7", "t8", 
                                                                      "t9", "t10", "t11", "t12", 
                                                                      "t13", "t14", "t15")))

prev1 <- prevalence.msm(x = output, 
                        timezero = 0,
                        times = seq(0,181,7),
                        piecewise.times = pci_splits[which(pci_splits<181 )], 
                        piecewise.covariates = list(timeperiodRSV = c("t0","t1","t2", "t3", 
                                                                      "t4", "t5"),
                                                    timeperiodFlu = c("t0","t1","t2", "t3", 
                                                                      "t4", "t5")))


prev1_Estimated <- data.table(prev1$Expected)
prev1_Estimated$time <- 1:nrow(prev1_Estimated)
prev1_Observed <- data.table(prev1$Observed)
prev1_Observed$time <- 1:nrow(prev1_Observed)

colnames(prev1_Observed) <- colnames(prev1_Estimated)

prev1_Observed[, Influenza := SI + II + PRI]
prev1_Estimated[, Influenza := SI + II + PRI]
prev1_Observed[, RSV := IS + II + IPR]
prev1_Estimated[, RSV := IS + II + IPR]

prev1_observed <- melt(prev1_Observed, id.vars = "time")
prev1_observed$type = "observed"
prev1_observed$season = "Season 1"
prev1_Estimated <- melt(prev1_Estimated, id.vars = "time")
prev1_Estimated$type = "estimated"
prev1_Estimated$season = "Season 1"



prev1_Estimated[, Date := (time*7) - 6]
prev1_Estimated[, Date := as.Date(Date, origin = "2016-05-02")]
prev1_observed[, Date := (time*7) - 6]
prev1_observed[, Date := as.Date(Date, origin = "2016-05-02")]

prev2_Estimated <- data.table(prev2$Expected)
prev2_Estimated$time <- 1:nrow(prev2_Estimated)
prev2_Observed <- data.table(prev2$Observed)
prev2_Observed$time <- 1:nrow(prev2_Observed)

colnames(prev2_Observed) <- colnames(prev2_Estimated)

prev2_Observed[, Influenza := SI + II + PRI]
prev2_Estimated[, Influenza := SI + II + PRI]
prev2_Observed[, RSV := IS + II + IPR]
prev2_Estimated[, RSV := IS + II + IPR]

prev2_observed <- melt(prev2_Observed, id.vars = "time")
prev2_observed$type = "observed"
prev2_observed$season = "Season 2"
prev2_Estimated <- melt(prev2_Estimated, id.vars = "time")
prev2_Estimated$type = "estimated"
prev2_Estimated$season = "Season 2"

prev2_observed[, Date := (time*7) - 6]
prev2_observed[, Date := as.Date(Date, origin = "2017-01-16")]
prev2_Estimated[, Date := (time*7) - 6]
prev2_Estimated[, Date := as.Date(Date, origin = "2017-01-16")]

prev3_Estimated <- data.table(prev3$Expected)
prev3_Estimated$time <- 1:nrow(prev3_Estimated)
prev3_Observed <- data.table(prev3$Observed)
prev3_Observed$time <- 1:nrow(prev3_Observed)

colnames(prev3_Observed) <- colnames(prev3_Estimated)

prev3_Observed[, Influenza := SI + II + PRI]
prev3_Estimated[, Influenza := SI + II + PRI]
prev3_Observed[, RSV := IS + II + IPR]
prev3_Estimated[, RSV := IS + II + IPR]

prev3_observed <- melt(prev3_Observed, id.vars = "time")
prev3_observed$type = "observed"
prev3_observed$season = "Season 3"
prev3_Estimated <- melt(prev3_Estimated, id.vars = "time")
prev3_Estimated$type = "estimated"
prev3_Estimated$season = "Season 3"

prev3_observed[, Date := (time*7) - 6]
prev3_observed[, Date := as.Date(Date, origin = "2018-01-15")]
prev3_Estimated[, Date := (time*7) - 6]
prev3_Estimated[, Date := as.Date(Date, origin = "2018-01-15")]

prev_combo <- rbind(prev1_observed, prev1_Estimated,
                     prev2_observed, prev2_Estimated,
                     prev3_observed, prev3_Estimated)



prev_combo[variable == "State 1", variable := "SS"]
prev_combo[variable == "State 2", variable := "IS"]
prev_combo[variable == "State 3", variable := "PS"]
prev_combo[variable == "State 4", variable := "RS"]
prev_combo[variable == "State 5", variable := "SI"]
prev_combo[variable == "State 6", variable := "II"]
prev_combo[variable == "State 7", variable := "PRI"]
prev_combo[variable == "State 8", variable := "SP"]
prev_combo[variable == "State 9", variable := "IPR"]
prev_combo[variable == "State 10", variable := "SR"]
prev_combo[variable == "State 11", variable := "RR"]

prev_combo$variable <- factor(prev_combo$variable, levels = 
                                c("Influenza", "RSV", "SI","PRI", "IS", "IPR", "II","PS", "SP", "RS", "SR", "SS", 
                                  "RR", "Total"))


MODEL_FIT <- ggplot(prev_combo[type == "observed"], aes(x = Date, y = value, colour = type, group = season)) + 
  facet_grid(variable ~ ., scales = "free_y") + geom_point(size =3) + 
  geom_line(data = prev_combo[type == "estimated"])+
  theme_linedraw() + 
  labs(x = "Date", y=  "Number of individuals in each compartment", colour = "Data") + 
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size = 14))+ 
  scale_colour_manual(values = c("gray20",colour3))

MODEL_FIT_INFECTIONS <- ggplot(prev_combo[type == "observed" & variable %in% c("Influenza", 
                                                                    "RSV", 
                                                                    "SI","PRI", "IS", "IPR", "II")], 
                               aes(x = Date, y = value, colour = type, group = season)) + 
  facet_grid(variable ~ ., scales = "free_y") +
  geom_point(size = 3) + 
  geom_line(data = prev_combo[type == "estimated"& variable %in% c("Influenza", 
                                                                   "RSV", 
                                                                   "SI","PRI", "IS", "IPR", "II")])+
  theme_linedraw() + 
  labs(x = "Date", y=  "Number of individuals in each compartment", 
       colour = "Data") + 
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size=15), 
        title = element_text(size=17), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size = 14)) + 
  scale_colour_manual(values = c("gray20",colour3))

tiff("Model_fit_Infections.TIFF", width=1200, height=800)
MODEL_FIT_INFECTIONS
dev.off()

tiff("Model_fit.TIFF", width=1200, height=900)
MODEL_FIT
dev.off()

