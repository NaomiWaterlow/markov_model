# Creating the data to be used

#### TODO - when RSV data is cleaned, confirm NAs vs nos. Currently assuming all NA is a Negative for RSV


#### load in the data and subset + basic manip ####

PHIRST_NP <- data.table(read.csv("~/Documents/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Nasopharyngeal specimens 2016-2018 Flu & RSV 2021-10-19.csv",
                                 na.strings = ""))
# format correctly
main_db <- PHIRST_NP[,c("indid", "npsdatecol", "flu", "npsrsva", "npsrsvb")]
main_db[npsrsva == 1 | npsrsvb == 1, rsv := 1]
main_db$indid <- as.factor(main_db$indid)
main_db$npsdatecol <- as.Date(main_db$npsdatecol, "%d-%m-%Y")
# earliest date recprded
start_date <- min(main_db$npsdatecol)
# convert to numeric time
main_db[,timestep := (npsdatecol - start_date)]
main_db[,timestep := as.numeric(timestep, units="days")]

##### Investigate the numbers that get infected multiple times ####
num_days_neg <- 14
#
flu_subset <- main_db[flu == "Pos"]
num_flu_infected <- length(unique(flu_subset$indid))
flu_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =indid]
flu_tab <- table(table(flu_subset[time_diff > num_days_neg]$indid))
print(sum(flu_tab[-1])/num_flu_infected*100)

### how many double positives for flu, assuming X days negative in between)
rsv_subset <- main_db[rsv == 1]
num_rsv_infected <- length(unique(rsv_subset$indid))
rsv_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =indid]
rsv_tab <- table(table(rsv_subset[time_diff > num_days_neg]$indid))
print(sum(rsv_tab[-1])/num_rsv_infected*100)


#### Look at data ####
#recode to same format
observed_only[is.na(rsv), rsv := 0]
observed_only[is.na(flu), flu2 := 0]
observed_only[flu=="Neg", flu2 := 0]
observed_only[flu=="Pos", flu2 := 1]
# total numbers of each by timestep
to_plot <- aggregate(observed_only[, c("flu2", "rsv")], by=list(observed_only$timestep), sum)
# use this if want to look at daily reports
to_plot_m1 <- melt(to_plot, id.vars = "Group.1")
ggplot(to_plot_m1, aes(x = Group.1, y = value, colour = variable)) +
  geom_line()
# calulate which month of the year it is (i.e. seasonality)
to_plot$month <- month(as.Date(to_plot$Group.1, origin = start_date))
to_plot <- aggregate(to_plot[, c("flu2", "rsv")], by = list(to_plot$month), sum)
to_plot_m2 <- melt(to_plot, id.vars = "Group.1")
# plot the month of year vs cases (i.e. seasonality)
ggplot(to_plot_m2, aes(x = Group.1, y = value, colour = variable)) +
  geom_line() + 
  labs(x = "MOnth of the year", y = "Total cases")


#### Assign the states ####
# calculate time since previous flu infection
main_db[rsv==1, time_diff_rsv := timestep-shift(timestep, 1L, type ="lag"), by=indid]
main_db[flu=="Pos", time_diff_flu := timestep-shift(timestep, 1L, type ="lag"), by=indid]
main_db[, state := "tbd"]
# for each person in the data
for(i in unique(main_db$indid)){
  # subset the relevant data
  temp <- main_db[indid == i,]
  # when do they test positive for flu or rsv.
  flu_time <- as.numeric(unlist(temp[flu == "Pos", "timestep"]))
  rsv_time <- as.numeric(unlist(temp[rsv == 1, "timestep"]))
  end_time <- as.numeric(temp[min(which(temp$time_diff_flu>12),
                                  which(temp$time_diff_rsv>12)),"timestep"])
  if(is.na(rsv_time[1])){rsv_time = Inf}
  if(is.na(flu_time[1])){flu_time = Inf}
  if(is.na(end_time)){end_time = Inf}
  no_info <- as.numeric(unlist(temp[is.na(flu), "timestep"]))
  
  temp[timestep < flu_time[1] & timestep < rsv_time[1], state := "SS"]
  temp[timestep > flu_time[1] & timestep < rsv_time[1], state := "SPR"]
  temp[timestep < flu_time[1] & timestep >rsv_time[1], state := "PRS"]
  temp[timestep > flu_time[1] & timestep > rsv_time[1], state := "RR"]
  temp[timestep %in% flu_time & timestep< rsv_time[1], state := "SI"]
  temp[timestep < flu_time[1] & timestep %in% rsv_time, state := "IS"]
  temp[timestep > flu_time[1] & timestep %in% rsv_time, state := "IPR"]
  temp[timestep %in% flu_time & timestep > rsv_time[1], state := "PRI"]
  temp[timestep %in% flu_time & timestep %in% rsv_time, state := "II"]
  temp[timestep >= end_time, state := NA]
  temp[timestep %in% no_info, state := NA] 
  
  main_db[indid == i,] <- temp
}

main_db$factor_id <- as.numeric(main_db$indid)

# Convert states to numeric
main_db[state == "SS", num_state := 1]
main_db[state == "IS", num_state := 2]
main_db[state == "PRS", num_state := 99]
main_db[state == "SI", num_state := 5]
main_db[state == "II", num_state := 6]
main_db[state == "PRI", num_state := 7]
main_db[state == "SPR", num_state := 101]
main_db[state == "IPR", num_state := 9]
main_db[state == "RR", num_state := 11]
main_db[is.na(state), num_state := NA]
#### Create transition matrix ####

Q_matrix <- matrix(0,nrow=11, ncol=11)
colnames(Q_matrix) <- rownames(Q_matrix) <- c("SS", "IS", "PS", "RS", 
                                              "SI", "II", "PRI", "SP",
                                              "IPR", "SR", "RR")
# lambda I
Q_matrix["SS", "SI"] <- 0.005
Q_matrix["RS", "PRI"] <- 0.005
Q_matrix["IS", "II"] <- 0.005
Q_matrix["PS", "PRI"] <- 0.005
# lambda R
Q_matrix["SS", "IS"] <- 0.005
Q_matrix["SR", "IPR"] <- 0.005
Q_matrix["SI", "II"] <- 0.005
Q_matrix["SP", "IPR"] <- 0.005
# gamma I
Q_matrix["SI", "SP"] <- 1/3.8
Q_matrix["II", "IPR"] <- 1/3.8
Q_matrix["PRI", "RR"] <- 1/3.8
# gamma R
Q_matrix["IS", "PS"] <- 1/9
Q_matrix["II", "PRI"] <- 1/9
Q_matrix["IPR", "RR"] <- 1/9
# rho 
Q_matrix["PS", "RS"] <- 0.2
Q_matrix["SP", "SR"] <- 0.2

##### Create the misclassifciation matrix #####

E_matrix <- matrix(0.01,nrow=11, ncol=11)
colnames(E_matrix) <- rownames(E_matrix) <- c("SS", "IS", "PS", "RS", 
                                              "SI", "II", "PRI", "SP",
                                              "IPR", "SR", "RR")
# if you observe an II, it has to be II
E_matrix[,"II"] <- 0
# if you observe an IX, it can only be other IX's
E_matrix[,"IS"] <- 0
E_matrix[c("II", "IPR"),"IS"] <- 0.01
E_matrix[,"SI"] <- 0
E_matrix[c("II", "PRI"),"SI"] <- 0.01
E_matrix[,"IPR"] <- 0
E_matrix[c("II", "IS"),"IPR"] <- 0.01
E_matrix[,"PRI"] <- 0
E_matrix[c("II", "SI"),"PRI"] <- 0.01
# E_matrix["SS", c("IS","SI","II")] <- 0.2
# 
# E_matrix["SI", c("II")] <- 0.2
# E_matrix["IS", c("II")] <- 0.2
# 
# E_matrix["SP", c("SI", "IPR")] <- 0.2
# E_matrix["SR", c("SI", "IPR")] <- 0.2
# 
# E_matrix["PS", c("IS", "PRI")] <- 0.2
# E_matrix["RS", c("IS", "PRI")] <- 0.2
# 
# E_matrix["RR", c("PRI", "IPR", "II")] <- 0.2
# 
# E_matrix["IPR", c("II")] <- 0.2
# E_matrix["PRI", c("II")] <- 0.2

# convert to right way round
#E_matrix <- t(E_matrix)

#### times for pci splits ####

pci_splits <- c(30,60,90,120,150,259,289,319,349,379,409,439,469,499,529,623,653,683,713,743,773,803,833,863)


# parameters to keep constant
# fix 3 and 5 (gammas)
a <- c(3,5)
# then each transition (of which there are 16) for each time period
# then the other parameters (sigma covariate etc.) all need to be variable
b <- c(8,9,10)+1
c <- c(8,9,10)+1
for(i in 1:(length(pci_splits)-1)){
  b<- c(b,c+(i*5))
}
fixed_pars <- c(a,b)


##### RUN THE MODEL #####


observed_only <- main_db[!is.na(num_state)]
observed_only$sigma <- 1
#run the model
output <- msm(formula = num_state~timestep,
              subject = indid, 
              data = observed_only,
              qmatrix = Q_matrix,
              obstype = 1, 
             censor = c(99,101),
              censor.states = list(
               c(3,4),
              c(8,10)),
               fixedpars = fixed_pars,
              qconstraint = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
              ematrix = E_matrix,
             econstraint = rep(1,68),
             initprobs = c(rep(0.5,11)),
             est.initprobs = T,
           covariates = list("2-6" = ~sigma, "3-7" = ~sigma,
           "5-6" = ~sigma, "8-9" = ~sigma),
            constraint = list(sigma = c(1,1,1,1),
                              "timeperiod[30,60)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[60,90)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[90,120)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[120,150)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[150,259)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[259,289)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[289,319)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[319,349)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[349,379)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[379,409)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[409,439)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[439,469)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[469,499)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[499,529)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[529,623)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[623,653)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[653,683)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[683,713)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[713,743)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[743,773)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[773,803)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[803,833)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[833,863)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                              "timeperiod[863,Inf)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1)
                  ),
              pci = pci_splits,
                opt.method="nlm",
          hessian=FALSE, iterlim = 10
              # control=list(fnscale = 5000, reltol= 1e-16, 
              #                                  maxit = 100000
        
             )

save(output, file= "msm_output.Rdata")
