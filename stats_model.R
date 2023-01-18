# self controlled case study

# Data setup, splicing seasons together
library(data.table)
library(ggplot2)
PHIRST_NP <- data.table(read.csv("~/Documents/OneDrive - London School of Hygiene and Tropical Medicine/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Nasopharyngeal specimens 2016-2018 Flu & RSV 2022-04-19.csv",
                                 na.strings = ""))
# format correctly
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

main_db[, factor_id := as.numeric(factor(indid))]
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

# remove unneccesary columns
main_db <- main_db[,c("npsdatecol","factor_id","rsv", "flu", "timestep", "flutype", 
                      "age_grp")]

num_days_neg <- 14
flu_subset <- main_db[flu == 1]
num_flu_infected <- length(unique(flu_subset$factor_id))
flu_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =factor_id]
flu_subset[, first_swab := F]

flu_subset[time_diff>num_days_neg | is.na(time_diff), first_swab := T]

main_db[,first_swab_flu := F]
main_db[flu_subset, on=c("npsdatecol", "factor_id"), first_swab_flu := i.first_swab]

rsv_subset <- main_db[rsv == 1]
num_rsv_infected <- length(unique(rsv_subset$factor_id))
rsv_subset[, time_diff := timestep-shift(timestep, 1L, type ="lag"), by =factor_id]
rsv_subset[, first_swab := F]

rsv_subset[time_diff>num_days_neg | is.na(time_diff), first_swab := T]

main_db[,first_swab_rsv := F]
main_db[rsv_subset, on=c("npsdatecol", "factor_id"), first_swab_rsv := i.first_swab]
# label all the exposure groups

exposure_period_1 <- 7
exposure_period_2 <- 14
exposure_period_3 <- 21
main_db$exposure_group <- "tbd"

# subset to ignore those who don't have both a flu and an rsv infection
 # - done in the loop below

# for the influenza background rate, just need to take account of time, not the number of flu cases
# as just putting that as an exposure grouping. 


for(i in unique(main_db$factor_id)){
  # subset for each person
  temp <- main_db[factor_id == i,]
  # calculate timings of infection
  flu_time <- as.numeric(unlist(temp[first_swab_flu == T, "timestep"]))
  rsv_time <- as.numeric(unlist(temp[first_swab_rsv == T, "timestep"]))
  if(is.na(rsv_time[1])){rsv_time = Inf}
  if(is.na(flu_time[1])){flu_time = Inf}
  if(any(flu_time== Inf) | any(rsv_time == Inf)){
    temp$exposure_group <- NA

  }else{
# now need to work out the exposure periods
  # here taking flu to be the event and rsv to be the exposure
    # so need to only look at the rsv exposure periods. after first exposure!
    temp[timestep < rsv_time[1], exposure_group := NA ]
    temp[timestep >= rsv_time[1] & timestep <= (rsv_time[1]+exposure_period_1), exposure_group := 1 ]
    temp[timestep > (rsv_time[1]+exposure_period_1) & timestep <= (rsv_time[1]+exposure_period_2), exposure_group := 2 ]
    temp[timestep > (rsv_time[1]+exposure_period_2) & timestep <= (rsv_time[1]+exposure_period_3), exposure_group := 3 ]
    temp[timestep >= rsv_time[2] & timestep <= (rsv_time[2]+exposure_period_1), exposure_group := 1 ]
    temp[timestep > (rsv_time[2]+exposure_period_1) & timestep <= (rsv_time[2]+exposure_period_2), exposure_group := 2 ]
    temp[timestep > (rsv_time[2]+exposure_period_2) & timestep <= (rsv_time[2]+exposure_period_3), exposure_group := 3 ]
    temp[timestep >= rsv_time[3] & timestep <= (rsv_time[3]+exposure_period_1), exposure_group := 1 ]
    temp[timestep > (rsv_time[3]+exposure_period_1) & timestep <= (rsv_time[3]+exposure_period_2), exposure_group := 2 ]
    temp[timestep > (rsv_time[3]+exposure_period_2) & timestep <= (rsv_time[3]+exposure_period_3), exposure_group := 3 ]

}
    main_db[factor_id == i,] <- temp

}


included <- na.omit(main_db[exposure_group != "tbd"])

#Specify groupings
start1 <- yday(min(included$npsdatecol))
end1 <- yday(max(included[npsdatecol < as.Date("2017-01-01")]$npsdatecol))
start2 <- yday(min(included[npsdatecol > as.Date("2017-01-01")]$npsdatecol))
end2 <- yday(max(included[npsdatecol < as.Date("2018-01-01")]$npsdatecol))
start3 <- yday(min(included[npsdatecol > as.Date("2018-01-01")]$npsdatecol))
end3 <- yday(max(included[npsdatecol < as.Date("2019-01-01")]$npsdatecol))

time_zones <- c(seq(from = min(start1,start2,start3), to = max(end1,end2,end3), by = 14),Inf)

# do want to do it time of the year?? Yes will try this.. 
included[,yday := yday(npsdatecol)]
for(i in 1:length(time_zones)){
   included[yday >= time_zones[i]& yday < time_zones[i+1], time_zone := paste0("t_",i) ]
  
}

time_zone_lookup <- data.table(
  time_zone = paste0("t_", c(1:length(time_zones))),
  parameter = 1
)

# overall incdence for this individual
calc_ll_individual <- function(temp, p_rsv_1,p_rsv_2, p_rsv_3, time_zone_lookup){
  # get the combinations of exposure groups
  intervals <- unique(temp[,c("exposure_group", "time_zone")])
  intervals[,interval := 1:nrow(intervals)]
  num_intervals <- nrow(intervals)
  temp[intervals, on = c("exposure_group", "time_zone"), interval := i.interval]
  # get the number of 
  output_sum <- temp[,sum(first_swab_flu), by =c("interval", "exposure_group", "time_zone")]
  # calculate the durations of the intervals
  output_sum$min <-  temp[,min(timestep), by =c("interval", "exposure_group", "time_zone")]$V1
  output_sum$max <-  temp[,max(timestep), by =c("interval", "exposure_group", "time_zone")]$V1
  output_sum[,duration := (1+max)-min]
  # match up the exposure groups with the parameters
  output_sum[ exposure_group == 1, exposure_param := p_rsv_1 ]
  output_sum[ exposure_group == 2, exposure_param := p_rsv_2 ]
  output_sum[ exposure_group == 3, exposure_param := p_rsv_3 ]
  output_sum[time_zone_lookup, on = c("time_zone"), time_param := i.parameter]
  # calculate the incidence in each interval
  output_sum[, incidence := duration * exposure_param* time_param]
  # total incidence
  individual_incidence <- sum(output_sum$incidence)
  # probability in each interval
  output_sum[,prob_interval := incidence/individual_incidence]
  # convert to log probability
  output_sum[,log_prob := log(prob_interval)]
  # multinomial contribution - would have be ^V1, but on log scale so *
  output_sum[, multinom_contribution := V1*log_prob]
  
  log_multinom_indiviual <- sum(output_sum$multinom_contribution)
  
  return(log_multinom_indiviual)
}


calc_ll_model <- function(input_params, input_data){
  # start off the likelihood
ll_all <- 0 

# convert input parameters to a table for later lookup
p_rsv_1 <- input_params[1]
p_rsv_2 <- input_params[2]
p_rsv_3 <- input_params[3]

time_zone_lookup$parameter <- input_params[4:25]

for(i in unique(input_data$factor_id)){
  # subset for each person
  temp <- input_data[factor_id == i,]
  # calculate the individual likelihoods
  ll_i <- calc_ll_individual(temp = temp, p_rsv_1 =p_rsv_1, 
                                  p_rsv_2 =p_rsv_2, p_rsv_3 =p_rsv_3,
                                  time_zone_lookup = time_zone_lookup )

# add
  ll_all <- ll_all + ll_i
}

return(ll_all)
}


calc_ll_all <- function(input_params){
  
  ll_model <- calc_ll_model(input_params = input_params,input_data=included)
  ll_prior <- 0#calc_ll_prior()
  
  ll_posterior <- ll_prior + ll_model
  
  return(ll_posterior)
  
}


init_theta <- rep(1,25)
names(init_theta) <- c("p_rsv_1", "p_rsv_2", "p_rsv_3", paste0("t_", c(1:length(time_zones))))
proposal_sd <- rep(0.1,25)
names(proposal_sd) <- c("p_rsv_1", "p_rsv_2", "p_rsv_3", paste0("t_", c(1:length(time_zones))))
lower_limits <- rep(0,25)
upper_limits <- rep(5,25)
names(lower_limits) = names(upper_limits) = names(proposal_sd)

trace <- mcmcMH(target = calc_ll_all, # target distribution
               init.theta = init_theta , # intial parameter guess
               proposal.sd = proposal_sd, 
               n.iterations = 10000, 
               print.info.every = 100,# standard deviation of
               limits = list(lower=lower_limits ,
                             upper=upper_limits),

)
  
trace_out <- data.table(trace$trace)
trace_out$timestep <- 1:nrow(trace_out)
trace_out_m <- melt(trace_out, id.vars = "timestep")

params <- unique(trace_out_m$variable)

ggplot(trace_out_m[variable %in% params[c(20:25,26)]], aes(x = timestep, y = value) ) + facet_grid(variable~., 
                                                                scales = "free_y") + 
  geom_line()

output_estimate <- optim(par = c(0.2,0.2,0.2,rep(1,22)),
      fn = calculate_ll_all,
      input_data = included, 
      hessian=T)

fisher_info<-solve(output_estimate$hessian)
prop_sigma<-sqrt(diag(fisher_info))
upper<-output_estimate$par+1.96*prop_sigma
lower<-output_estimate$par-1.96*prop_sigma
interval<-unlist(data.frame(value=output_estimate$par, upper=upper, lower=lower))

output_samples[samp+(j-1)*n_samps,] <- 
  c(samp, output_estimate$value, rhos[j], interval[1:3])
print(output_samples)
  





library(beepr)
beep(4)

output_samples <- data.table(output_samples)
together <- output_samples[,mean(estimate), by = "rho"]
together$log_lik <- output_samples[,mean(log_lik), by = "rho"]$V1
colnames(together)[2] <- "interaction"

