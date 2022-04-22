# Evaluate the output
library(ggplot2)
# see the overall output 
output

# get the transition rates / hazard ratios 

# what is the impact of interaction?
# strength
# 
when_inf_first <- qratio.msm(output, ind1 = c(5,6), ind2 = c(1,2))
when_rsv_first <- qratio.msm(output, ind1 = c(2,6), ind2 = c(1,5))
print(paste0("Transmission from S to I with RSV is ", when_inf_first[1]," (", when_inf_first[3],
             " - ", when_inf_first[4],
             ") more likely when already infected with flu"))
print(paste0("Transmission from S to I with flu is ", when_rsv_first[1]," (", when_rsv_first[3],
             " - ", when_rsv_first[4],
             ") more likely when already infected with rsv"))
# duration



# What is the seasonality. Can plot over time? 
# hazard.msm creates a list of the hazard rates for each pci value
# For each time there are two parameters estimated
# extract from model
timing_hazards_list <- hazard.msm(output)
# combine into data table
timing_hazards <- data.table(from = pci_splits, 
                             to = c(pci_splits[-1],NA))
timing_hazards$flu_foi <- sapply(timing_hazards_list, `[[`, 2)
timing_hazards$rsv_foi <- sapply(timing_hazards_list, `[[`, 1)
timing_hazards[,Date := as.Date(from, origin = "2016-05-02")]
# make long
timing_hazards_m <- melt.data.table(timing_hazards, 
                                    id.vars = c("from", "to", "Date"), 
                                    measure.vars = c("flu_foi", "rsv_foi"))
# add cis
timing_hazards_m$lower <- c(sapply(timing_hazards_list, function(x) x[[2,2]]),sapply(timing_hazards_list, function(x) x[[1,2]]))
timing_hazards_m$upper <- c(sapply(timing_hazards_list, function(x) x[[2,3]]),sapply(timing_hazards_list, function(x) x[[1,3]]))
# plot the output 
ggplot(timing_hazards_m, aes(x = Date, y = value, colour = variable)) + 
  geom_line()+
  geom_point() + 
  theme_linedraw() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), alpha = 0.3, colour =NA) + 
  labs(y= "Hazard ratio compared to first month of data")

# What is the misclassifciation

# goodness of model?
plot.prevalence.msm(output)




