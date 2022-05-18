# Evaluate the output
library(ggplot2)
library(msm)
library(data.table)
# see the overall output 
output

# get the transition rates / hazard ratios 

# what is the impact of interaction?
# strength
# 
when_inf_first <- qratio.msm(output, ind1 = c(5,6), ind2 = c(1,2))
when_rsv_first <- qratio.msm(output, ind1 = c(2,6), ind2 = c(1,5))
print(paste0("Transmission from S to I with RSV is ", round(when_inf_first[1],2)," (", round(when_inf_first[3],2),
             " - ", round(when_inf_first[4],2),
             ") more likely when already infected with flu"))
print(paste0("Transmission from S to I with flu is ", round(when_rsv_first[1],2)," (", round(when_rsv_first[3],2),
             " - ", round(when_rsv_first[4],2),
             ") more likely when already infected with rsv"))



# What is the seasonality. Can plot over time? 
# hazard.msm creates a list of the hazard rates for each pci value
# For each time there are two parameters estimated
# extract from model
timing_hazards_list <- hazard.msm(output)
# combine into data table
timing_hazards <- data.table(from = pci_splits, 
                             to = c(pci_splits[-1],NA))
timing_hazards$Influenza <- sapply(timing_hazards_list, `[[`, 2)
timing_hazards$RSV <- sapply(timing_hazards_list, `[[`, 1)
timing_hazards[,Date := as.Date(from, origin = "2016-05-02")]
# make long
timing_hazards_m <- melt.data.table(timing_hazards, 
                                    id.vars = c("from", "to", "Date"), 
                                    measure.vars = c("Influenza", "RSV"))
# add cis
timing_hazards_m$lower <- c(sapply(timing_hazards_list, function(x) x[[2,2]]),sapply(timing_hazards_list, function(x) x[[1,2]]))
timing_hazards_m$upper <- c(sapply(timing_hazards_list, function(x) x[[2,3]]),sapply(timing_hazards_list, function(x) x[[1,3]]))
# plot the output 
ggplot(timing_hazards_m, aes(x = Date, y = value, colour = variable)) + 
  geom_line()+
  geom_point() + 
  theme_linedraw() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), alpha = 0.3, colour =NA) + 
  labs(y= "Hazard ratio compared to first month of data", fill = "Force of Infection", 
       colour = "Force of Infection")

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
        legend.text = element_text(size=15))



BLANK <- ggplot() + theme_minimal() + labs(title ="C") +
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
