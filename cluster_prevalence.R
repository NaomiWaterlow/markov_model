library(msm)
library(data.table)

load("july_fixed_7_rep.Rdata")

prev3 <- prevalence.msm(x = output, 
                        timezero = 468,
                        times = seq(468, 757,7),
                        piecewise.times = pci_splits[which(pci_splits>497)], 
                        piecewise.covariates = list(timeperiodRSV = c("t16","t17", "t18", 
                                                                      "t19", "t20", "t21", "t22", 
                                                                      "t23", "t24", "t25"),
                                                    timeperiodFlu = c("t16","t17", "t18", 
                                                                      "t19", "t20", "t21", "t22", 
                                                                      "t23", "t24", "t25")),
                        ci = "normal")

save(prev3, file= "prev3.Rdata")

prev2 <- prevalence.msm(x = output, 
                        timezero = 182,
                        times= seq(182, 467, 7),
                        piecewise.times = c(211,241,271,301,331,361,391,421,451), 
                        piecewise.covariates = list(timeperiodRSV = c("t6","t7", "t8", 
                                                                      "t9", "t10", "t11", "t12", 
                                                                      "t13", "t14", "t15"),
                                                    timeperiodFlu = c("t6","t7", "t8", 
                                                                      "t9", "t10", "t11", "t12", 
                                                                      "t13", "t14", "t15")),
                        ci = "normal")

save(prev2, file= "prev2.Rdata")

prev1 <- prevalence.msm(x = output, 
                        timezero = 0,
                        times = seq(0,181,7),
                        piecewise.times = pci_splits[which(pci_splits<181 )], 
                        piecewise.covariates = list(timeperiodRSV = c("t0","t1","t2", "t3", 
                                                                      "t4", "t5"),
                                                    timeperiodFlu = c("t0","t1","t2", "t3", 
                                                                      "t4", "t5")),
                        ci = "normal")

save(prev1, file= "prev1.Rdata")