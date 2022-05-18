## Data setup, all in one with waning
library(data.table)
library(ggplot2)
library(msm)


load("~/Documents/GitHub/markov_model/data_waning.Rdata")

# the Q matrix
Q_matrix <- matrix(0,nrow=16, ncol=16)
colnames(Q_matrix) <- rownames(Q_matrix) <- c("SS", "IS", "PS", "RS", 
                                              "SI", "II", "PI","RI",
                                              "SP", "IP", "PP", "RP",
                                              "SR", "IR", "PR" ,"RR")
# lambda I
Q_matrix["SS", "SI"] <- 0.001
Q_matrix["RS", "RI"] <- 0.01
Q_matrix["IS", "II"] <- 0.01
Q_matrix["PS", "PI"] <- 0.001
# lambda R
Q_matrix["SS", "IS"] <- 0.005
Q_matrix["SR", "IR"] <- 0.005
Q_matrix["SI", "II"] <- 0.006
Q_matrix["SP", "IP"] <- 0.006
# gamma I
Q_matrix["SI", "SP"] <- 1/3.8
Q_matrix["II", "IP"] <- 1/3.8
Q_matrix["PI", "PP"] <- 1/3.8
Q_matrix["RI", "RP"] <- 1/3.8
# gamma R
Q_matrix["IS", "PS"] <- 1/9
Q_matrix["II", "PI"] <- 1/9
Q_matrix["IP", "PP"] <- 1/9
Q_matrix["IR", "PR"] <- 1/9
# rho 
Q_matrix["PS", "RS"] <- 0.1
Q_matrix["PI", "RI"] <- 0.1
Q_matrix["PP", "RP"] <- 0.1
Q_matrix["PR", "RR"] <- 0.1
Q_matrix["SP", "SR"] <- 0.1
Q_matrix["IP", "IR"] <- 0.1
Q_matrix["PP", "PR"] <- 0.1
Q_matrix["RP", "RR"] <- 0.1

# waning R
Q_matrix["RS", "SS"] <- 1/365
Q_matrix["RI", "SI"] <- 1/365
Q_matrix["RP", "SP"] <- 1/365
Q_matrix["RR", "SR"] <- 1/365

# waning I
Q_matrix["SR", "SS"] <- 1/365
Q_matrix["IR", "IS"] <- 1/365
Q_matrix["PR", "PS"] <- 1/365
Q_matrix["RR", "RS"] <- 1/365

### THE E matrix (misclassifcation)
E_matrix <- matrix(0.01,nrow=16, ncol=16)
colnames(E_matrix) <- rownames(E_matrix) <- c("SS", "IS", "PS", "RS", 
                                              "SI", "II", "PI","RI",
                                              "SP", "IP", "PP", "RP",
                                              "SR", "IR", "PR" ,"RR")
# if you observe an II, it has to be II
E_matrix[,"II"] <- 0

# if you observe an IX, it can only be other IX's
E_matrix[,"IS"] <- 0
E_matrix[c("II", "IP", "IR"),"IS"] <- 0.01
E_matrix[,"SI"] <- 0
E_matrix[c("II", "PI", "RI"),"SI"] <- 0.01

E_matrix[,"IP"] <- 0
E_matrix[c("II", "IS", "IR"),"IP"] <- 0.01
E_matrix[,"PI"] <- 0
E_matrix[c("II", "SI", "RI"),"PI"] <- 0.01

E_matrix[,"IR"] <- 0
E_matrix[c("II", "IS", "IP"),"IR"] <- 0.01
E_matrix[,"RI"] <- 0
E_matrix[c("II", "SI", "PI"),"RI"] <- 0.01

# Parameters (based on order in Q matrix)
# 1 = FOI rsv
# 2 = FOI flu
# 3 = gamma R
# 4 = sigma foi I
# 5 = rho
# 6 = waning R
# 7 = sigma foi R
# 8 = gamma I
# 9 = waning I

# fixed should be (for now): gammaR (3), gammaI (8), rho (5), waningI (6), waningR (9)
fixed_pars <- c(3,5,6,8,9)

# which parameters in the matrix are the same?
q_constraints <- c(1,2,3,4,5,4,6,2,
                   7,8,3,8,5,8,6,8,
                   7,5,3,5,5,5,6,5,
                   1,9,3,9,5,9,6,9)
#run the model
output <- msm(formula = state~timestep,
              subject = factor_id, 
              data = input_data,
              qmatrix = Q_matrix,
              obstype = 1, 
              censor = c(101,102,103,104),
              censor.states = list(
                c(1,3,4,9,11,12,13,15,16), # 0 rsv, 0 flu
                c(2,10,14), # 1 rsv, 0 flu
                c(5,7,8), # 0 rsv, 1 flu
                c(6) # 1 rsv, 1 flu
                ),
              fixedpars = fixed_pars,
              qconstraint = q_constraints,
              ematrix = E_matrix,
              econstraint = rep(1,153),
              initprobs = c(1,rep(0.1,15)), # probability of starting in each state
              est.initprobs = T,
              covariates = list(  "1-2" = ~age_cat + rsv_week,
                                  "1-5" = ~age_cat + flu_week, 
                                  "2-6" = ~age_cat + flu_week,
                                  "3-7" = ~age_cat + flu_week,
                                  "4-8" = ~age_cat + flu_week,
                                  "5-6" = ~age_cat + rsv_week,
                                  "9-10" = ~age_cat + rsv_week,
                                  "13-14" = ~age_cat + rsv_week
                ),
              constraint = list(age_cat5_to_18 = c(1,2,2,2,2,1,1,1),
                                age_cat19_to_65 = c(1,2,2,2,2,1,1,1),
                                age_catover_65 = c(1,2,2,2,2,1,1,1),
                                rsv_week = c(1,1,1,1), 
                                flu_week = c(1,1,1,1)
                                ) ,
              method = "BFGS", 
              control=list(maxit=10000)
              
)

save(output, file= "msm_output.Rdata")
