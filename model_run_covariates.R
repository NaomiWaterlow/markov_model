
library(msm)
library(data.table)
library(minqa)

load("main_db_seasons_april_covariates.Rdata")

# run this line to run the sensitivity analysis
main_db<- main_db[track_colour %in% c("second", "third")]

#### times for pci splits ####

#pci_splits <-  seq(from = 30, by = 30, to = 289)
pci_splits <- c(30,60,90,120,150,181,
                211,241,271,301,331,361,391,421,451,467,
                497,527,557,587,617,647,677,707,737)

# run this for pci splits
pci_splits <- c(
                211,241,271,301,331,361,391,421,451,467,
                497,527,557,587,617,647,677,707,737)

# parameters to keep constant
# fix 3 and 7 (gammas)
a <- c(3,7) # 5 includes rho
# then each transition (of which there are 16) for each time period
# c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1) <- 1s and 2s variable, the rest fixed
# then the other parameters (sigma covariate etc.) all need to be variable
# b <- c(10,11,12)
# c <- c(10,11,12)
# for(i in 1:(length(pci_splits)-1)){
#   b<- c(b,c+(i*5))
# }
fixed_pars <- c(a)


#### Create transition matrix ####

Q_matrix <- matrix(0,nrow=11, ncol=11)
colnames(Q_matrix) <- rownames(Q_matrix) <- c("SS", "IS", "PS", "RS",
                                              "SI", "II", "PRI", "SP",
                                              "IPR", "SR", "RR")
# lambda I
Q_matrix["SS", "SI"] <- 0.01
Q_matrix["RS", "PRI"] <- 0.01
Q_matrix["IS", "II"] <- 0.001
Q_matrix["PS", "PRI"] <- 0.001
# lambda R
Q_matrix["SS", "IS"] <- 0.01
Q_matrix["SR", "IPR"] <- 0.01
Q_matrix["SI", "II"] <- 0.001
Q_matrix["SP", "IPR"] <- 0.001
# gamma I
Q_matrix["SI", "SP"] <- 1/3.8
Q_matrix["II", "IPR"] <- 1/3.8
Q_matrix["PRI", "RR"] <- 1/3.8
# gamma R
Q_matrix["IS", "PS"] <- 1/9
Q_matrix["II", "PRI"] <- 1/9
Q_matrix["IPR", "RR"] <- 1/9
# rho 
Q_matrix["PS", "RS"] <- 1/30
Q_matrix["SP", "SR"] <- 1/30

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


#run the model
output <- msm(formula = num_state~timestep,
              subject = factor_id,
              data = main_db,
              qmatrix = Q_matrix,
              obstype = 1,
              censor = c(99,101),
              censor.states = list(
                c(3,4),
                c(8,10)),
              fixedpars = fixed_pars,
              qconstraint = c(1,2,3,4,5,4,2,6,7,3,7,7,6,5,3,1),
              ematrix = E_matrix,
              econstraint = rep(1,68),
              initprobs = c(0.95,rep(0.01,10)),
              est.initprobs = T,
              covariates = list(  "1-2" = ~age_grp + timeperiodRSV,
                                  "1-5" = ~age_grp + flutype + timeperiodFlu,
                                  "2-6" = ~age_grp + flutype + timeperiodFlu,
                                  "3-7" = ~age_grp + flutype + timeperiodFlu,
                                  "4-7" = ~age_grp + flutype + timeperiodFlu,
                                  "5-6" = ~age_grp + timeperiodRSV,
                                  "8-9" = ~age_grp + timeperiodRSV,
                                  "10-9" = ~age_grp + timeperiodRSV
              ),
              constraint = list(
                age_grpage2 = c(1,2,2,2,2,1,1,1),
                                age_grpage3 = c(1,2,2,2,2,1,1,1),
                                age_grpage4 = c(1,2,2,2,2,1,1,1),
                                flutypeB = c(1,1,1,1),
                timeperiodRSVt1  = c(1,1,1,1),
                timeperiodRSVt2  = c(1,1,1,1),
                timeperiodRSVt3  = c(1,1,1,1),
                timeperiodRSVt4  = c(1,1,1,1),
                timeperiodRSVt5  = c(1,1,1,1),
                timeperiodRSVt6  = c(1,1,1,1),
                timeperiodRSVt7  = c(1,1,1,1),
                timeperiodRSVt8  = c(1,1,1,1),
                timeperiodRSVt9  = c(1,1,1,1),
                timeperiodRSVt10 = c(1,1,1,1),
                timeperiodRSVt11 = c(1,1,1,1),
                timeperiodRSVt12 = c(1,1,1,1),
                timeperiodRSVt13 = c(1,1,1,1),
                timeperiodRSVt14 = c(1,1,1,1),
                timeperiodRSVt15 = c(1,1,1,1),
                timeperiodRSVt16 = c(1,1,1,1),
                timeperiodRSVt17 = c(1,1,1,1),
                timeperiodRSVt18 = c(1,1,1,1),
                timeperiodRSVt19 = c(1,1,1,1),
                timeperiodRSVt20 = c(1,1,1,1),
                timeperiodRSVt21 = c(1,1,1,1),
                timeperiodRSVt22 = c(1,1,1,1),
                timeperiodRSVt23 = c(1,1,1,1),
                timeperiodRSVt24 = c(1,1,1,1),
                timeperiodRSVt25 = c(1,1,1,1),
                timeperiodFlut1  = c(1,1,1,1),
                timeperiodFlut2  = c(1,1,1,1),
                timeperiodFlut3  = c(1,1,1,1),
                timeperiodFlut4  = c(1,1,1,1),
                timeperiodFlut5  = c(1,1,1,1),
                timeperiodFlut6  = c(1,1,1,1),
                timeperiodFlut7  = c(1,1,1,1),
                timeperiodFlut8  = c(1,1,1,1),
                timeperiodFlut9  = c(1,1,1,1),
                timeperiodFlut10 = c(1,1,1,1),
                timeperiodFlut11 = c(1,1,1,1),
                timeperiodFlut12 = c(1,1,1,1),
                timeperiodFlut13 = c(1,1,1,1),
                timeperiodFlut14 = c(1,1,1,1),
                timeperiodFlut15 = c(1,1,1,1),
                timeperiodFlut16 = c(1,1,1,1),
                timeperiodFlut17 = c(1,1,1,1),
                timeperiodFlut18 = c(1,1,1,1),
                timeperiodFlut19 = c(1,1,1,1),
                timeperiodFlut20 = c(1,1,1,1),
                timeperiodFlut21 = c(1,1,1,1),
                timeperiodFlut22 = c(1,1,1,1),
                timeperiodFlut23 = c(1,1,1,1),
                timeperiodFlut24 = c(1,1,1,1),
                timeperiodFlut25 = c(1,1,1,1)
    
                  # old saved c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
              ),
            # pci = pci_splits
               method= "BFGS"
              ,control=list(fnscale = 57500, maxit = 5)
)

save(output, file= "sesitivity.Rdata")
