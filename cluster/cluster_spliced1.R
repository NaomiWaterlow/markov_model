
library(msm)
library(data.table)
library(minqa)

load("main_db_seasons.Rdata")

#### times for pci splits ####

pci_splits <- c(30,60,90,120,150,181,
                211,241,271,301,331,361,391,421,451,467,
                497,527,557,587,617,647,677,707,737)
# pci_splits <- seq(from = 30, by = 30, to = 289)

# parameters to keep constant
# fix 3 and 5 (gammas)
a <- c(3,7)
# then each transition (of which there are 16) for each time period
# c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1) <- 1s and 2s variable, the rest fixed
# then the other parameters (sigma covariate etc.) all need to be variable
b <- c(10,11,12)
c <- c(10,11,12)
for(i in 1:(length(pci_splits)-1)){
  b<- c(b,c+(i*5))
}

fixed_pars <- c(a,b)

# fixed_pars <- c(fixed_pars, 500)
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
              qconstraint = c(1,2,3,4,5,4,2,6,7,3,7,7,6,5,3,1), # different FOIs for if sigma present or not
              ematrix = E_matrix,
              econstraint = rep(1,68),
              initprobs = c(1,rep(0.01,10)),
              est.initprobs = T,
              constraint = list(
                # same impact of time on the paired FOIs with and without sigma
                "timeperiod[30,60)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[60,90)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[90,120)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[120,150)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[150,181)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[181,211)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[211,241)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[241,271)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1),
                "timeperiod[271,301)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[301,331)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[331,361)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[361,391)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[391,421)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[421,451)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[451,467)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[467,497)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[497,527)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[527,557)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[557,587)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[587,617)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[617,647)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[647,677)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[677,707)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[707,737)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1), 
                "timeperiod[737,Inf)" = c(1,2,3,2,4,2,2,1,5,3,5,5,1,4,3,1)
                
              ),
              pci = pci_splits,
              control=list(fnscale = 5000)
              
)
save(output, file= "msm_output.Rdata")
