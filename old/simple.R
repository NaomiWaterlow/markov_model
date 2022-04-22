# simple SIR

main_17

flu_only <- main_17[, c("indid", "npsdatecol", "flu", "state")]
flu_only[state == "XX" | state == "XI", state := 99]
flu_only[state =="IX", state := 2]

flu_only$timestamp <- rep(1:82,577)

for(i in unique(flu_only$indid)){
  
  temp <- flu_only[indid == i,]
  
  flu_time <- as.numeric(unlist(temp[flu == "Pos", "timestamp"]))

  if(is.na(flu_time[1])){flu_time = Inf}
  
  temp[flu == "Neg" , state := 99]
  temp[flu == "Pos", state := 2]
  temp[is.na(flu), state := NA]
  
  flu_only[indid == i,] <- temp
}

start_date <- min(flu_only$npsdatecol)

flu_only$state <- as.numeric(flu_only$state)
flu_only[,timestep := (npsdatecol - start_date)]
flu_only[,timestep := as.numeric(timestep, units="days")]

Q_matrix <- rbind(c(0,0.02,0),
                  c(0,0,1/3.8), 
                  c(0.8,0,0))
# rows represent underlying states, columns represent observed states
E_matrix <- rbind(c(0,0,0), 
                  c(0,0,0.1), 
                  c(0,0,0))

colnames(Q_matrix) <- rownames(Q_matrix) <- c("S", "I", "R")
colnames(E_matrix) <- rownames(E_matrix) <- c("S", "I", "R")

output <- msm(formula = state~timestep,
              subject = indid, 
              data = flu_only,
              qmatrix = Q_matrix,
              ematrix = E_matrix,
              obstype = 1, 
              censor = c(99),
              censor.states = list(
                c(1,3)),
              opt.method="bobyqa", 
              fixedpars = c(2,3,5),
           #  covariates = list("1-2" = timeperiod),
              pci = c(100))



to_plot <- flu_only[flu=="Pos", .N,  by=npsdatecol ]

ggplot(to_plot, aes(x = npsdatecol, y=N )) + 
  geom_point()
