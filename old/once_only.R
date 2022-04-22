 # with the assumption that can only get each infection once per season


main_17$timestamp <- rep(1:82,577)

for(i in unique(main_17$indid)){

temp <- main_17[indid == i,]

flu_time <- as.numeric(unlist(temp[flu == "Pos", "timestamp"]))
rsv_time <- as.numeric(unlist(temp[rsv == 1, "timestamp"]))
if(is.na(rsv_time[1])){rsv_time = Inf}
if(is.na(flu_time[1])){flu_time = Inf}

temp[timestamp < flu_time[1] & timestamp < rsv_time[1], state := "SS"]
temp[timestamp > flu_time[1] & timestamp < rsv_time[1], state := "SPR"]
temp[timestamp < flu_time[1] & timestamp >rsv_time[1], state := "PRS"]
temp[timestamp > flu_time[1] & timestamp > rsv_time[1], state := "RR"]
temp[timestamp == flu_time[1] & timestamp< rsv_time[1], state := "SI"]
temp[timestamp < flu_time[1] & timestamp == rsv_time[1], state := "IS"]
temp[timestamp > flu_time[1] & timestamp == rsv_time[1], state := "IPR"]
temp[timestamp == flu_time[1] & timestamp > rsv_time[1], state := "PRI"]
temp[timestamp == flu_time[1] & timestamp == rsv_time[1], state := "II"]
temp[timestamp > tail(flu_time,1) & timestamp > tail(rsv_time,1), state := "RR"]

main_17[indid == i,] <- temp
}

main_17[state == "SS", num_state :=1]
main_17[state == "IS", num_state :=2]
main_17[state == "PRS", num_state :=99]
main_17[state == "SI", num_state :=5]
main_17[state == "II", num_state :=6]
main_17[state == "PRI", num_state :=7]
main_17[state == "SPR", num_state :=101]
main_17[state == "IPR", num_state :=9]
main_17[state == "RR", num_state :=11]


Q_matrix <- matrix(0,nrow=11, ncol=11)
colnames(Q_matrix) <- rownames(Q_matrix) <- c("SS", "IS", "PS", "RS", 
                                              "SI", "II", "PRI", "SP",
                                              "IPR", "SR", "RR")

# lambda I
Q_matrix["SS", "SI"] <- 0.05
Q_matrix["RS", "PRI"] <- 0.05
# lambda R
Q_matrix["SS", "IS"] <- 0.05
Q_matrix["SR", "IPR"] <- 0.05
# sigma I
Q_matrix["IS", "II"] <- 0.01
Q_matrix["PS", "PRI"] <- 0.01
# sigma R
Q_matrix["SI", "II"] <- 0.01
Q_matrix["SP", "IPR"] <- 0.01
# gamma I
Q_matrix["SI", "SP"] <- 0.5
Q_matrix["II", "IPR"] <- 0.5
Q_matrix["PRI", "RR"] <- 0.5
# gamma R
Q_matrix["IS", "PS"] <- 0.5
Q_matrix["II", "PRI"] <- 0.5
Q_matrix["IPR", "RR"] <- 0.5
# rho 
Q_matrix["PS", "RS"] <- 0.2
Q_matrix["SP", "SR"] <- 0.2


output <- msm(formula = num_state~npsdatecol,
              subject = indid, 
              data = main_17,
              qmatrix = Q_matrix,
              obstype = 1, 
              censor = c(99,101),
              censor.states = list(
                c(3,4),
                c(8,10)),
              qconstraint = c(1,2,3,4,5,4,2,6,7,3,7,7,6,8,3,1),
              opt.method="optim", control=list(fnscale = 5000, reltol= 1e-16, 
                                               maxit = 100000))


