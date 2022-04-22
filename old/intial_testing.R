# install.packages("msm")
library(msm)
library(data.table)
library(ggplot2)
library(minqa)

PHIRST_Master<- data.table(read.csv("~/Documents/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Master 2016-2018 Flu & RSV 2021-10-19.csv"))
View(PHIRST_Master)

PHIRST_NP <- data.table(read.csv("~/Documents/Data/PHIRST Nasopharyngeal specimens and Master 2021-10-19/PHIRST Nasopharyngeal specimens 2016-2018 Flu & RSV 2021-10-19.csv", 
                                 na.strings = c(""))) 
View(PHIRST_NP)

main_db <- PHIRST_NP[,c("indid", "npsdatecol", "flu", "npsrsva", "npsrsvb")]
main_db[npsrsva == 1 | npsrsvb == 1, rsv := 1]
main_db$indid <- as.factor(main_db$indid)
main_db$npsdatecol <- as.Date(main_db$npsdatecol, "%d-%m-%Y")

main_16 <- copy(main_db[npsdatecol < as.Date("2017-01-01")])
main_17 <- copy(main_db[npsdatecol < as.Date("2018-01-01") & npsdatecol > as.Date("2017-01-01")])
main_18 <- copy(main_db[ npsdatecol > as.Date("2018-01-01")])

##### just working with 2017 for now #####

# check the maximum time between samples
main_17[, diff := npsdatecol-shift(npsdatecol, 1L, type= "lag"), by =indid]
max(main_17$diff,na.rm=T)
# 6 days- can I assume that noone got infected and cleared within the 6 day swabs? assume so for now

main_17[(flu != "Pos"| is.na(rsv)) & (rsv !=1 | is.na(rsv)) , state := "XX"]
main_17[(flu == "Pos") & (rsv !=1| is.na(rsv)), state := "IX"]
main_17[(flu != "Pos"| is.na(rsv)) & (rsv ==1), state := "XI"]
main_17[(flu == "Pos") & (rsv ==1), state := "XX"]

main_17[state == "XX",num_state := 100]
main_17[state == "XI",num_state := 1000]
main_17[state == "IX",num_state := 2000]
main_17[state == "II",num_state := 6]

main_17[PHIRST_Master, on = c(indid="ind_id"),  ]

# where are transitions allowed?
Q_matrix <- rbind(c(0,0.02,0,0,0.02,0,0,0,0,0,0,0,0,0,0,0),#0.2
                  c(0,0,0.2,0,0,0.01,0,0,0,0,0,0,0,0,0,0),#2
                  c(0,0,0,0.2,0,0,0.01,0,0,0,0,0,0,0,0,0),#3
                  c(0.05,0,0,0,0,0,0,0.02,0,0,0,0,0,0,0,0),#4
                  c(0,0,0,0,0,0.01,0,0,0.2,0,0,0,0,0,0,0),#5
                  c(0,0,0,0,0,0,0.2,0,0,0.2,0,0,0,0,0,0),#6
                  c(0,0,0,0,0,0,0,0.2,0,0,0.2,0,0,0,0,0),#7
                  c(0,0,0,0,0.05,0,0,0,0,0,0,0.2,0,0,0,0),#8
                  c(0,0,0,0,0,0,0,0,0,0.01,0,0,0.2,0,0,0),#9
                  c(0,0,0,0,0,0,0,0,0,0,0.2,0,0,0.2,0,0),#10
                  c(0,0,0,0,0,0,0,0,0,0,0,0.2,0,0,0.2,0),#11
                  c(0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0.2),#12
                  c(0.05,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0,0),#13
                  c(0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0),#14
                  c(0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0.2),#15
                  c(0,0,0,0.05,0,0,0,0,0,0,0,0.05,0,0,0,0)#16
       )


output <- msm(formula = num_state~npsdatecol,
    subject = indid, 
    data = main_17,
    qmatrix = Q_matrix,
    # ematrix = E_matrix,00
    obstype = 1, 
    censor = c(100,1000,2000),
    censor.states = list(
      c(1,3,4,9,11,12,13,15,16),
      c(2,10,14), 
      c(5,7,8)),
    qconstraint = c(1,2,3,4,5,4,6,2,7,8,3,8,5,8,6,8,
                    7,5,3,5,5,5,6,5,1,9,3,9,5,9,6,9),
  opt.method="bobyqa")


# E_matrix <- rbind(c(1,0,0,0), #SS
#                   c(0,1,0,0), #SI
#                   c(1,0,0,0), #SP
#                   c(1,0,0,0), #SR
#                   c(0,0,1,0), #IS
#                   c(0,0,0,1), #II
#                   c(0,0,1,0), #IP
#                   c(0,0,1,0), #IR
#                   c(1,0,0,0), #PS
#                   c(0,1,0,0), #PI
#                   c(1,0,0,0), #PP
#                   c(1,0,0,0), #PR
#                   c(1,0,0,0), #RS
#                   c(0,1,0,0), #RI
#                   c(1,0,0,0), #RP
#                   c(1,0,0,0) #RR
#                   )
colnames(E_matrix) <- c("XX", "XI", "IX", "II")


####### NO LONGER USING #####


# then need to convert them all into different states
temp <- main_17[indid == "A100-005",]
temp$timestamp <- c(1:nrow(temp))
flu_time <- as.numeric(unlist(temp[flu == "Pos", "timestamp"]))
rsv_time <- as.numeric(temp[rsv == 1, "timestamp"])

if(is.na(rsv_time[1])){rsv_time = Inf}
if(is.na(flu_time[1])){flu_time = Inf}
temp[timestamp == flu_time[1] & timestamp< rsv_time[1], state := "SI"]
temp[timestamp < flu_time[1] & timestamp == rsv_time[1], state := "IS"]
temp[timestamp > flu_time[1] & timestamp < rsv_time[1], state := "SPR"]
temp[timestamp < flu_time[1] & timestamp >rsv_time[1], state := "PRS"]
temp[timestamp == flu_time[1] & timestamp == rsv_time[1], state := "II"]
temp[timestamp > flu_time[1] & timestamp == rsv_time[1], state := "IPR"]
temp[timestamp == flu_time[1] & timestamp > rsv_time[1], state := "PRI"]
temp[timestamp > flu_time[1] & timestamp > rsv_time[1], state := "RR"]





main_17[,timestep := c(1:82), by=indid]



