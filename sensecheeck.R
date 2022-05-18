# checkin things
library(data.table)
library(plyr)
main_db <- data.table(main_db)
main_db$flu <- as.numeric(main_db$flu)
main_db$rsv <- as.numeric(main_db$rsv)
main_db[rsv==1 &flu==1]
main_db[, rounded := round_any(timestep, 7, f = ceiling)]

weekly <-data.table( main_db[,sum(rsv), by=c("rounded", "factor_id")])
weekly_flu<- data.table(main_db[,sum(flu), by=c("rounded", "factor_id")])
weekly[weekly_flu, flu := i.V1, on=c("rounded", "factor_id")]
colnames(weekly)[3] <- "rsv"

total <- 0
nextstep <- 0
nextstep2 <- 0
nextstep3 <- 0
previous <- 0
concurrent <- 0

for(i in 1:length(unique(weekly$factor_id))){
  temp <- weekly[factor_id == i]
  rsv_weeks <- which(temp$rsv >0)
  if(length(rsv_weeks)>1){{rsv_weeks <- rsv_weeks[1]}}
  last_step <- nrow(temp)
  if(!(1 %in% rsv_weeks | last_step %in% rsv_weeks | (last_step-1) %in% rsv_weeks| (last_step-2) %in% rsv_weeks)){
    for(k in rsv_weeks){
      if(temp[k,]$flu > 0){
        concurrent = concurrent+1
      }
      if(temp[k+1,]$flu > 0){
        nextstep = nextstep+1
      }
      if(temp[k+2,]$flu > 0){
        nextstep2 = nextstep2+1
      }
      if(temp[k+3,]$flu > 0){
        nextstep3 = nextstep3+1
      }
      if(temp[k-1,]$flu > 0){
        previous = previous+1
      }      
      if(temp[k-2,]$flu > 0){
        previous2 = previous2+1
      }
      
      total <- total+1
    }
  }
}

precent_flu_same_week <- ((concurrent+nextstep)/total)*100
percent_flu_next_week <- ((nextstep2 + nextstep3)/total)*100
percent_flu_previous_week <- (previous + previous2/total)*100

precent_flu_same_week
percent_flu_next_week
percent_flu_next_week2
percent_flu_next_week3
percent_flu_previous_week


