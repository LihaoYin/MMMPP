##################Chase Credit Card Dataset
path = paste(getwd(), "/transactions.csv",sep="")
file = read.csv(path)

Months = c(0,31,60,91,121,152,182,213,244,274,305,335)

timeStamps = c(); days = c()
for(i in 1:nrow(file)){
  t = as.character(file$transactionDateTime[i])
  day = strsplit(t,"T")[[1]][1]
  day = strsplit(day, '-')[[1]]
  time = strsplit(t,"T")[[1]][2]
  time = strsplit(time,":")[[1]]
  days = c(days, Months[as.integer(day[2])] + as.integer(day[3]))
  timeStamps = c(timeStamps, as.integer(time[1])/24+as.integer(time[2])/1440+as.integer(time[3])/86400)
}

file = data.frame(accountNumber=file$accountNumber, timeStamps=timeStamps, days=days, category=rep(1,length(days))
write.csv(file, file = "dataset.csv", row.names = FALSE)
