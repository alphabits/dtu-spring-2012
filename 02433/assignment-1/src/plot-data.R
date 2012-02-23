# All paths are relative to assignment root dir
source('src/utils.R')

data = read.table("data/soap.txt")
x = data[,1]
mean.x = mean(x)
sd.x = sd(x)

plot(x, type="l", xlab="Week", ylab="Weekly soap sale",
     main="Sales of soap product for 242 consecutive weeks")
