data = read.table("data/soap.txt")
x = data[,1]
mean.x = mean(x)
sd.x = sd(x)
acf(x)
