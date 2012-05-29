measurements.pr.day = 288
num.training.samples = 8000

dat = read.csv('data/myklimdata-noheader.csv', header=FALSE)

data.set.size = length(dat[,1])
num.test.samples = data.set.size - num.training.samples

train.idx = 1:num.training.samples
test.idx = (num.training.samples+1):data.set.size

dat.train = dat[train.idx,]
dat.test = dat[test.idx,]

times = dat.train[,1]
day.of.year = dat.train[,2]
dates = as.Date("2002-01-01") + day.of.year
power.production = dat.train[,3]

times.test = dat.test[,1]
day.of.year.test = dat.test[,2]
dates.test = as.Date("2002-01-01") + day.of.year.test
power.production.test = dat.test[,3]

day.shifts = which(diff(ceiling(day.of.year))==1)
first.day.shift = day.shifts[1]
first.day.shift.test = which(diff(ceiling(day.of.year.test))==1)[1]
