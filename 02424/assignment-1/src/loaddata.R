dat = read.table("data/dioxin.csv", sep=",", head=TRUE)
dat$OXYGEN = relevel(dat$OXYGEN, ref="N")
dat$LOAD = relevel(dat$LOAD, ref="N")
dat$PRSEK = relevel(dat$PRSEK, ref="N")
dat$TIME = as.factor(dat$TIME)

dat2 = dat[dat$OXYGEN!="N",]
