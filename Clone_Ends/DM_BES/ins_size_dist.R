
my.dat <-
  read.table("puke")

head(my.dat)
nrow(my.dat)

my.dat$V3 <- my.dat$V2-my.dat$V1



hist(log10(my.dat$V3))
rug (log10(my.dat$V3))

my.lower <- 2500

abline(v=log10(my.lower), col=2)

hist(log10(my.dat[my.dat$V3>my.lower,"V3"]))
rug (log10(my.dat[my.dat$V3>my.lower,"V3"]))

my.upper <- 400000

abline(v=log10(my.upper), col=2)



hist(log10(my.dat[my.dat$V3>my.lower & my.dat$V3<my.upper,"V3"]))
rug (log10(my.dat[my.dat$V3>my.lower & my.dat$V3<my.upper,"V3"]))

hist(my.dat[my.dat$V3>my.lower & my.dat$V3<my.upper,"V3"])

