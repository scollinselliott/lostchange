library(coda)
library(ggplot2)

CNcoins = read.csv("finalresults-['CN'].csv")
PTcoins = read.csv("finalresults-['PT'].csv")
PIcoins = read.csv("finalresults-['PI'].csv")
MZcoins = read.csv("finalresults-['MZ'].csv")

CNrmean = read.csv("rmeancollected-['CN'].csv")
PTrmean = read.csv("rmeancollected-['PT'].csv")
PIrmean = read.csv("rmeancollected-['PI'].csv")
MZrmean = read.csv("rmeancollected-['MZ'].csv")


MZcoinsdata = MZcoins
colnames(MZcoinsdata) <- c("year","mu", "musub", "nu", "xi")
as.data.frame.matrix(MZcoinsdata) 
x <- seq(0, 2.5, length = 1000)
lower <- c()
upper <- c()
for (i in 1:dim(MZcoins)[1]){
	pdf1 <- rgamma(x, shape  = MZcoinsdata[i,]$nu, scale = MZcoinsdata[i,]$xi)
	z <- HPDinterval(as.mcmc(pdf1), prob=0.5)
	lower <- c(lower, z[1])
	upper <- c(upper, z[2]) 
}
MZcoinsdata["upper"] <- upper
MZcoinsdata["lower"] <- lower

CNcoinsdata = CNcoins
colnames(CNcoinsdata) <- c("year","mu", "musub", "nu", "xi")
as.data.frame.matrix(CNcoinsdata) 
x <- seq(0, 2.5, length = 1000)
lower <- c()
upper <- c()
for (i in 1:dim(CNcoins)[1]){
	pdf1 <- rgamma(x, shape  = CNcoinsdata[i,]$nu, scale = CNcoinsdata[i,]$xi)
	z <- HPDinterval(as.mcmc(pdf1), prob=0.5)
	lower <- c(lower, z[1])
	upper <- c(upper, z[2]) 
}
CNcoinsdata["upper"] <- upper
CNcoinsdata["lower"] <- lower

PIcoinsdata = PIcoins
colnames(PIcoinsdata) <- c("year","mu", "musub", "nu", "xi")
as.data.frame.matrix(PIcoinsdata) 
x <- seq(0, 2.5, length = 1000)
lower <- c()
upper <- c()
for (i in 1:dim(CNcoins)[1]){
	pdf1 <- rgamma(x, shape  = PIcoinsdata[i,]$nu, scale = PIcoinsdata[i,]$xi)
	z <- HPDinterval(as.mcmc(pdf1), prob=0.5)
	lower <- c(lower, z[1])
	upper <- c(upper, z[2]) 
}
PIcoinsdata["upper"] <- upper
PIcoinsdata["lower"] <- lower


PTcoinsdata = PTcoins
colnames(PTcoinsdata) <- c("year","mu", "musub", "nu", "xi")
as.data.frame.matrix(PTcoinsdata) 
x <- seq(0, 2.5, length = 1000)
lower <- c()
upper <- c()
for (i in 1:dim(PTcoins)[1]){
	pdf1 <- rgamma(x, shape  = PTcoinsdata[i,]$nu, scale = PTcoinsdata[i,]$xi)
	z <- HPDinterval(as.mcmc(pdf1), prob=0.5)
	lower <- c(lower, z[1])
	upper <- c(upper, z[2]) 
}
PTcoinsdata["upper"] <- upper
PTcoinsdata["lower"] <- lower


ggplot(data=MZcoinsdata, aes(year)) +  
geom_ribbon(aes(ymin=MZcoinsdata$lower,ymax=MZcoinsdata$upper), fill="#000000", alpha="0.10") +
geom_ribbon(aes(ymin=CNcoinsdata$lower,ymax=CNcoinsdata$upper), fill="#CC0000", alpha="0.10") +
geom_ribbon(aes(ymin=PIcoinsdata$lower,ymax=PIcoinsdata$upper), fill="#0000FF", alpha="0.10") +
geom_ribbon(aes(ymin=PTcoinsdata$lower,ymax=PTcoinsdata$upper), fill="#009933", alpha="0.10") +
  theme_classic()



z <- 0
for (i in 1:dim(MZcoinsdata)[1]){
zz <- MZcoinsdata$mu[i]
if (zz > 0) { z <- z + 1 }
}
MZbias <- (1 / z) * sum(MZcoinsdata$mu - MZcoinsdata$musub)

z <- 0
for (i in 1:dim(CNcoinsdata)[1]){
zz <- CNcoinsdata$mu[i]
if (zz > 0) { z <- z + 1 }
}
CNbias <- (1 / z) * sum(CNcoinsdata$mu - CNcoinsdata$musub)

z <- 0
for (i in 1:dim(PTcoinsdata)[1]){
zz <- PTcoinsdata$mu[i]
if (zz > 0) { z <- z + 1 }
}
PTbias <- (1 / z) * sum(PTcoinsdata$mu - PTcoinsdata$musub)

z <- 0
for (i in 1:dim(PIcoinsdata)[1]){
zz <- PIcoinsdata$mu[i]
if (zz > 0) { z <- z + 1 }
}
PIbias <- (1 / z) * sum(PIcoinsdata$mu - PIcoinsdata$musub)


