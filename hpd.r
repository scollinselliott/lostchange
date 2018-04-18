library(coda)
library(ggplot2)
library(fitdistrplus)
library(MCMCpack)
library(fitdistrplus)

#setwd("C:/...") 
#import bootstrap and subsampling results from python

CN.boot = read.csv("rmeancollected-boot-['CN'].csv")
PT.boot = read.csv("rmeancollected-boot-['PT'].csv")
PI.boot = read.csv("rmeancollected-boot-['PI'].csv")
MZ.boot = read.csv("rmeancollected-boot-['MZ'].csv")

CN.sub = read.csv("rmeancollected-subsampled-['CN'].csv")
PT.sub = read.csv("rmeancollected-subsampled-['PT'].csv")
PI.sub = read.csv("rmeancollected-subsampled-['PI'].csv")
MZ.sub = read.csv("rmeancollected-subsampled-['MZ'].csv")

daterange <- MZcoins[,1]

#create bootstrap and subsampling dataframes
MZ.b = data.frame(daterange)
MZ.s = data.frame(daterange)
PI.b = data.frame(daterange)
PI.s = data.frame(daterange)
CN.b = data.frame(daterange)
CN.s = data.frame(daterange)
PT.b = data.frame(daterange)
PT.s = data.frame(daterange)

lower.b <- c()
upper.b <- c()
lower.s <- c()
upper.s <- c()

#calculate hpd for both bootstrap and subsampling results
for (i in 1:dim(MZ.sub)[1]){
	samples.b <- as.numeric(as.vector( MZ.boot[i,] ))
	samples.s <- as.numeric(as.vector( MZ.sub[i,] ))

	if(sum(samples.b) > 0){
	nonzerosamples.b <- samples.b[samples.b > 0]
	x.b <- nonzerosamples.b
	den.b <- density(x.b)
	dat.b <- data.frame(x = den.b$x, y = den.b$y)
	fit.gamma.b <- fitdist(x.b, "gamma", lower = c(0, 0))
	shape.b <- fit.gamma.b$estimate[1]
	rate.b <- fit.gamma.b$estimate[2]
	res.b <- rgamma(100000, shape = shape.b, rate = rate.b)
	hpd.b <- HPDinterval(as.mcmc(res.b), prob = 0.5)

	lower.b <- c(lower.b, hpd.b[1])
	upper.b <- c(upper.b, hpd.b[2]) }

	if(sum(samples.b) == 0){
	lower.b <- c(lower.b, 0)
	upper.b <- c(upper.b, 0) 
	}

	if (sum(samples.s) > 0){
	nonzerosamples.s <- samples.s[samples.s > 0]
	x.s <- nonzerosamples.s
	den.s <- density(x.s)
	dat.s <- data.frame(x = den.s$x, y = den.s$y)
	fit.gamma.s <- fitdist(x.s, "gamma", lower = c(0, 0))
	shape.s <- fit.gamma.s$estimate[1]
	rate.s <- fit.gamma.s$estimate[2]
	res.s <- rgamma(100000, shape = shape.s, rate = rate.s)
	hpd.s <- HPDinterval(as.mcmc(res.s), prob = 0.5)

	lower.s <- c(lower.s, hpd.s[1])
	upper.s <- c(upper.s, hpd.s[2]) }


	if(sum(samples.s) == 0){
	lower.s <- c(lower.s, 0)
	upper.s <- c(upper.s, 0) 
	}
}
MZ.b["upper"] <- upper.b
MZ.b["lower"] <- lower.b
MZ.s["upper"] <- upper.s
MZ.s["lower"] <- lower.s
lower.b <- c()
upper.b <- c()
lower.s <- c()
upper.s <- c()

for (i in 1:dim(PT.sub)[1]){
	samples.b <- as.numeric(as.vector( PT.boot[i,] ))
	samples.s <- as.numeric(as.vector( PT.sub[i,] ))

	if(sum(samples.b) > 0){
	nonzerosamples.b <- samples.b[samples.b > 0]
	x.b <- nonzerosamples.b
	den.b <- density(x.b)
	dat.b <- data.frame(x = den.b$x, y = den.b$y)
	fit.gamma.b <- fitdist(x.b, "gamma", lower = c(0, 0))
	shape.b <- fit.gamma.b$estimate[1]
	rate.b <- fit.gamma.b$estimate[2]
	res.b <- rgamma(100000, shape = shape.b, rate = rate.b)
	hpd.b <- HPDinterval(as.mcmc(res.b), prob = 0.5)

	lower.b <- c(lower.b, hpd.b[1])
	upper.b <- c(upper.b, hpd.b[2]) }

	if(sum(samples.b) == 0){
	lower.b <- c(lower.b, 0)
	upper.b <- c(upper.b, 0) 
	}

	if (sum(samples.s) > 0){
	nonzerosamples.s <- samples.s[samples.s > 0]
	x.s <- nonzerosamples.s
	den.s <- density(x.s)
	dat.s <- data.frame(x = den.s$x, y = den.s$y)
	fit.gamma.s <- fitdist(x.s, "gamma", lower = c(0, 0))
	shape.s <- fit.gamma.s$estimate[1]
	rate.s <- fit.gamma.s$estimate[2]
	res.s <- rgamma(100000, shape = shape.s, rate = rate.s)
	hpd.s <- HPDinterval(as.mcmc(res.s), prob = 0.5)

	lower.s <- c(lower.s, hpd.s[1])
	upper.s <- c(upper.s, hpd.s[2]) }

	if(sum(samples.s) == 0){
	lower.s <- c(lower.s, 0)
	upper.s <- c(upper.s, 0) 
	}
}
PT.b["upper"] <- upper.b
PT.b["lower"] <- lower.b
PT.s["upper"] <- upper.s
PT.s["lower"] <- lower.s
lower.b <- c()
upper.b <- c()
lower.s <- c()
upper.s <- c()

for (i in 1:dim(PI.sub)[1]){
	samples.b <- as.numeric(as.vector( PI.boot[i,] ))
	samples.s <- as.numeric(as.vector( PI.sub[i,] ))

	if(sum(samples.b) > 0){
	nonzerosamples.b <- samples.b[samples.b > 0]
	x.b <- nonzerosamples.b
	den.b <- density(x.b)
	dat.b <- data.frame(x = den.b$x, y = den.b$y)
	fit.gamma.b <- fitdist(x.b, "gamma", lower = c(0, 0))
	shape.b <- fit.gamma.b$estimate[1]
	rate.b <- fit.gamma.b$estimate[2]
	res.b <- rgamma(100000, shape = shape.b, rate = rate.b)
	hpd.b <- HPDinterval(as.mcmc(res.b), prob = 0.5)

	lower.b <- c(lower.b, hpd.b[1])
	upper.b <- c(upper.b, hpd.b[2]) }

	if(sum(samples.b) == 0){
	lower.b <- c(lower.b, 0)
	upper.b <- c(upper.b, 0) 
	}

	if (sum(samples.s) > 0){
	nonzerosamples.s <- samples.s[samples.s > 0]
	x.s <- nonzerosamples.s
	den.s <- density(x.s)
	dat.s <- data.frame(x = den.s$x, y = den.s$y)
	fit.gamma.s <- fitdist(x.s, "gamma", lower = c(0, 0))
	shape.s <- fit.gamma.s$estimate[1]
	rate.s <- fit.gamma.s$estimate[2]
	res.s <- rgamma(100000, shape = shape.s, rate = rate.s)
	hpd.s <- HPDinterval(as.mcmc(res.s), prob = 0.5)

	lower.s <- c(lower.s, hpd.s[1])
	upper.s <- c(upper.s, hpd.s[2]) }


	if(sum(samples.s) == 0){
	lower.s <- c(lower.s, 0)
	upper.s <- c(upper.s, 0) 
	}
}
PI.b["upper"] <- upper.b
PI.b["lower"] <- lower.b
PI.s["upper"] <- upper.s
PI.s["lower"] <- lower.s
lower.b <- c()
upper.b <- c()
lower.s <- c()
upper.s <- c()

for (i in 1:dim(CN.sub)[1]){
	samples.b <- as.numeric(as.vector( CN.boot[i,] ))
	samples.s <- as.numeric(as.vector( CN.sub[i,] ))

	if(sum(samples.b) > 0){
	nonzerosamples.b <- samples.b[samples.b > 0]
	x.b <- nonzerosamples.b
	den.b <- density(x.b)
	dat.b <- data.frame(x = den.b$x, y = den.b$y)
	fit.gamma.b <- fitdist(x.b, "gamma", lower = c(0, 0))
	shape.b <- fit.gamma.b$estimate[1]
	rate.b <- fit.gamma.b$estimate[2]
	res.b <- rgamma(100000, shape = shape.b, rate = rate.b)
	hpd.b <- HPDinterval(as.mcmc(res.b), prob = 0.5)

	lower.b <- c(lower.b, hpd.b[1])
	upper.b <- c(upper.b, hpd.b[2]) }

	if(sum(samples.b) == 0){
	lower.b <- c(lower.b, 0)
	upper.b <- c(upper.b, 0) 
	}

	if (sum(samples.s) > 0){
	nonzerosamples.s <- samples.s[samples.s > 0]
	x.s <- nonzerosamples.s
	den.s <- density(x.s)
	dat.s <- data.frame(x = den.s$x, y = den.s$y)
	fit.gamma.s <- fitdist(x.s, "gamma", lower = c(0, 0))
	shape.s <- fit.gamma.s$estimate[1]
	rate.s <- fit.gamma.s$estimate[2]
	res.s <- rgamma(100000, shape = shape.s, rate = rate.s)
	hpd.s <- HPDinterval(as.mcmc(res.s), prob = 0.5)

	lower.s <- c(lower.s, hpd.s[1])
	upper.s <- c(upper.s, hpd.s[2]) }


	if(sum(samples.s) == 0){
	lower.s <- c(lower.s, 0)
	upper.s <- c(upper.s, 0) 
	}
}
CN.b["upper"] <- upper.b
CN.b["lower"] <- lower.b
CN.s["upper"] <- upper.s
CN.s["lower"] <- lower.s

