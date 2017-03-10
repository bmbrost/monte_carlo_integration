###################################################################
### Population Monte Carlo implementation of hidden semi-Markov model of Cappe et al. 2004
###################################################################

###
### Simulate Negative binomially-distributed soujourn times of latten process
###

s <- c(3,5)  # size parameter of negative binomial distribution
p <- 0.5  # probability of success for negative binomial distribution
hist(rnbinom(1000,size=s[1],prob=p))
hist(rnbinom(1000,size=s[2],prob=p))

J <- 100  # number of sojourns
idx <- seq(1,J,by=2)
d <- numeric(J)  # latent sojourn times
d[idx] <- rnbinom(J/2,s[1],p)  # sojourn times in first latent state
d[-idx] <- rnbinom(J/2,s[2],p)  # sojourn times in second latent state
hist(d[idx],prob=TRUE,ylim=c(0,1),xlim=c(0,max(d)+1),col=rgb(1,0,0,0.25),breaks=5)
hist(d[-idx],prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE,breaks=10)

# T <- sum(d) + sum(d==0)  # total number of time steps (including sojourns of duration 0)
idx <- rep(0:1,J/2)  # latent state indicator variable
z <- unlist(sapply(1:J,function(x) rep(idx[x],d[x])))  # latent states for times t in 1:T
d <- d[d!=0]  # sojourn durations (excluding sojourns of duration 0)
T <- sum(d)  # total number of time steps (excluding sojourns of duration 0)
length(z)
T


###
### Simulate Gamma-distributed soujourn times of latten process
###

s <- c(3,3)  # shape parameter of gamma distribution
lambda <- c(0.05,0.1)  #  rate parameter of gamma distribution
mean(rgamma(1000,shape=s[2],rate=lambda[2]))  # equals s/lambda
hist(rgamma(1000,shape=s[1],rate=lambda[1]))
hist(rgamma(1000,shape=s[2],rate=lambda[2]))

J <- 100  # number of sojourns
idx <- seq(1,J,by=2)
d <- numeric(J)  # latent sojourn times
d[idx] <- rgamma(J/2,shape=s[1],rate=lambda[1])  # sojourn times in first latent state
d[-idx] <- rgamma(J/2,shape=s[2],rate=lambda[2])  # sojourn times in second latent state
hist(d[idx],prob=TRUE,ylim=c(0,0.05),xlim=c(0,max(d)+1),col=rgb(1,0,0,0.25),breaks=5)
hist(d[-idx],prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE,breaks=10)

sojourn <- data.frame(start=0,end=0,duration=d,state=rep(1:2,J/2))
sojourn$start[-1] <- cumsum(sojourn$d[-J])
sojourn$end[-J] <- sojourn$start[-1]
sojourn$end[J] <- sojourn$start[J]+sojourn$duration[J]
sojourn

###
### Simulate true latent state
###

T <- max(sojourn$end)
t <- 1:T
z <- sapply(t,function(x) findInterval(x,sojourn$start))
z <- sojourn$state[z]


###
### Simulate observation process
###

mu.true <- c(0,2)
sigma.true <- 0.5

y <- rnorm(T,mu.true[z],sigma.true)
boxplot(y~z)
plot(1:T,y,type="l")
abline(h=mu.true,col=2,lty=2)


###
### Prior distribution parameters
###

theta <- 0
tau <- 10


