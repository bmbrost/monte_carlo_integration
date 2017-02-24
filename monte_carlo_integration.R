rm(list=ls())

source('~/Git/utilities/rescale.R')  # function to rescale values to interval [a,b]

k <- function(x,a=0.4,b=0.08){  # arbitrary target distribution
	exp(a*(x-a)^2-b*x^4)
}

curve(k(x),from=-4,to=4)  # plot of unnormalized target distribution



###
### Classical Monte Carlo integration
###

n <- 100000  # Monte Carlo sample size
a <- -5  # lower bound of uniform random variable
b <- 5  # upper bound of univorm random variable
x <- runif(n,a,b)  # uniform random variable
w <- (b-a)*k(x)  # importance weights
# w <- k(x)/dunif(x,a,b)  # importance weights; same as above
plot(cumsum(w)/c(1:n),type="l")  # integral relative to sample size
abline(h=mean(w),col=2,lty=2)
mean(w)  # Monte Carlo integral
integrate(k,a,b)  # calculate integral via numerical quadrature

q <- w/sum(w)  # standardized importance weights
sum(x*w/sum(w))  # expected value of target distribution


###
### Importance sampling
###

# Guassian distribution as importance function
mu <- -1.5  # mean of importance function
sd <- 10  # sd of importance function
x <- rnorm(n,mu,sd)  # importance function (or instrumental distribution)

w <- k(x)/dnorm(x,mu,sd)  # importance weights
plot(cumsum(w)/c(1:n),type="l")  # integral relative to sample size
abline(h=mean(w),col=2,lty=2)
mean(w)  # Monte Carlo integral
integrate(k,-10,10)  # calculate integral via numerical quadrature

q <- w/sum(w)  # standardized importance weights
sum(x*w/sum(w))  # expected value of target distribution


# Gaussian mixture as importance function

dmix <- function(x,mu,sd,p){  # Guassian mixture distribution
	p*dnorm(x,mu[1],sd[1])+(1-p)*dnorm(x,mu[2],sd[2])
}

curve(k(x),from=-4,to=4)  # plot of unnormalized target distribution

p <- 0.7  # mixture proportion
mu <- c(-1.9,1.4)  # means of mixture components
sd <- c(1.0,1.0)  # sd of mixture components
idx <- sample(1:2,n,replace=TRUE,prob=c(p,1-p))  # mixture component indicator variable
x <- rnorm(n,mu[idx],sd[idx])  # sample from importance function
hist(x,prob=TRUE)
lines(seq(-10,10,0.1),rescale(k(seq(-10,10,0.1)),0,max(hist(x,prob=TRUE)$density)),col=2)

w <- k(x)/dmix(x,mu,sd,p)  # importance weights
plot(cumsum(w)/c(1:n),type="l")  # integral relative to sample size
abline(h=mean(w),col=2,lty=2)
mean(w)  # Monte Carlo integral
integrate(k,-10,10)  # calculate integral via numerical quadrature

q <- w/sum(w)  # standardized importance weights
sum(x*w/sum(w))  # expected value of target distribution


###
### Sampling importance resampling
###

y <- sample(x,n,q,replace=TRUE)  # simulate realizations from target distribution
hist.y <- hist(y,breaks=100,prob=TRUE)
lines(seq(-10,10,0.1),rescale(k(seq(-10,10,0.1)),0,max(hist.y$density)),col=2)
mean(y)  # expected value
