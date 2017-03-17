rm(list=ls())

source('~/Git/utilities/rescale.R')  # function to rescale values of x to interval [a,b]

####################################################################
### Simple Guassian example
####################################################################

wts <- function(y,mu,sigma,q=list(mu,sigma),priors=list(mu,lambda)){  
	# calculate importance weights
	prior <- dnorm(mu,priors[[1]],sigma/priors[[2]],log=TRUE)  # prior
	lik <- sapply(mu,function(x) sum(dnorm(y,x,sigma,log=TRUE)))  # likelihood
	q <- dnorm(mu,q[[1]],q[[2]],log=TRUE)  # pdf of importance function
	m <- max(lik)  # for "log-sum-exp" trick to avoid numerical underflow
	w <- exp(prior+lik-q-m)  # importance weights
	w <- w/sum(w)  # normalize importance weights
	w
}

mcmc <- function(y,sigma,lambda,n.mcmc=1000){  # estimate mu via MCMC
	n <- length(y)
	A.inv <- (n/sigma^2+1/(sigma^2/lambda))^-1
	b <- sum(y)/sigma^2
	rnorm(n.mcmc,A.inv*b,sqrt(A.inv))
}


###
### Simulate observations
###

n <- 100  # sample size
mu.true <- 0  # mean
sigma.true <- 1  # standarad deviation
y <- rnorm(n,mu.true,sigma.true)

# Prior specification (a la Cappe et al. 2004)
theta <- 0  # prior mean
lambda <- 0.1  # hyperparameter for prior sd


###
### Step 0: Standard importance sampling scheme
###

n.samp <- 10000  # number of samples from importance function
mu.q <- 0  # mean of importance function
sigma.q <- 10  # sd of importance function
mu <- rnorm(n.samp,mu.q,sigma.q)  # sample from importance function
w <- wts(y,mu,sigma.true,list(mu.q,sigma.q),priors=list(theta,lambda))  # importance weights
plot(mu,w,pch=19,col=rgb(0,0,0,0.35))

# plot(cumsum(w)/c(1:n.samp),type="l")  # Monte Carlo integration

# Initial IS approximation to posterior distribution
mu.is <- sample(mu,10000,w,replace=TRUE)  
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
mean(mu.is)

# Compare to posterior obtained via MCMC
mu.mcmc <- mcmc(y,sigma.true,lambda,n.mcmc=10000)  # posterior of mu via MCMC
lines(density(mu.mcmc),col=4)
hist(mu.mcmc,breaks=100,prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE)
mean(mu.mcmc)


###
### Stesp 1-K: Update importance function, sample mu, calculate importance weights...
###

K <- 10  # number of iterations
for(k in 1:K){ # loop through iterations
	mu.q <- sum(w*mu)  # update mean of importance function
	sigma.q <- sqrt(sum(w*(mu-mu.q)^2))  # update sd of importance function
	mu <- rnorm(n.samp,mu.q,sigma.q)  # sample from importance function
	w <- wts(y,mu,sigma.true,list(mu.q,sigma.q),priors=list(theta,lambda))  # importance weights
}

# Posterior via importance sampling
mu.is <- sample(mu,10000,w,replace=TRUE)
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
mean(mu.is)

# Add posterior via MCMC
hist(mu.mcmc,breaks=100,prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE)
mean(mu.mcmc)

# Add truth
abline(v=mu.true,lty=2,lwd=2)
abline(v=mean(y),lty=2,lwd=2,col=2)


####################################################################
### Guassian example with mixture importance function a la Wraith et al. (2009)
####################################################################

rm(list=ls())

rmix <- function(n,mu,sigma,p){  # simulate random draws from mixture of Guassians
# browser()
	d <- length(mu)  # number of mixture components
	idx <- sample(1:d,n,replace=TRUE,prob=p)  # mixture component idx
	out <- rnorm(n,mu[idx],sigma[idx])
	out
}

dmix <- function(x,mu,sigma,p,sum=TRUE,log=TRUE){  # PDF of Guassian mixture distribution
# browser()
	d <- length(mu)  # number of mixture components
	out <- sapply(1:d,function(y) p[y]*dnorm(x,mu[y],sigma[y]))
	if(sum==TRUE) out <- rowSums(out)
	if(log==TRUE) out <- log(out)
	out
}

wts <- function(y,mu,sigma,  # calculate importance weights
	q=list(mu,sigma,p),priors=list(mu,sigma)){  
			
	prior <- dnorm(mu,priors[[1]],sigma/priors[[2]],log=TRUE)  # prior
	lik <- sapply(mu,function(x) sum(dnorm(y,x,sigma,log=TRUE)))  # likelihood
	q <- dmix(mu,q[[1]],q[[2]],q[[3]],sum=TRUE,log=TRUE)  # pdf of importance function
	m <- max(lik)  # for "log-sum-exp" trick to avoid numerical underflow
	w <- exp(prior+lik-q-m)  # importance weights
	w <- w/sum(w)  # normalize importance weights
	w
}

mcmc <- function(y,sigma,lambda,n.mcmc=1000){   # estimate mu via MCMC
	n <- length(y)
	A.inv <- (n/sigma^2+1/(sigma^2/lambda))^-1
	b <- sum(y)/sigma^2
	rnorm(n.mcmc,A.inv*b,sqrt(A.inv))
}


###
### Simulate observations
###

n <- 1000  # sample size
mu.true <- 0  # mean
sigma.true <- 1  # standarad deviation
y <- rnorm(n,mu.true,sigma.true)  # simulated observations
hist(y,breaks=20)


###
### Prior specification (a la Cappe et al. 2004)
###

theta <- 0  # prior mean
lambda <- 0.1  # hyperparameter for prior sd


###
### Specify parameters in mixture importance function
###

d <- 3  # number of mixture components
mu.q <- seq(-1,1)  # mean of mixture components
sigma.q <- rep(0.25,d)  # standard deviation of mixture components
p.q <- rep(1/d,d)  # mixture probabilities


###
### Step 0: Standard importance sampling scheme
###

n.samp <- 10000  # number of samples from importance function
mu <- rmix(n.samp,mu.q,sigma.q,p.q)  # sample from importance function
# hist(mu,breaks=100);abline(v=mu.q,col=2)
w <- wts(y,mu,sigma.true,q=list(mu.q,sigma.q,p.q),priors=list(theta,lambda))  # importance weights
plot(mu,w,pch=19,col=rgb(0,0,0,0.35))

# plot(cumsum(w)/c(1:n.samp),type="l")

# Initial IS approximation to posterior distribution
mu.is <- sample(mu,10000,w,replace=TRUE)
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
mean(mu.is)

# Posterior of mu obtained via MCMC
mu.mcmc <- mcmc(y,sigma.true,lambda,n.mcmc=10000)
lines(density(mu.mcmc),col=4)
hist(mu.mcmc,breaks=100,prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE)
mean(mu.mcmc)


###
### Steps 1-K: Update importance function, sample mu, calculate importance weights
###

K <- 10  # number of iterations
for(k in 1:K){ # loop through iterations
	rho <- dmix(mu,mu.q,sigma.q,p.q,sum=FALSE,log=FALSE)   
	# rho <- rho/matrix(rowSums(rho),nrow(rho),ncol(rho))  
	rho <- exp(log(rho)-log(rowSums(rho)))  # avoid numerical underflow
	w <- wts(y,mu,sigma.true,  # importance weights
		q=list(mu.q,sigma.q,p.q),priors=list(theta,lambda))  
	# p.q <- colSums(w*rho)  # update mixture probabilities
	p.q <- colSums(exp(log(w)+log(rho)))  # avoid numerical underflow
	mu.q <- colSums((w*mu)*rho)/p.q  # update means of mixture components
	sigma.q <- (mu-matrix(mu.q,length(mu),length(mu.q),byrow=TRUE))^2
	sigma.q <- sqrt(colSums(w*sigma.q*rho)/p.q)  # update sd of mixture components
	mu <- rmix(n.samp,mu.q,sigma.q,p.q)  # sample from importance function
}

# Posterior via importance sampling
mu.is <- sample(mu,10000,w,replace=TRUE)
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
mean(mu.is)

# Add posterior via MCMC
hist(mu.mcmc,breaks=100,prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE)
mean(mu.mcmc)

# Add truth
abline(v=mu.true,lty=2,lwd=2)
abline(v=mean(y),lty=2,lwd=2,col=2)


####################################################################
### Mixture of Guassians example of Cappe et al. (2004)
####################################################################

rm(list=ls())

library(akima)

dmix <- function(x,mu,sigma,p,sum=TRUE,log=TRUE){  # PDF of Guassian mixture distribution
# browser()
	d <- length(mu)  # number of mixture components
	out <- sapply(1:d,function(y) p[y]*dnorm(x,mu[y],sigma))
	if(sum==TRUE) out <- rowSums(out)
	if(log==TRUE) out <- log(out)
	out
}

get.wts <- function(y,mu,sigma,p,  # calculate importance weights
	q=list(mu,sigma),priors=list(mu,lambda)){  
	prior <- rowSums(dnorm(mu,priors[[1]],sigma/priors[[2]],log=TRUE))  # prior
	lik <-  apply(mu,1,function(x) sum(dmix(y,x,sigma,p,sum=TRUE,log=TRUE)))  # likelihood	
	d.q <- dnorm(mu[,1],q[[1]][1],q[[2]][1],log=TRUE)+ # pdf of importance function
		dnorm(mu[,2],q[[1]][2],q[[2]][2],log=TRUE)
	m <- max(lik)  # for "log-sum-exp" trick to avoid numerical underflow
	w <- exp(prior+lik-d.q-m)  # importance weights
	w <- w/sum(w)  # normalize importance weights
	w
}


###
### Simulate observations
###

n <- 1000 # sample size
mu.true <- c(0,2)  # mean
sigma.true <- 0.5  # standarad deviation
p.true <- c(0.2,1-0.2)  # mixture proportions
J <- length(mu.true)  # number of mixture components in model

idx <- sample(1:J,n,replace=TRUE,prob=p.true)  # mixture component idx
y <- rnorm(n,mu.true[idx],sigma.true)  # simulated observations
hist(y,breaks=100)


###
### Prior specification (a la Cappe et al. 2004)
###

theta <- 0  # prior mean
tau <- 10  # hyperparameter for prior sd

a <- 0
b <- 10

###
### Step 0: Standard importance sampling scheme
###

n.samp <- 10000  # sample size from importance function

q.mu <- list(mean=rep(0,J),sd=rep(10,J))  # parameters of importance fxn for mu
q.sigma <- list(mean=1,sd=0.5)  # parameters of importance fxn for sigma


###
### Stesp 1-K: Update importance function, sample mu, calculate importance weights
###

K <- 10  # number of iterations
for(k in 1:K){ # loop through iterations
	
	###
	### Update mu
	### 

	# Sample mu
	mu <- sapply(1:J,function(x) rnorm(n.samp,q.mu$mean[x],q.mu$sd[x]))

	# Calculate importance weights for mu
	loglik <- apply(mu,1,function(x) sum(dmix(y,x,q.sigma$mean,p.true,sum=TRUE,log=TRUE)))
	prior <- rowSums(dnorm(mu,theta,tau,log=TRUE))
	q.density <- rowSums(sapply(1:J,function(x) dnorm(mu[,x],q.mu$mean[x],q.mu$sd[x],log=TRUE)))
	mu.wts <- exp(prior+loglik-q.density-max(loglik))
	mu.wts <- mu.wts/sum(mu.wts)

	# contour(interp(x=mu[,1],y=mu[,2],z=wts))
	# abline(v=mu.true[1],h=mu.true[2],col=2)

	# Update parameters of importance function for mu
	q.mu <- list(mean=colSums(mu.wts*mu),
		sd=sqrt(colSums(mu.wts*t(apply(mu,1,function(x) x-q.mu$mean))^2)))


	###
	### Update sigma
	###

	# Sample sigma
	sigma <- exp(rnorm(n.samp,log(q.sigma$mean),q.sigma$sd))
	# mean(sigma)
	# hist(sigma,breaks=100)
	
	# Calculate importance weights for sigma
	loglik <- sapply(sigma,function(x) sum(dmix(y,q.mu$mean,x,p.true,sum=TRUE,log=TRUE)))
	prior <- dunif(sigma,a,b,log=TRUE)
	q.density <- dnorm(log(sigma),log(q.sigma$mean),q.sigma$sd,log=TRUE)
	sigma.wts <- exp(prior+loglik-q.density-max(loglik))
	sigma.wts <- sigma.wts/sum(sigma.wts)

	# Update parameters of importance function for sigma
	q.sigma <- list(mean=sum(sigma.wts*sigma), 
		sd=sqrt(sum(sigma.wts*(sigma-q.sigma$sd)^2)))
}



# Posterior via importance sampling
mu.is <- sample(1:n.samp,10000,mu.wts,replace=TRUE)
mu.is <- c(mu[mu.is,])
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25));abline(v=mu.true,lty=2)

sigma.is <- sample(1:n.samp,10000,sigma.wts,replace=TRUE)
sigma.is <- sigma[sigma.is]
hist(sigma.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25));abline(v=sigma.true,lty=2)

# Class memberships
lik <- dmix(y,q.mu$mean,q.sigma$mean,p.true,sum=FALSE,log=FALSE)
p <- lik/rowSums(lik)
idx <- apply(p,1,function(x) sample(1:J,1,prob=x))
boxplot(y~idx);abline(h=mu.true,col=2)
plot(y,p[,1]);points(y,p[,2],col=2)
