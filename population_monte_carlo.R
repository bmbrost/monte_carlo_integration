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

dmix <- function(x,mu,sigma,p,sum=TRUE,log=TRUE){  # PDF of Guassian mixture distribution
# browser()
	d <- length(mu)  # number of mixture components
	out <- sapply(1:d,function(y) p[y]*dnorm(x,mu[y],sigma[y]))
	if(sum==TRUE) out <- rowSums(out)
	if(log==TRUE) out <- log(out)
	out
}

wts <- function(y,mu,sigma,p,  # calculate importance weights
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

n <- 20 # sample size
mu.true <- c(-1,1)  # mean
sigma.true <- c(0.5,0.5)  # standarad deviation
p.true <- c(0.33,1-0.33)  # mixture proportions

idx <- sample(1:length(mu.true),n,replace=TRUE,prob=p.true)  # mixture component idx
y <- rnorm(n,mu.true[idx],sigma.true[idx])  # simulated observations
hist(y,breaks=20)


###
### Prior specification (a la Cappe et al. 2004)
###

theta <- 0  # prior mean
lambda <- 0.1  # hyperparameter for prior sd


###
### Step 0: Standard importance sampling scheme
###

n.samp <- 10000  # sample size from importance function
mu1.q <- mu2.q <- 0  # means of importance functions
sigma1.q <- sigma2.q <- 10  # sd of importance functions
mu1 <- rnorm(n.samp,mu1.q,sigma1.q)  # sample mu1 from importance function
mu2 <- rnorm(n.samp,mu2.q,sigma2.q)  # sample mu2 from importance function
w <- wts(y,cbind(mu1,mu2),sigma.true,p.true,  # importance weights
	list(c(mu1.q,mu2.q),c(sigma1.q,sigma2.q)),priors=list(theta,lambda))

# plot(cumsum(w)/c(1:n.samp),type="l")

# Initial IS approximation to posterior distribution
mu.is <- sample(1:n.samp,10000,w,replace=TRUE)
mu.is <- c(mu1[mu.is],mu2[mu.is])
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=mu.true,lty=2)

###
### Stesp 1-K: Update importance function, sample mu, calculate importance weights
###

K <- 10  # number of iterations
for(k in 1:K){ # loop through iterations
	mu1.q <- sum(w*mu1)  # update mean of importance function
	mu2.q <- sum(w*mu2)  # update mean of importance function
	sigma1.q <- sqrt(sum(w*(mu1-mu1.q)^2))  # update sd of importance function
	sigma2.q <- sqrt(sum(w*(mu2-mu2.q)^2))  # update sd of importance function
	mu1 <- rnorm(n.samp,mu1.q,sigma1.q)  # sample mu1 from importance function
	mu2 <- rnorm(n.samp,mu2.q,sigma2.q)  # sample mu2 from importance function
	w <- wts(y,cbind(mu1,mu2),sigma.true,p.true,  # importance weights
		list(c(mu1.q,mu2.q),c(sigma1.q,sigma2.q)),priors=list(theta,lambda))
}

# Posterior via importance sampling
mu.is <- sample(1:n.samp,10000,w,replace=TRUE)
mu.is <- c(mu1[mu.is],mu2[mu.is])
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=mu.true,lty=2)






