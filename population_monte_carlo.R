
source('~/Git/utilities/rescale.R')  # function to rescale values of x to interval [a,b]

####################################################################
### Simple Guassian example
####################################################################

# Model statement:
# y ~ N(mu,sigma^2)
# mu ~ N(theta,tau^2)
# sigma ~ unif(a,b)

rm(list=ls())

source('~/Git/monte_carlo_integration/mcmc/gaussian_mcmc.R')  # estimate mu and sigma via mcmc


###
### Simulate observations
###

n <- 100  # sample size
mu.true <- 2  # mean
sigma.true <- 5  # standarad deviation
y <- rnorm(n,mu.true,sigma.true)
hist(y)

# Prior specifications
theta <- 0  # prior mean for mu
tau <- 10  # prior sd for mu
a <- 0  # uniform prior lower bound for sigma 
b <- 20  # uniofrm prior upper bound for sigma 

# Starting parameters of importance functions
q.mu <- list(mean=0,sd=10)  # importance function for mu
q.sigma <- list(mean=3,sd=3)  # importance function for sigma


###
### Adapt importance function
###

n.samp <- 1000  # number of samples from importance function

T <- 10  # number of iterations
for(t in 1:T){  # iterations to adapt importance function

	if(t==T) n.samp <- 100000

	###
	### Sample mu
	### 
	
	# Update parameters of importance function for mu
	if(t>1){  
		q.mu$sd <- sqrt(sum(mu.wts*(mu-q.mu$mean)^2)) 
		q.mu$mean <- sum(mu.wts*mu)
	}

	# Sample mu from importance function
	mu <- rnorm(n.samp,q.mu$mean,q.mu$sd)
	
	# Calculate log-likelihood, prior, and density of importance samples
	log.lik <- sapply(mu,function(x) sum(dnorm(y,x,q.sigma$mean,log=TRUE)))  # log-likelihood
	prior <- dnorm(mu,theta,tau,log=TRUE)  # prior
	q.density <- dnorm(mu,q.mu$mean,q.mu$sd,log=TRUE)  # pdf of importance function

	# Calculate importance weights
	mu.wts <- exp(log.lik+prior-q.density-  # importance weights
		max(log.lik))  # "log-sum-exp" trick to avoid numerical overflow
	mu.wts <- mu.wts/sum(mu.wts)


	###
	### Update sigma
	###

	# Update parameters of importance function for sigma
	if(t>1){  
		q.sigma$sd <- sqrt(sum(sigma.wts*(sigma-q.sigma$mean)^2))
		q.sigma$mean <- sum(sigma.wts*sigma)
	}

	# Sample sigma from importance function
	sigma <- exp(rnorm(n.samp,log(q.sigma$mean),q.sigma$sd))
	# hist(sigma,breaks=100)
	
	# Calculate log-likelihood, prior, and density of importance samples
	log.lik <- sapply(sigma,function(x) sum(dnorm(y,q.mu$mean,x,log=TRUE)))
	prior <- dunif(sigma,a,b,log=TRUE)
	q.density <- dnorm(log(sigma),log(q.sigma$mean),q.sigma$sd,log=TRUE)
	
	# Calculate importance weights
	sigma.wts <- exp(log.lik+prior-q.density-max(log.lik))
	sigma.wts <- sigma.wts/sum(sigma.wts)
}


###
### Posterior distributions
###

# Fit model via MCMC for comparison
out <- guassian.mcmc(y,theta,tau,a,b,start=list(mu=mu.true,sigma=sigma.true),tune=0.35,10000)

# Posterior distribution of mu via importance sampling
mu.is <- sample(mu,10000,mu.wts,replace=TRUE)  
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=mu.true,lty=2,lwd=2)
abline(v=mean(y),lty=2,lwd=2,col=2)

# Compare to posterior distribution of mu obtained via MCMC
lines(density(out$mu),col=4)
# hist(out$mu,breaks=100,prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE)

# Posterior distribution of sigma via importance sampling
sigma.is <- sample(sigma,10000,sigma.wts,replace=TRUE)  
hist(sigma.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=sigma.true,lty=2,lwd=2)
abline(v=sd(y),lty=2,lwd=2,col=2)

# Compare to posterior distribution of mu obtained via MCMC
lines(density(out$sigma),col=4)
hist(out$sigma,breaks=100,prob=TRUE,col=rgb(0,0,1,0.25),add=TRUE)
# hist(out$sigma,breaks=100,prob=TRUE,col=rgb(0,0,1,0.25),add=FALSE)



####################################################################
### Poisson example with fixed effects
####################################################################

# Model statement:
# y ~ Pois(lambda)
# log(lambda) = X%*%beta
# beta ~ N(theta,tau^2)

rm(list=ls())

###
### Simulate observations
###

n <- 100  # sample size
beta.true <- c(2,1)  # fixed effects
X <- cbind(1,rnorm(n))  # design matrix
qX <- ncol(X)
lambda.true <- exp(X%*%beta.true)  # Poisson rate
y <- rpois(n,lambda.true)
hist(y)

# Prior specifications
theta <- 0  # prior mean for beta
tau <- 10  # prior sd for beta

# Starting parameters of importance functions
q.beta <- list(mean=rep(0,qX),sd=rep(3,qX))  # importance function for mu


###
### Adapt importance function
###

n.samp <- 1000  # number of samples from importance function

T <- 10  # number of iterations
for(t in 1:T){  # iterations to adapt importance function

	if(t==T) n.samp <- 100000

	###
	### Sample mu
	### 
	
	# Update parameters of importance function for mu
	if(t>1){  
		q.beta$sd <- sqrt(colSums(beta.wts*(t(t(beta)-q.beta$mean))^2))
		q.beta$mean <- colSums(beta.wts*beta)
	}

	# Sample mu from importance function
	beta <- sapply(1:qX,function(x) rnorm(n.samp,q.beta$mean[x],q.beta$sd[x]))

	# Calculate likelihood, prior, and density of importance samples
	log.lik <- apply(beta,1,function(x) sum(dpois(y,exp(X%*%x),log=TRUE)))  # log-likelihood
	prior <- rowSums(dnorm(beta,theta,tau,log=TRUE))  # prior	
	q.density <- rowSums(sapply(1:qX,function(x) 
		dnorm(beta[,x],q.beta$mean[x],q.beta$sd[x],log=TRUE)))	

	# Calculate importance weights using 'log-sum-exp' trick to avoid numerical overflow
	beta.wts <- exp(log.lik+prior-q.density-max(log.lik)) 
	beta.wts <- beta.wts/sum(beta.wts)

# plot(beta,cex=sqrt(beta.wts),asp=1);abline(v=beta.true[1],h=beta.true[2],lty=2)

}

###
### Posterior distributions
###

# Posterior of alpha via importance sampling
beta.is <- sample(1:n.samp,10000,beta.wts,replace=TRUE)
beta.is <- beta[beta.is,]
plot(beta.is);abline(v=beta.true[1],h=beta.true[2],col=2)
hist(c(beta.is),breaks=100,prob=TRUE,col=rgb(1,0,0,0.25));abline(v=beta.true,lty=2)

# Compare to posterior distribution of beta obtained via MCMC
source('~/Git/glm/poisson.glm.mcmc.R', chdir = TRUE)
out1 <- poisson.glm.mcmc(y,X,priors=list(sigma.beta=10),start=list(beta=beta.true),
	tune=list(beta=0.02),n.mcmc=10000)
hist(c(out1$beta),breaks=100,prob=TRUE,add=TRUE,col=rgb(0,0,1,0.25))



####################################################################
### Mixture of Gaussians example using integrated likelihood
####################################################################

# Model statement:
# y ~ p_1*N(mu_1,sigma^2)+p_2*N(mu_2,sigma^2)
# mu_1 ~ N(mu.mu,sigma.mu^2)
# mu_2 ~ N(mu.mu,sigma.mu^2)
# logit(p) = beta
# beta ~ N(mu.beta,sigma.beta^2)
# Note that sigma is treated as known to simplify this example

rm(list=ls())

dnorm.mix <- function(y,mu,sigma,p,sum=TRUE,log=TRUE){  # mixture of Guassians density
	K <- length(mu)
	out <- sapply(1:K,function(x) p[x]*dnorm(y,mu[x],sigma,log=FALSE))
	if(sum==TRUE) out <- rowSums(out)
	if(log==TRUE) out <- log(out)
	out
}

logit <- function(x){
	log(x/(1-x))
}

expit <- function(x){
	exp(x)/(1+exp(x))
}

source('~/Git/monte_carlo_integration/mcmc/gaussian_mixture_mcmc.R')  # estimate mu via mcmc

###
### Simulate observations
###

n <- 100  # sample size
mu.true <- c(0,2)  # mean of mixture components
sigma.true <- 0.5  # standard deviation of mixture components
beta.true <- 0.75
p.true <- expit(beta.true)
K <- length(mu.true)  # number of mixture components

z <- sample(1:K,n,replace=TRUE,prob=c(p.true,1-p.true))  # mixture component idx
y <- rnorm(n,mu.true[z],sigma.true)  # simulated observations
hist(y,breaks=100)

# Prior specifications
mu.mu <- 0  # prior mean for mu
sigma.mu <- 10  # prior sd for mu
mu.beta <- 0  # prior mean for beta
sigma.beta <- 2  # prior sd for beta

# Starting parameters of importance functions
q.mu <- list(mean=rep(0,K),sd=rep(3,K))  # importance function for mu
q.beta <- list(mean=0,sd=1)
p <- expit(q.beta$mean)


###
### Adapt importance function
###

n.samp <- 1000  # number of samples from importance function

T <- 10  # number of iterations
for(t in 1:T){  # iterations to adapt importance function

	if(t==T) n.samp <- 100000

	###
	### Sample mu
	### 
	
	# Update parameters of importance function for mu
	if(t>1){  
		q.mu$sd <- sqrt(colSums(mu.wts*(t(t(mu)-q.mu$mean))^2))
		q.mu$mean <- colSums(mu.wts*mu)
	}

	# Sample mu from importance function
	mu <- sapply(1:K,function(x) rnorm(n.samp,q.mu$mean[x],q.mu$sd[x]))
	
	# Calculate likelihood, prior, and density of importance samples
	log.lik <- apply(mu,1,function(x) sum(dnorm.mix(y,x,sigma.true,c(p,1-p))))  # log-likelihood
	prior <- rowSums(dnorm(mu,mu.mu,sigma.mu,log=TRUE))  # prior
	q.density <- sapply(1:K,function(x) dnorm(mu[,x],q.mu$mean[x],q.mu$sd[x],log=TRUE))
	q.density <- rowSums(q.density)  # density of importance samples

	# Calculate importance weights
	mu.wts <- exp(log.lik+prior-q.density-  # importance weights
		max(log.lik))  # "log-sum-exp" trick to avoid numerical overflow
	mu.wts <- mu.wts/sum(mu.wts)

# plot(mu,cex=sqrt(mu.wts),asp=1);abline(v=mu.true[1],h=mu.true[2],lty=2)
# abline(v=mean(y[z==1]),h=mean(y[z==2]),col=2,lty=2)

	###
	### Sample beta
	###
	
	# Update parameters of importance function for mu
	if(t>1){  
		q.beta$sd <- sqrt(sum(beta.wts*(beta-q.beta$mean)^2))
		q.beta$mean <- sum(beta.wts*beta)
		p <- expit(q.beta$mean)
	}

	# Sample beta from importance function; calculate p
	beta <- rnorm(n.samp,q.beta$mean,q.beta$sd)
	p.tmp <- expit(beta)
	
	# Calculate likelihood and prior
	log.lik <- sapply(p.tmp,function(x)  # log-likelihood 
		sum(dnorm.mix(y,q.mu$mean,sigma.true,c(x,1-x))))
	prior <- dnorm(beta,mu.beta,sigma.beta,log=TRUE)
	q.density <- dnorm(beta,q.beta$mean,q.beta$sd,log=TRUE)
	
	# Calculate importance weights
	beta.wts <- exp(log.lik+prior-q.density-  # importance weights
		max(log.lik))  # "log-sum-exp" trick to avoid numerical overflow
	beta.wts <- beta.wts/sum(beta.wts)
# plot(beta,beta.wts);abline(v=beta.true)

}

###
### Posterior distributions
###

# Compare to posterior distribution of mu obtained via MCMC
source('~/Git/monte_carlo_integration/mcmc/gaussian_mixture_mcmc.R')  # estimate mu via mcmc
out <- gaussian.mixture.mcmc(y,theta,tau,start=list(mu=mu.true,sigma=sigma.true,beta=beta.true),
	tune=list(mu1=0.1,mu2=0.1,beta=0.6),n.mcmc=10000,z.truth=z)

# Posterior distribution of mu via importance sampling
mu.is <- sample(1:nrow(mu),10000,mu.wts,replace=TRUE)  
mu.is <- c(mu[mu.is,])
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=mu.true,col=1,lty=2)
abline(v=tapply(y,idx,mean),col=2,lty=2)

hist(c(out$mu),breaks=100,prob=TRUE,add=TRUE,col=rgb(0,0,1,0.25))


# Posterior distribution of mu via importance sampling
beta.is <- sample(1:length(beta),10000,beta.wts,replace=TRUE)  
beta.is <- beta[beta.is]
hist(expit(beta.is),breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=beta.true,col=1,lty=2)

hist(expit(out$beta),breaks=100,prob=TRUE,add=TRUE,col=rgb(0,0,1,0.25))



####################################################################
### Mixture of Gaussians example using indicator variables
####################################################################

# Model statement:
# y ~ N(mu_1,sigma^2), z=0
# y ~ N(mu_2,sigma^2), z=1
# z ~ Bern(p)
# mu_1 ~ N(mu.mu,sigma.mu^2)
# mu_2 ~ N(mu.mu,sigma.mu^2)
# logit(p) = beta
# beta ~ N(mu.beta,sigma.beta^2)

# Note that sigma is treated as known to simplify this example

rm(list=ls())

dnorm.mix <- function(y,mu,sigma,p,sum=TRUE,log=TRUE){  # mixture of Guassians density
	K <- length(mu)
	out <- sapply(1:K,function(x) p[x]*dnorm(y,mu[x],sigma,log=FALSE))
	if(sum==TRUE) out <- rowSums(out)
	if(log==TRUE) out <- log(out)
	out
}

logit <- function(x){
	log(x/(1-x))
}

expit <- function(x){
	exp(x)/(1+exp(x))
}

source('~/Git/monte_carlo_integration/mcmc/gaussian_mixture_mcmc.R')  # estimate mu via mcmc

###
### Simulate observations
###

n <- 100  # sample size
mu.true <- c(0,2)  # mean of mixture components
sigma.true <- 0.5  # standard deviation of mixture components
beta.true <- 0.75
p.true <- expit(beta.true)
K <- length(mu.true)  # number of mixture components

z.true <- sample(1:K,n,replace=TRUE,prob=c(p.true,1-p.true))  # mixture component idx
y <- rnorm(n,mu.true[z.true],sigma.true)  # simulated observations
hist(y,breaks=100)

# Prior specifications
mu.mu <- 0  # prior mean for mu
sigma.mu <- 10  # prior sd for mu
mu.beta <- 0  # prior mean for beta
sigma.beta <- 2  # prior sd for beta

# Starting parameters of importance functions
q.mu <- list(mean=rep(0,K),sd=rep(3,K))  # importance function for mu
q.beta <- list(mean=0,sd=1)
p <- expit(q.beta$mean)
# z <- rbinom(n,1,p)

###
### Adapt importance function
###

n.samp <- 1000  # number of samples from importance function

T <- 10  # number of iterations
for(t in 1:T){  # iterations to adapt importance function

	if(t==T) n.samp <- 100000

	
	###
	### Sample mu
	### 
	
	# Update parameters of importance function for mu
	if(t>1){  
		q.mu$sd <- sqrt(colSums(mu.wts*(t(t(mu)-q.mu$mean))^2))
		q.mu$mean <- colSums(mu.wts*mu)
	}

	# Sample mu from importance function
	mu <- sapply(1:K,function(x) rnorm(n.samp,q.mu$mean[x],q.mu$sd[x]))
	
	# Calculate likelihood, prior, and density of importance samples
	log.lik <- apply(mu,1,function(x) sum(dnorm(y,x[z],sigma.true,log=TRUE))) 
	prior <- rowSums(dnorm(mu,mu.mu,sigma.mu,log=TRUE))  # prior
	q.density <- sapply(1:K,function(x) dnorm(mu[,x],q.mu$mean[x],q.mu$sd[x],log=TRUE))
	q.density <- rowSums(q.density)  # density of importance samples

	# Calculate importance weights
	mu.wts <- exp(log.lik+prior-q.density-  # importance weights
		max(log.lik))  # "log-sum-exp" trick to avoid numerical overflow
	mu.wts <- mu.wts/sum(mu.wts)

# plot(mu,cex=sqrt(mu.wts),asp=1);abline(v=mu.true[1],h=mu.true[2],lty=2)
# abline(v=mean(y[z==1]),h=mean(y[z==2]),col=2,lty=2)

	###
	### Sample beta
	###
	
	# Update parameters of importance function for mu
	if(t>1){  
		q.beta$sd <- sqrt(sum(beta.wts*(beta-q.beta$mean)^2))
		q.beta$mean <- sum(beta.wts*beta)
		p <- expit(q.beta$mean)
	
	}

	# Sample beta from importance function; calculate p
	beta <- rnorm(n.samp,q.beta$mean,q.beta$sd)
	p.tmp <- expit(beta)	
	
	# Calculate likelihood and prior
	log.lik <- sapply(p.tmp,function(x)  # log-likelihood 
		sum(dnorm.mix(y,q.mu$mean,sigma.true,c(x,1-x))))
	prior <- dnorm(beta,mu.beta,sigma.beta,log=TRUE)
	q.density <- dnorm(beta,q.beta$mean,q.beta$sd,log=TRUE)
	
	# Calculate importance weights
	beta.wts <- exp(log.lik+prior-q.density-  # importance weights
		max(log.lik))  # "log-sum-exp" trick to avoid numerical overflow
	beta.wts <- beta.wts/sum(beta.wts)
# plot(beta,beta.wts);abline(v=beta.true)

	###
	### Sample z 
	###
	
	p1 <- p*dnorm(y,q.mu$mean[1],sigma.true)
	p0 <- (1-p)*dnorm(y,q.mu$mean[2],sigma.true)
	p.tilde <- p1/(p1+p0)
	z <- replicate(n.samp,rbinom(n,1,p.tilde))
z <- ifelse(z==0,2,1)
# boxplot(p.tilde~z)

	# Calculate likelihood and prior
	log.lik <- apply(z,2,function(x) sum(dnorm(y,q.mu$mean[x],sigma.true))) # log-likelihood 	
	
	prior <- dnorm(beta,mu.beta,sigma.beta,log=TRUE)
	apply(z,2,function(x) sum(dbinom(z,p)
	q.density <- dnorm(beta,q.beta$mean,q.beta$sd,log=TRUE)
	
	# Calculate importance weights
	beta.wts <- exp(log.lik+prior-q.density-  # importance weights
		max(log.lik))  # "log-sum-exp" trick to avoid numerical overflow
	beta.wts <- beta.wts/sum(beta.wts)
# plot(beta,beta.wts);abline(v=beta.true)


}

###
### Posterior distributions
###

# Compare to posterior distribution of mu obtained via MCMC
source('~/Git/monte_carlo_integration/mcmc/gaussian_mixture_mcmc.R')  # estimate mu via mcmc
out <- gaussian.mixture.mcmc(y,theta,tau,start=list(mu=mu.true,sigma=sigma.true,beta=beta.true),
	tune=list(mu1=0.1,mu2=0.1,beta=0.6),n.mcmc=10000,z.truth=z)

# Posterior distribution of mu via importance sampling
mu.is <- sample(1:nrow(mu),10000,mu.wts,replace=TRUE)  
mu.is <- c(mu[mu.is,])
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=mu.true,col=1,lty=2)
abline(v=tapply(y,idx,mean),col=2,lty=2)

hist(c(out$mu),breaks=100,prob=TRUE,add=TRUE,col=rgb(0,0,1,0.25))


# Posterior distribution of mu via importance sampling
beta.is <- sample(1:length(beta),10000,beta.wts,replace=TRUE)  
beta.is <- beta[beta.is]
hist(expit(beta.is),breaks=100,prob=TRUE,col=rgb(1,0,0,0.25))
abline(v=beta.true,col=1,lty=2)

hist(expit(out$beta),breaks=100,prob=TRUE,add=TRUE,col=rgb(0,0,1,0.25))




















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
### Mixture of Guassians example a la Cappe et al. (2004) with mixure importance function
####################################################################

rm(list=ls())

library(akima)
library(MCMCpack)

dnorm.mix <- function(y,mu,sigma,p,sum=TRUE,log=TRUE){  # PDF of Guassian mixture distribution
# browser()
	J <- length(p)  # number of mixture components
	out <- sapply(1:J,function(x) p[x]*dnorm(y,mu[x],sigma,log=FALSE))
	if(sum==TRUE) out <- rowSums(out)
	if(log==TRUE) out <- log(out)
	out
}


###
### Simulate observations
###

n <- 1000 # sample size
mu.true <- c(0,2)  # mean
sigma.true <- 0.5  # standarad deviation
p.true <- c(0.2,1-0.2)  # mixture proportions
K <- length(mu.true)  # number of mixture components in model

idx <- sample(1:J,n,replace=TRUE,prob=p.true)  # mixture component idx
y <- rnorm(n,mu.true[idx],sigma.true)  # simulated observations
hist(y,breaks=100)


###
### Prior specification (a la Cappe et al. 2004)
###

theta <- 0  # prior mean for mu
tau <- 10  # prior standard deviation for mu

a <- 0  # lower bound for uniform prior on sigma
b <- 10  # upper bound for uniform prior on sigma


###
### Step 0: Standard importance sampling scheme
###

n.samp <- 10000  # sample size from importance function

q.mu <- list(mean=rep(0,K),sd=rep(2,K),p=rep(1/K,K))  # parameters of importance fxn for mu
# q.mu <- list(mean=mu.true,sd=rep(sigma.true,K),p=p.true)  # parameters of importance fxn for mu
q.sigma <- list(mean=1,sd=0.5)  # parameters of importance fxn for sigma


###
### Stesp 1-T: Update importance function, sample mu, calculate importance weights
###

T <- 10  # number of iterations
for(t in 1:T){ # loop through iterations
	
	###
	### Update mu
	### 

	# Sample mu from importance function
	mu <- sapply(1:K,function(x) rnorm(n.samp,q.mu$mean[x],q.mu$sd[x]))

	# Calculate posterior probability of mixture component membership
	rho.tmp <- sapply(1:K,function(x) p.true[x]*dnorm(mu[,x],q.mu$mean[x],q.mu$sd[x],log=FALSE))
	q.density <- log(rowSums(rho.tmp))
	rho <- exp(log(rho.tmp)-q.density)
	
	# Calculate importance weights for mu
	log.lik <- apply(mu,1,function(x) sum(dnorm.mix(y,x,q.sigma$mean,p.true,sum=TRUE,log=TRUE)))
	prior <- rowSums(dnorm(mu,theta,tau,log=TRUE))
	mu.wts <- exp(log.lik+prior-q.density-max(log.lik))
	mu.wts <- mu.wts/sum(mu.wts)

	# contour(interp(x=mu[,1],y=mu[,2],z=mu.wts))
	# abline(v=mu.true[1],h=mu.true[2],col=2)

	# Update parameters of importance function for mu
	q.mu$p <- colSums(mu.wts*rho)	
	q.mu$mean <- colSums((mu.wts*mu)*rho)/q.mu$p
	q.mu$sd <- sqrt(sapply(1:K,function(x) sum(mu.wts*(mu[,x]-q.mu$mean[x])^2*rho[,x])/q.mu$p[x]))
# p.tmp
# q.mu	

	###
	### Update p
	### 
	
	#???
	
	###
	### Update sigma
	###

	# Sample sigma
	sigma <- exp(rnorm(n.samp,log(q.sigma$mean),q.sigma$sd))
	# mean(sigma)
	# hist(sigma,breaks=100)
	
	# Calculate importance weights for sigma
	log.lik <- sapply(sigma,function(x) sum(dnorm.mix(y,q.mu$mean,x,p.true,sum=TRUE,log=TRUE)))
	prior <- dunif(sigma,a,b,log=TRUE)
	q.density <- dnorm(log(sigma),log(q.sigma$mean),q.sigma$sd,log=TRUE)
	sigma.wts <- exp(prior+log.lik-q.density-max(log.lik))
	sigma.wts <- sigma.wts/sum(sigma.wts)

	# Update parameters of importance function for sigma
	q.sigma <- list(mean=sum(sigma.wts*sigma), 
		sd=sqrt(sum(sigma.wts*(sigma-q.sigma$sd)^2)))

}

# Posterior of mu via importance sampling
mu.is <- sample(1:n.samp,10000,mu.wts,replace=TRUE)
mu.is <- c(mu[mu.is,])
hist(mu.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25));abline(v=mu.true,lty=2)
hist(mu.is,breaks=1000,prob=TRUE,col=rgb(1,0,0,0.25),xlim=c(-0.25,0.25));abline(v=mu.true,lty=2)
hist(mu.is,breaks=1000,prob=TRUE,col=rgb(1,0,0,0.25),xlim=c(1.75,2.25));abline(v=mu.true,lty=2)

# Posterior of sigma via importance sampling
sigma.is <- sample(1:n.samp,10000,sigma.wts,replace=TRUE)
sigma.is <- sigma[sigma.is]
hist(sigma.is,breaks=100,prob=TRUE,col=rgb(1,0,0,0.25));abline(v=sigma.true,lty=2)

# Class memberships???
lik <- sapply(1:K,function(x) dnorm(y,q.mu$mean[x],q.sigma$mean))
test <- lik/rowSums(lik)
idx <- apply(test,1,function(x) sample(1:J,1,prob=x))
table(idx)
boxplot(y~idx);abline(h=mu.true,col=2)
plot(y,test[,1]);points(y,test[,2],col=2)



source('~/Documents/personal/coursework/csu/FW 673/biweek 4/code/mixture/finch.MCMC.R')
mcmc.out=finch.MCMC(y,10000,.1,z.known=ifelse(idx==1,1,0))
hist(c(mcmc.out$mu.mat),breaks=1000,prob=TRUE,add=TRUE)
hist(c(mcmc.out$s),breaks=1000,prob=TRUE,add=TRUE,col="gray")