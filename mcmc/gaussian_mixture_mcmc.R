gaussian.mixture.mcmc <- function(y,theta,tau,start,tune,n.mcmc,
	z.truth=NULL,p.truth=NULL,adapt=TRUE){

get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
	# a <- min(0.01,1/sqrt(k))
	a <- min(0.025,1/sqrt(k))
	exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
}

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

# browser()

n <- length(y)
mu.save <- matrix(0,n.mcmc,2)
beta.save <- numeric(n.mcmc)

# p.save <- rep(0,n.mcmc)
# z.mean <- rep(0,n)

mu1 <- start$mu[1]
mu2 <- start$mu[2]
beta <- start$beta
p <- expit(beta)
sigma <- start$sigma

z <- start$z
if(!is.null(z.truth)) z <- z.truth

# p <- start$p
# if(!is.null(p.truth)) p <- p.truth


keep <- list(mu1=0,mu2=0,beta=0)
keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
Tb <- 50  # frequency of adaptive tuning

for(k in 1:n.mcmc){

	if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
		keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
		tune$mu1 <- get.tune(tune$mu1,keep.tmp$mu1,k)
		tune$mu2 <- get.tune(tune$mu2,keep.tmp$mu2,k)
		tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
		keep.tmp <- lapply(keep.tmp,function(x) x*0)
	} 	


	###############################################################################
	### Updates using Metropolis-Hastings and an integrated likelihood
	###############################################################################

# browser()	
	###
	### Update mu1 and mu2
	###
	
	mu1.star <- rnorm(1,mu1,tune$mu1)
	mu2.star <- rnorm(1,mu2,tune$mu2)
	mh.0 <- sum(dnorm.mix(y,c(mu1,mu2),sigma=rep(sigma,2),p=c(p,1-p)))
	mh.star <- sum(dnorm.mix(y,c(mu1.star,mu2.star),sigma=rep(sigma,2),p=c(p,1-p)))	
	# mh.0 <- sum(log(p*dnorm(y,mu1,sigma,log=FALSE)+
		# (1-p)*dnorm(y,mu2,sigma,log=FALSE)))
	# mh.star <- sum(log(p*dnorm(y,mu1.star,sigma,log=FALSE)+
		# (1-p)*dnorm(y,mu2.star,sigma,log=FALSE)))	
	if(exp(mh.star-mh.0)>runif(1)){
		mu1 <- mu1.star
		mu2 <- mu2.star
		keep$mu1 <- keep$mu1+1
		keep$mu2 <- keep$mu2+1
		keep.tmp$mu1 <- keep.tmp$mu1+1
		keep.tmp$mu2 <- keep.tmp$mu2+1
	}
	
	###
	### Update beta
	###

	beta.star <- rnorm(1,beta,tune$beta)	
	p.star <- expit(beta.star)  # binomial probability
  	mh.0 <- sum(dnorm.mix(y,c(mu1,mu2),sigma=rep(sigma,2),p=c(p,1-p)))
	mh.star <- sum(dnorm.mix(y,c(mu1,mu2),sigma=rep(sigma,2),p=c(p.star,1-p.star)))	
  	if(exp(mh.star-mh.0)>runif(1)){
		beta <- beta.star
		p <- p.star
		keep$beta <- keep$beta+1
		keep.tmp$beta <- keep.tmp$beta+1
	}

	

	###############################################################################
	### Update mu1 and mu2 using Gibbs sampling and indicator variables
	###############################################################################

	###
	### Sample z
	### 

	# if(is.null(z.truth)){
		# p.tmp <- p*dnorm(y,mu1,s)/(p*dnorm(y,mu1,s)+(1-p)*dnorm(y,mu2,s))
 		# z <- rbinom(n,1,p.tmp)		
	# }

	# if(!is.null(z.truth)) z <- z.truth

	###
	### Sample mu1
	### 

	# A.inv <- (sum(z==1)/sigma^2+1/tau^2)^-1
	# bb <- sum(y[z==1])/sigma^2
	# mu1 <- rnorm(1,A.inv*bb,sqrt(A.inv))


	###
	### Sample mu1
	### 
	
	# A.inv <- (sum(z==2)/sigma^2+1/tau^2)^-1
	# bb <- sum(y[z==2])/sigma^2
	# mu2 <- rnorm(1,A.inv*bb,sqrt(A.inv))
	

	###
	### Sample p
	### 
	 
	# p=rbeta(1,sum(z)+1,sum(1-z)+1)
 
	###
	###  Save Samples 
 	###

	mu.save[k,] <- c(mu1,mu2)
 	beta.save[k] <- beta
 	# s.vec[k]=s
	# p.vec[k]=p

}

###
### Write output
### 
	
list(mu=mu.save,beta=beta.save,tune=tune,keep=lapply(keep,function(x) x/n.mcmc))

}
