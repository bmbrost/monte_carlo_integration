guassian.mcmc <- function(y,theta,tau,a,b,start,tune,n.mcmc,adapt=TRUE){
# browser()
	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
			# a <- min(0.01,1/sqrt(k))
			a <- min(0.025,1/sqrt(k))
			exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	n <- length(y)
	mu.save <- rep(0,n.mcmc)
	sigma.save <- rep(0,n.mcmc)
	mu <- start$mu
	sigma <- start$sigma

	keep <- list(sigma=0)
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	for(k in 1:n.mcmc){

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			tune <- get.tune(tune,keep.tmp$sigma,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	
# browser()	
		#  Sample sigma
		sigma.star <- rnorm(1,sigma,tune)
		if(sigma.star>a & sigma.star<b){
		  	mh.0 <- sum(dnorm(y,mu,sigma,log=TRUE))
			mh.star <- sum(dnorm(y,mu,sigma.star,log=TRUE))
			if(exp(mh.star-mh.0)>runif(1)){
				sigma <- sigma.star
				keep$sigma <- keep$sigma+1
				keep.tmp$sigma <- keep.tmp$sigma+1
			}
		}

 		#  Sample mu 
		A.inv <- (n/sigma^2+1/tau^2)^-1
		bb <- sum(y)/sigma^2
		mu <- rnorm(1,A.inv*bb,sqrt(A.inv))

		#  Save Samples 
		mu.save[k] <- mu
  		sigma.save[k] <- sigma

	}

	list(mu=mu.save,sigma.save=sigma.save,n.mcmc=n.mcmc,keep=keep$sigma/n.mcmc,tune=tune)

}
