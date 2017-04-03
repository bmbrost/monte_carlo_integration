###
### Douc et al. 2007 - Example 1
###

pmc3norind2=function(nsimu=10000,niter=25,var1=1/3,var2=2/3,var3=1,m1=-1,m2=1,m3=2)
{
v=c(var1,var2,var3)
m=c(m1,m2,m3)

pp=matrix(0,niter,3)
pp[1,]=c(0.05,0.05,0.9)

for (i in 1:niter)
  {
  if (i!=1)
    {
    pp[i,1]=sum(w*(sigma==v[1]))
    pp[i,2]=sum(w*(sigma==v[2]))
    pp[i,3]=sum(w*(sigma==v[3]))
    }
  K=sample(1:3,nsimu,prob=pp[i,],rep=T)
  sigma=v[K]
  mu=m[K]
  xx=rnorm(nsimu,mu,sqrt(sigma))
  lw=log(1/3*dnorm(xx,m[1],sqrt(v[1]))+1/3*dnorm(xx,m[2],sqrt(v[2]))+
  1/3*dnorm(xx,m[3],sqrt(v[3])))-log(pp[i,1]*dnorm(xx,m[1],sqrt(v[1]))+
  pp[i,2]*dnorm(xx,m[2],sqrt(v[2]))+pp[i,3]*dnorm(xx,m[3],sqrt(1)))
  w=exp(lw-max(lw))
  w=w/sum(w)
  }
list(p=pp,x=xx,w=w)
}

res=pmc3norind2()
plot(res$p[,1],ylim=c(0,1),type="l",col="sienna3",lwd=2,ylab="p")
lines(res$p[,1]+res$p[,2],col="sienna3",lw=2)
lines(rep(.333,length(res$p[,1])),col="steelblue2",lty=2)
lines(rep(.666,length(res$p[,1])),col="steelblue2",lty=2)


###
### Douc et al. 2007 - Example 2
###

pmc3nor=function(nsimu=5000,niter=25,var1=18,var2=4,var3=.25,pp1=0.33,pp2=0.33,pp3=0.33)
{
v=c(var1,var2,var3)
pp=matrix(0,niter,3)

pp[1,]=c(pp1,pp2,pp3)/(pp1+pp2+pp3)
xx=rt(nsimu,df=10)
lw=dnorm(xx,log=T)-dt(xx,df=10,log=T)
w=exp(lw-max(lw))
w=w/sum(w)
J=sample(1:nsimu,nsimu,prob=w,rep=T)
xxtilde=xx[J]

for (i in 1:niter){

  if (i!=1){
	pp[i,1]=sum(sigmatilde==var1)/nsimu
	pp[i,2]=sum(sigmatilde==var2)/nsimu
	pp[i,3]=sum(sigmatilde==var3)/nsimu
	}

  K=sample(1:3,nsimu,prob=pp[i,],rep=T)
  sigma=v[K]
  T=sum(K>1)
  xx[K==1]=rt(nsimu-T,df=2)+xxtilde[K==1]
  xx[K>1]=rnorm(T,mean=xxtilde[K>1],sd=sqrt(sigma[K>1]))

  lw=dnorm(xx,log=T)-log(pp[i,1]*dt((xx-xxtilde),df=2)+
    pp[i,2]*dnorm(xx,mean=xxtilde,sd=sqrt(var2))+
    pp[i,3]*dnorm(xx,mean=xxtilde,sd=sqrt(var3)))
  w=exp(lw-max(lw))
  w=w/sum(w)
  J=sample(1:nsimu,nsimu,prob=w,rep=T)
  xxtilde=xx[J]
  sigmatilde=sigma[J]
  }
list(p=pp,x=xx,w=w)
}

plotur=function(T=25000){
#Internal 3 component mixture
xx=rnorm(T)
dif=rnorm(T,sd=sqrt(2))#rnorm(T)-xx
f1=dt(dif,df=2.0)
f2=dnorm(dif,sd=2)
f3=dnorm(dif,sd=.5)
d1=f1-f3
d2=f2-f3

D=50
p=seq(.001,.999,length=D)
q=seq(.001,.999,length=D)
int=matrix(0,ncol=D,nrow=D)

ent=mean(log(dnorm(xx)))
for (i in 1:D){
b1=p[i]*d1 + f3
for (j in (1:D)[q<1-p[i]])
  int[i,j]=mean(log( b1+q[j]*d2 ))
  }
int[int!=0]=ent-int[int!=0]

mi=max(int)
int[int==0]=mi
image(p,q,int,xlab=expression(alpha[1]),ylab=expression(alpha[2]))
contour(p,q,int,nlevels=300,add=T)
}

plotur()
res=pmc3nor(niter=250,pp1=.2,pp2=.25,pp3=.55)
lines(res$p[,1:2],pch=19,cex=.2,col="blue4")
res=pmc3nor(niter=250,pp1=.6,pp2=.05,pp3=.35)
lines(res$p[,1:2],pch=19,cex=.2,col="yellow")
res=pmc3nor(niter=250,pp1=.05,pp2=.05,pp3=.9)
lines(res$p[,1:2],pch=19,cex=.2,col="white")
res=pmc3nor(niter=250,pp1=.49,pp2=.49,pp3=.02)
lines(res$p[,1:2],pch=19,cex=.2,col="green2")




###
### Douc et al. 2007 - Example 3
###

like=function(par,ob){
# log likelihood function
  ob[1,1]*par[2]-exp(par[2])+
  ob[1,2]*par[3]-exp(par[3])+
  ob[2,1]*(par[1]+par[2])-exp(par[1]+par[2])+
  ob[2,2]*(par[1]+par[3])-exp(par[1]+par[3])
  }

pmc4=function(nsimu=50000,niter=5,nvar=10,ranJ=c(-10,10),pp0=(rep(1,nvar)/nvar),grafison=T)
# nsimu is the number of IS particles
# niter is the number of PMC iterations
# nvar is the number of kernels
# ranJ is the range of the variances
# pp0 is the vector of weights
# grafison stands for graph is `on'

{
# Rash model
obs=matrix(c(60,364,36,240),ncol=2,byrow=T)

v=sum(obs)*exp(seq(ranJ[1],ranJ[2],length=nvar)) #choice of variances
pp=matrix(pp0,nrow=niter+1,ncol=nvar,byrow=T)

#Fish info
library(mgcv)
library(MASS)
mle=c(-0.43,4.06,5.9)
Ii=matrix(c(403,37.8,238,37.8,96,0,238,0,604),ncol=3)
SqI=ginv(mroot(Ii))

# Sample of parameters
pareto=matrix(rnorm(3*nsimu,sd=.001),ncol=3)
parms=matrix(mle,ncol=3,nrow=nsimu,byrow=T)+1000*t(SqI%*%t(pareto))
lw=apply(parms,1,like,ob=obs)-apply(dnorm(pareto,log=T,sd=.001),1,sum)

# Importance weights
w=exp(lw-max(lw))
w=w/sum(w)
partilde=parms[sample(1:nsimu,nsimu,prob=w,rep=T),]

if (grafison){
  postscript(file="iter0.eps")
  nems=c(expression(alpha[1]),expression(beta[0]),expression(beta[1]))
  par(mfrow=c(2,2),mar=c(4,2,4,2))
  for (j in 1:3)
    hist(partilde[,j],nclass=135,col="gold",main=nems[j],xlab="")
  hist(apply(partilde,1,like,ob=obs),col="wheat",main="likelihood",xlab="")
  par(new=T)
  plot(sort(w),type="l",col="sienna4",axes=F,xlab="",ylab="",lwd=2)
  dev.off()
  }

if (niter>0){

for (i in 1:niter){

  K=sample(1:nvar,nsimu,prob=pp[i,],rep=T)
  sigma=v[K]
  # normals
  dife=matrix(rnorm(3*nsimu),ncol=3)
  parms=partilde+sqrt(sigma)*t(SqI%*%t(dife))
  dife=sigma*dife^2
  # Rao Blackwellisation
  rbl=pp[i,1]*exp(-0.5*apply(dife,1,sum)/v[1])/v[1]^(1.5)
  for (t in 2:nvar)
   rbl=rbl+pp[i,t]*exp(-0.5*apply(dife,1,sum)/v[t])/v[t]^(1.5)

  lw=apply(parms,1,like,ob=obs)-log(rbl)
  w=exp(lw-max(lw))
  w=w/sum(w)
  partilde=parms[sample(1:nsimu,nsimu,prob=w,rep=T),]

  #Graphics
  if (grafison){
    file=paste("iter",i,".eps",sep="")
    postscript(file=file)
    par(mfrow=c(2,2),mar=c(4,2,4,2))
    for (j in 1:3)
      hist(partilde[,j],nclass=135,col="gold",main=nems[j],xlab="")
    hist(apply(partilde,1,like,ob=obs),col="wheat",main="likelihood & weights",xlab="")
    par(new=T)
    plot(sort(w),type="l",col="sienna4",axes=F,xlab="",ylab="",lwd=2)
    dev.off()
    }

  # More R&B !
  for (j in 1:nvar)
          pp[(i+1),j]=sum(w[K==j])

  pp[(i+1),] = (pp[(i+1),]) #+.0001)/1.001
  pp[(i+1),]
  }
}
list(p=pp,sig=v,is=w,llike=apply(partilde,1,like,ob=obs),par=parms)
}

exploit=function(res){
# Graphical representation of the pmc4 output
# postscript(file="sumacon.eps")

par(mfrow=c(3,3),mar=c(4,2,4,2))

#subsample resampled
mle=c(-0.43,4.06,5.9)
selecto=sample(1:length(res$is),5000,prob=res$is,rep=T)
pur=res$par[selecto,]

# Rash model
obs=matrix(c(60,364,36,240),ncol=2,byrow=T)

hist(pur[,1],nclass=250,col="gold3",xlab=expression(alpha[1]),main="",ylab="",proba=T)
hist(pur[,2],nclass=250,col="gold3",xlab=expression(beta[0]),main="",ylab="",proba=T)
hist(pur[,3],nclass=250,col="gold3",xlab=expression(beta[1]),main="",ylab="",proba=T)

#likelihood sliz
alpha1=seq(min(pur[,1]),max(pur[,1]),length=75)
beta0=seq(min(pur[,2]),max(pur[,2]),length=75)
beta1=seq(min(pur[,3]),max(pur[,3]),length=75)

sliz1=alpha1%*%t(beta0)
par1=matrix(0,ncol=3,nrow=75)
par1[,1]=alpha1
par1[,3]=mle[3]
for (i in 1:75){
 par1[,2]=beta0[i]
 sliz1[,i]=apply(par1,1,like,ob=obs)
}
image(alpha1,beta0,sliz1,xlab=expression(alpha[1]),ylab=expression(beta[0]),col = heat.colors(123))
points(pur[,1],pur[,2],cex=.3,pch=19)

sliz2=sliz1
par2=par1
par2[,1]=alpha1
par2[,2]=mle[2]
for (i in 1:75){
 par2[,3]=beta1[i]
 sliz2[,i]=apply(par2,1,like,ob=obs)
}
image(alpha1,beta1,sliz2,xlab=expression(alpha[1]),ylab=expression(beta[1]),col = heat.colors(123))
points(pur[,1],pur[,3],cex=.3,pch=19)

sliz3=sliz1
par3=par1
par3[,1]=mle[1]
par3[,2]=beta0
for (i in 1:75){
 par3[,3]=beta1[i]
 sliz3[,i]=apply(par3,1,like,ob=obs)
}
image(beta0,beta1,sliz3,xlab=expression(beta[0]),ylab=expression(beta[1]),col = heat.colors(123))
points(pur[,2],pur[,3],cex=.3,pch=19)
#likelihood dots
cala=heat.colors(123)[round(123*(max(res$llike[selecto])-res$llike[selecto])/diff(range(res$llike[selecto])))]

plot(pur[1,1],pur[1,2],cex=.3,pch=19,col=cala[1],xlim=range(alpha1),ylim=range(beta0),
xlab=expression(alpha[1]),ylab=expression(beta[0]))
for (i in 2:5000)
 points(pur[i,1],pur[i,2],cex=.3,pch=19,col=cala[i])

plot(pur[1,1],pur[132],cex=.3,pch=19,col=cala[1],xlim=range(alpha1),ylim=range(beta1),
xlab=expression(alpha[1]),ylab=expression(beta[1]))
for (i in 2:5000)
 points(pur[i,1],pur[i,3],cex=.3,pch=19,col=cala[i])

plot(pur[1,2],pur[1,3],cex=.3,pch=19,col=cala[1],xlim=range(beta0),ylim=range(beta1),
xlab=expression(beta[0]),ylab=expression(beta[1]))
for (i in 2:5000)
 points(pur[i,2],pur[i,3],cex=.3,pch=19,col=cala[i])
#dev.off()
}
