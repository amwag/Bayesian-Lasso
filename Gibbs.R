library(MASS) #For mvrnorm
library(statmod) #For rinvgauss
library(glmnet) #For setting initial Lasso estimates
library(lars) #Has Diabetes data
library(Matrix) #for de-sparsing glmnet coefficients

#LOAD DIABETES DATA
data(diabetes)
attach(diabetes)
x.diab=x
y.diab=y-mean(y)
detach(diabetes)

sample_betas<-function(y,X,sigma2,tau2s){
 #Generates a sample from the beta posterior
  A=gen_A(X,tau2s)
  return(mvrnorm(1,mu=solve(A)%*%t(X)%*%y,Sigma=sigma2*solve(A)))
}

sample_sigma2<-function(y,X,betas,tau2s){
 #Generates a sample from the sigma2 posterior
 #shape=(n-1)/2+p/2
  shape=(length(y)-1)/2+length(betas)/2
 #scale=||y-Xbeta||_2/2+betat*Dtau*beta
  scale=t(y-X%*%betas)%*%(y-X%*%betas)/2+t(betas)%*%solve(diag(tau2s))%*%betas/2
  return(1/rgamma(1,shape=shape,scale=1/scale))
}

sample_tau2s<-function(betas,sigma2,lambda){
 #Generates a sample from the tau2 posterior
  means=sqrt(lambda^2*sigma2/betas^2)
  return(1/rinvgauss(length(betas),means,lambda^2))  
}

gen_A<-function(X,tau2s){
 #Generates the matrix A used in posterior distributions
  return(t(X)%*%X+solve(diag(tau2s)))
}

GibbsSample<-function(y,x,lambda,samples=10000,burnin=1000){
 #y       - Array of length n
 #x       - Matrix of dimensions nxp
 #lambda  - 
 #samples - Number of (non-burn in) Gibbs samples
 #burnin  - Count of burn in samples

 #Create output/saving arrays
  p=dim(x)[2]
  betas=matrix(0,nrow=samples+burnin,ncol=p)
  sigma2s=matrix(0,nrow=samples+burnin,ncol=1)
  tau2s=matrix(0,nrow=samples+burnin,ncol=p)
 
 #Initialize starting values
  start=glmnet(x,y,lambda=c(lambda))
  betas[1,]=Matrix(start$beta)[,1]
  sigma2s[1,]=var(y-x%*%betas[1,])
  tau2s[1,]=1/sqrt((lambda^2)*sigma2s[1,]/(betas[1,]^2)+lambda^2)

 #Run Gibbs sampler
  for(i in 2:(samples+burnin)){
    betas[i,]=sample_betas(y,x,sigma2s[i-1,],tau2s[i-1,])
    sigma2s[i,]=sample_sigma2(y,x,betas[i,],tau2s[i-1,])
    tau2s[i,]=sample_tau2s(betas[i,],sigma2s[i,],lambda)
  }
 
 #Remove burn in samples and return
  inds=c((burnin+1):(samples+burnin))
  return(list("betas"=betas[inds,],"sigma2s"=sigma2s[inds,],"tau2s"=tau2s[inds,]))
}

MMLGibbsSample<-function(y,x,kmax=100,samples=10000,burnin=1000){
  #Uses MML to find an estimate for lambda
  #NOTE - REQUIRES OLS ESTIMATES (n>p, t(X)%*%X full rank, etc etc)
  lambdas=rep(0,kmax)  

  p=dim(x)[2]
  ols=lm(y~x)
  lambdas[1]=p*sqrt(var(ols$residuals))/sum(abs(ols$coefficients[2:(1+p)]))
  
  for(k in 2:kmax){
    gsample=GibbsSample(y,x,lambdas[k-1],samples,burnin)
    lambdas[k]=sqrt(2*p/sum(apply(gsample$tau2s,2,mean)))
    print(lambdas[k])
  }
  gsample=GibbsSample(y,x,lambdas[kmax])
  
  return(list("betas"=gsample$betas,"sigma2s"=gsample$sigma2s,"tau2s"=gsample$tau2s, "lambdas"=lambdas))
}

test=MMLGibbsSample(y.diab,x.diab,kmax=50,samples=1000,burnin=500)

sorted=apply(test$betas,2,sort)

#pdf('Figure2Recreation.pdf')
plot(c(min(sorted[501,]),max(sorted[9500,])),c(1,10),type='n',xlab="Standardized Coefficients",ylab="Variable Number",yaxt='n')
axis(2,at=c(1:10), labels=c(10:1))
for(i in 1:10){
      lines(sorted[c(501,9500),i],c(11-i,11-i))
      lines(sorted[c(501,501),i],c(11-i-.1,11-i+.1))
      lines(sorted[c(9500,9500),i],c(11-i-.1,11-i+.1))
      points(sorted[5000,i],11-i,pch='X')
}
abline(v=0)
#dev.off()