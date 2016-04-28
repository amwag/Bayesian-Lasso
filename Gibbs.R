library(MASS) #For mvrnorm
library(statmod) #For rinvgauss
library(glmnet) #For setting initial Lasso estimates
library(lars)

data(diabetes)
attach(diabetes)
#x.diab=x*1/sd(x[,1])
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

#GENERATE TEST DATA, SCALE + CENTER AS NEEDED
Xraw<-matrix(rnorm(700),ncol=7)
Xraw[,5]=rnorm(100,Xraw[,4])
yuncent=Xraw[,1]+5*Xraw[,2]-5*Xraw[,3]+Xraw[,4]+rnorm(100)
y=yuncent-mean(yuncent)
X=scale(Xraw)

test=glmnet(x.diab,y.diab)

#INITIALIZE PARAMETERS
n=floor(length(test$lambda)/2)
#lambdastart=test$lambda[n]
lambdastart=.237
betastart=array(test$beta[,n])
sigma2start=var(y.diab-x.diab%*%betastart)[1,1]
tau2start=1/sqrt((lambdastart^2)*sigma2start/betastart^2)+lambdastart^2

#These don't seem to be working as intended... Check tau2s in particular
betas=matrix(0,nrow=10000,ncol=10);betas[1,]=betastart
tau2s=matrix(0,nrow=10000,ncol=10);tau2s[1,]=tau2start
sigma2s=matrix(0,nrow=10000,ncol=1);sigma2s[1,]=sigma2start

for(i in 2:10000){
	betas[i,]=sample_betas(y.diab,x.diab,sigma2s[i-1,],tau2s[i-1,])
	sigma2s[i,]=sample_sigma2(y.diab,x.diab,betas[i,],tau2s[i-1,])
	tau2s[i,]=sample_tau2s(betas[i,],sigma2s[i,],lambdastart)
}

hist(betas[,1])#0
hist(betas[,2])#-200
hist(betas[,3])#525
hist(betas[,4])#300
hist(betas[,5])#-100
hist(betas[,6])#0
hist(betas[,7])#-180
hist(betas[,8])#100
hist(betas[,9])#500
hist(betas[,10])#70

print("Stop pressing ctrl-r")
