library(MASS) #For mvrnorm
library(statmod) #For rinvgauss

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
print("HAS MEAN")
print(scale/(shape-1))
print(scale)
print(shape)
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

#LOAD LIBRARY FOR PARAMETER INITIALIZATION
library(glmnet)
test=glmnet(X,y)

#INITIALIZE PARAMETERS
n=floor(length(test$lambda)/2)
lambdastart=test$lambda[n]
betastart=array(test$beta[,n])
sigma2start=var(y-X%*%betastart)[1,1]
tau2start=1/sqrt((lambdastart^2)*sigma2start/betastart^2)+lambdastart^2


#These don't seem to be working as intended... Check tau2s in particular
wee=sample_betas(y,X,sigma2start,tau2start)
wee2=sample_sigma2(y,X,wee,tau2start)
wee3=sample_tau2s(wee,wee2,lambdastart)