library(MASS) #For mvrnorm
library(statmod) #For rinvgauss
library(glmnet) #For setting initial Lasso estimates
library(lars) #Has Diabetes data
library(Matrix) #for de-sparsing glmnet coefficients
library(coda)



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



sample_lambdas<-function(p,r,tau2s,delta){
 #Generates a sample from the lambda posterior
  return(sqrt(rgamma(1,shape=p+r,rate=sum(tau2s)/2+delta)))
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

  subsample=FALSE
  if(subsample==TRUE){
    samples=samples*5
  }

 #Create output/saving arrays
  p=dim(x)[2]
  betas=matrix(0,nrow=samples+burnin,ncol=p)
  sigma2s=matrix(0,nrow=samples+burnin,ncol=1)
  tau2s=matrix(0,nrow=samples+burnin,ncol=p)
 
 #Initialize starting values
  start=glmnet(x,y,lambda=c(lambda))
  betas[1,]=Matrix(start$beta)[,1]
  sigma2s[1,]=var(y-x%*%betas[1,])
  tau2s[1,]=1/sqrt(lambda^2*sigma2s[1,]/betas[1,]^2+lambda^2)
  tau2s[1,][tau2s[1,]==0]=min(tau2s[1,][tau2s[1,]!=0])

 #Initialize wild starting values
  betas[1,]=matrix(rep(1000000,10))
  sigma2s[1,]=1000000
  tau2s[1,]=matrix(rep(1000000,10))

 #Run Gibbs sampler
  for(i in 2:(samples+burnin)){
    betas[i,]=sample_betas(y,x,sigma2s[i-1,],tau2s[i-1,])
    sigma2s[i,]=sample_sigma2(y,x,betas[i,],tau2s[i-1,])
    tau2s[i,]=sample_tau2s(betas[i,],sigma2s[i,],lambda)
  }

  inds=c((burnin+1):(samples+burnin))

  if(subsample==TRUE){
    inds=inds[seq(1,length(inds),5)]
  } 

  #inds=c(1:(burnin+samples))

  return(list("betas"=betas[inds,],"sigma2s"=sigma2s[inds,],"tau2s"=tau2s[inds,]))
}



GibbsSampleBadStart<-function(y,x,lambda,samples=10000,burnin=1000){
 #y       - Array of length n
 #x       - Matrix of dimensions nxp
 #lambda  - 
 #samples - Number of (non-burn in) Gibbs samples
 #burnin  - Count of burn in samples

  subsample=FALSE
  if(subsample==TRUE){
    samples=samples*5
  }

 #Create output/saving arrays
  p=dim(x)[2]
  betas=matrix(0,nrow=samples+burnin,ncol=p)
  sigma2s=matrix(0,nrow=samples+burnin,ncol=1)
  tau2s=matrix(0,nrow=samples+burnin,ncol=p)
 
 #Initialize starting values
  #start=glmnet(x,y,lambda=c(lambda))
  #betas[1,]=Matrix(start$beta)[,1]
  #sigma2s[1,]=var(y-x%*%betas[1,])
  #tau2s[1,]=1/sqrt(lambda^2*sigma2s[1,]/betas[1,]^2+lambda^2)
  #tau2s[1,][tau2s[1,]==0]=min(tau2s[1,][tau2s[1,]!=0])

 #Initialize wild starting values
  betas[1,]=matrix(rep(1000000,10))
  sigma2s[1,]=1000000
  tau2s[1,]=matrix(rep(1000000,10))

 #Run Gibbs sampler
  for(i in 2:(samples+burnin)){
    betas[i,]=sample_betas(y,x,sigma2s[i-1,],tau2s[i-1,])
    sigma2s[i,]=sample_sigma2(y,x,betas[i,],tau2s[i-1,])
    tau2s[i,]=sample_tau2s(betas[i,],sigma2s[i,],lambda)
  }

  inds=c((burnin+1):(samples+burnin))

  if(subsample==TRUE){
    inds=inds[seq(1,length(inds),5)]
  } 

  #inds=c(1:(burnin+samples))

  return(list("betas"=betas[inds,],"sigma2s"=sigma2s[inds,],"tau2s"=tau2s[inds,]))
}



hyperpriorGibbsSample<-function(y,x,r=1,delta=1.78,samples=10000,burnin=1000){
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
  lambdas=matrix(0,nrow=samples+burnin,ncol=1)
 
 #Initialize starting values
  lambdas[1]=sqrt(r/delta)
  start=glmnet(x,y,lambda=c(lambdas[1]))
  betas[1,]=Matrix(start$beta)[,1]
  sigma2s[1,]=var(y-x%*%betas[1,])
  tau2s[1,]=1/sqrt(lambdas[1]^2*sigma2s[1,]/betas[1,]^2+lambdas[1]^2)
  tau2s[1,][tau2s[1,]==0]=min(tau2s[1,][tau2s[1,]!=0])

 #Run Gibbs sampler
  for(i in 2:(samples+burnin)){
    betas[i,]=sample_betas(y,x,sigma2s[i-1,],tau2s[i-1,])
    sigma2s[i,]=sample_sigma2(y,x,betas[i,],tau2s[i-1,])
    tau2s[i,]=sample_tau2s(betas[i,],sigma2s[i,],lambdas[i-1,])
    lambdas[i,]=sample_lambdas(p,r,tau2s[i,],delta)
  }
 
  inds=c((burnin+1):(samples+burnin))
  #inds=c(1:(samples+burnin))
  return(list("betas"=betas[inds,],"sigma2s"=sigma2s[inds,],"tau2s"=tau2s[inds,],"lambdas"=lambdas[inds,]))
}



MMLGibbsSample<-function(y,x,kmax=100,samples=10000,burnin=1000){
  #Uses MML to find an estimate for lambda
  #NOTE - REQUIRES OLS ESTIMATES (n>p, t(X)%*%X full rank, etc etc)
  lambdas=rep(0,kmax)  

  p=dim(x)[2]
  ols=lm(y~x)
  lambdas[1]=p*sqrt(var(ols$residuals))/sum(abs(ols$coefficients[2:(1+p)]))
  
  for(k in 2:kmax){
    print(k)
    gsample=GibbsSample(y,x,lambdas[k-1],samples,burnin)
    lambdas[k]=sqrt(2*p/sum(apply(gsample$tau2s,2,mean)))
  }
  gsample=GibbsSample(y,x,mean(lambdas[round(kmax/2):kmax]))
  
  return(list("betas"=gsample$betas,"sigma2s"=gsample$sigma2s,"tau2s"=gsample$tau2s, "lambdas"=lambdas))
}



makeFig2<-function(y,x,lambda=0.237,samples=10000,burnin=1000,save=FALSE){
  if(save){pdf('SaveFigure2.pdf')}

 #Prepare OLS estimates
  p=dim(x)[2]
  ols=lm(y~x)
  
 #Prepare Bayesian Lasso estimates
  gsample=GibbsSample(y,x,lambda,samples,burnin)
  sorted=apply(gsample$betas,2,sort)
  low=floor(.025*samples); high=ceiling(.975*samples);

 #Prepare n-fold Lasso estimates
  nfold=cv.glmnet(x,y,nfolds=dim(x)[1])
  nfoldcoefs=Matrix(nfold$glmnet.fit$beta[,which(nfold$lambda==nfold$lambda.min)])[,1]  

 #Prepare L1-matched Lasso estimates
  L1norm=sum(abs(apply(gsample$betas,2,median)))
  #find closest fit
  nfoldnorms=apply(abs(nfold$glmnet.fit$beta),2,sum)
  d1=nfoldnorms-L1norm
  d2=L1norm-nfoldnorms
  lp1=nfold$glmnet.fit$lambda[d1>0][which.min(d1[d1>0])]
  lp2=nfold$glmnet.fit$lambda[d2>0][which.min(d2[d2>0])]
  #run second Lasso
  l1=min(lp1,lp2)
  l2=max(lp1,lp2,l1+.1)
  L1fitlassos=glmnet(x,y,lambda=seq(l1,l2,length.out=1000))
  L1fitnorms=apply(abs(L1fitlassos$beta),2,sum)
  d3=L1fitnorms-L1norm
  L1fitbetas=Matrix(L1fitlassos$beta[,which.min(abs(d3))])[,1]

 #Set up the plot 
  plot(c(min(c(sorted[low,],ols$coefficients[2:(1+p)])),max(c(sorted[high,],ols$coefficients[2:(1+p)]))),c(1,p)/2,type='n',xlab="Standardized Coefficients",ylab="Variable Number",yaxt='n')
  abline(v=0,lty=3)

 #Plot Bayesian Lasso estimates
  axis(2,at=c(1:p)/2, labels=c(p:1))
  for(i in 1:p){
    lines(sorted[c(low,high),i],c(p+1-i,p+1-i)/2)
    lines(sorted[c(low,low),i],c(p+1-i-.1,p+1-i+.1)/2)
    lines(sorted[c(high,high),i],c(p+1-i-.1,p+1-i+.1)/2)
    points(median(sorted[,i]),(p+1-i)/2,cex=1.2,pch=10)
  }
 
 #Plot OLS estimates
  points(ols$coefficients[2:(1+p)],c(p:1)/2,cex=2,pch=4)

 #Plot Lasso n-fold
  points(nfoldcoefs,(c(p:1)-.2)/2,pch=2)
 
 #Plot Lasso matched L1
  points(L1fitbetas,(c(p:1)+.2)/2,pch=6)

  if(save){dev.off()}
}



makeFig1<-function(y,x,samples=10000,burnin=1000,save=FALSE){
  if(save){pdf('SaveFigure1.pdf')}
  k=40
  lambdamax=45
  lambdas=1.6^c((1+20):(k+20))/1.6^(k+20)*lambdamax
    
 #Take Lasso over range of lambdas
  lbetas=t(Matrix(glmnet(x,y,lambda=lambdas)$beta))
  lloc=apply(abs(lbetas),1,sum)/max(apply(abs(lbetas),1,sum))
  lorder=sort.int(lloc,index.return=TRUE)$ix    

 #Take Ridge over range of lambdas
  rridge=lm.ridge(y~0+x,lambda=lambdas*13000)
  rbetas=t(Matrix(rridge$coef/rridge$scales))
  rloc=apply(abs(rbetas),1,sum)/max(apply(abs(rbetas),1,sum))
  rorder=sort.int(rloc,index.return=TRUE)$ix    

 #Take Gibbs sample over range of lambdas
  gbetas=matrix(0,nrow=k,ncol=dim(x)[2])
  for(i in 1:k){
    gsample=GibbsSample(y,x,lambdas[i],samples,burnin)
    gbetas[i,]=apply(gsample$betas,2,median)
  }
  gloc=apply(abs(gbetas),1,sum)/max(apply(abs(gbetas),1,sum))
  gorder=sort.int(gloc,index.return=TRUE)$ix

 #MML and n-fold setup
  gsample=GibbsSample(y,x,.237,samples,burnin)
  MMLloc=sum(abs(apply(gsample$betas,2,median)))/max(apply(abs(gbetas),1,sum))
  nfold=cv.glmnet(x,y,nfolds=dim(x)[1])
  nfoldbetas=nfold$glmnet.fit$beta[nfold$glmnet.fit$lambda==nfold$lambda.min]
  nfoldloc=sum(abs(nfoldbetas))/max(apply(abs(lbetas),1,sum))

 #Make three plots
  par(mfrow=c(1,3))
  inds=list(c(10,3,7),c(9,1,4,2),c(8,6,5))
  for(j in 1:3){
   #Initialize plot
    if(j==1){
      plot(c(0.03,.97),c(max(max(lbetas),max(gbetas),max(rbetas)),min(min(lbetas),min(gbetas),min(rbetas))),type='n',ylab="Bayesian Lasso",xlab=expression(abs(abs(hat(beta)))[1]/max(abs(abs(hat(beta)))[1])) )
      legend('topleft',lty=c(1,2,3),legend=c("Bayesian Estimates","Lasso Estimates","Ridge Estimates"))
    } else { 
      plot(c(0.03,.97),c(max(max(lbetas),max(gbetas),max(rbetas)),min(min(lbetas),min(gbetas),min(rbetas))),type='n',ylab='',xlab=expression(abs(abs(hat(beta)))[1]/max(abs(abs(hat(beta)))[1])) )
    }
    abline(h=0)

   #Plot Lasso samples
    for(i in inds[[j]]){
      lines(lloc[lorder],lbetas[lorder,i],lty=2)
    }

   #Plot Ridge samples
    for(i in inds[[j]]){
      lines(rloc[rorder],rbetas[rorder,i],lty=3)
    }

   #Plot Gibbs samples
    for(i in inds[[j]]){
      lines(gloc[gorder],gbetas[gorder,i])
      axis(4,at=gbetas[gorder[k],i],labels=c(i))
    }
    
   #Plot MML line and n-fold line
    if(j==1){
      lines(c(MMLloc,MMLloc),c(-2000,660))
      lines(c(nfoldloc,nfoldloc),c(-2000,660),lty=2)
    } else {
      abline(v=MMLloc)
      abline(v=nfoldloc,lty=2)
    }
  }	

  if(save){dev.off()}
}



lambdaConvergencePlot<-function(y,x,kmax,samples=10000,burnin=1000,save=FALSE){
  if(save){pdf('SaveLambdaConverge.pdf')}
  test=MMLGibbsSample(y,x,kmax,samples,burnin)
  plot(test$lambda,type='l',xlab='Iteration',ylab=expression(lambda^(k)))
  abline(h=.237,lty=2)
  abline(h=median(test$lambda[100:kmax]))
  print(median(test$lambda[100:kmax]))
  if(save){dev.off()}
}



hyperpriorCI<-function(y,x,r=1,delta=1.78,samples=10000,burnin=1000){
  sink("HyperpriorCI.txt")
  test=hyperpriorGibbsSample(y,x,r=1,delta=1.78,samples=10000,burnin=1000)
  print(quantile(test$lambdas,c(.025,.975)))
  print(median(test$lambdas))
  sink()
}



printEffectiveSizes<-function(y,x,lambda=.237){
  sink("EffectiveSizes.txt")
  for(i in 1:10){
    print(effectiveSize(mcmc(GibbsSample(y.diab,x.diab,lambda=.237)$betas)))
  }
  sink()
}



beta5ACFPlot<-function(y,x,lambda=.237){
  pdf('Beta5ACF.pdf')
  Beta5=GibbsSample(y.diab,x.diab,lambda)$betas[,5]
  acf(Beta5)
  dev.off()
}



genLogPosterior<-function(b){
  sigma=density(b)$bw
  logpost=rep(0,length(b))
  for(i in 1:length(logpost)){
    logpost[i]=mean(dnorm(b[i],b,sigma))
  }
  return(logpost)
}

makeFig2(y.diab,x.diab)

test=GibbsSample(y.diab,x.diab,lambda=.237,samples=10000,burnin=1000)
test2=GibbsSampleBadStart(y.diab,x.diab,lambda=.237,samples=10000,burnin=1000)

#testh=hyperpriorGibbsSample(y.diab,x.diab,samples=10000,burnin=1000)

#hyperpriorCI(y.diab,x.diab,r=1,delta=1.78,samples=10000,burnin=1000)
#printEffectiveSizes(y.diab,x.diab)
#beta5ACFPlot(y.diab,x.diab)
#lambdaConvergencePlot(y.diab,x.diab,kmax=1000,samples=10000,burnin=1000,save=TRUE)
#makeFig1(y.diab,x.diab,samples=10000,burnin=1000,save=TRUE)
#makeFig2(y.diab,x.diab,lambda=0.237,samples=10000,burnin=1000,save=TRUE)

#RUN ALL THESE THEN SET TO SUBSAMPLE
