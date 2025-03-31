library(Rcpp)
library(MASS)
library(ggplot2)
library(usmap)
library(fossil)
library(Matrix)
library(LaplacesDemon)
library(GGally)
sourceCpp("MFM.cpp")
source("function.R")
##Load Data
USgdist=readRDS("USgdist.rds")
USState=as.matrix(USgdist)
beta_posterior=array(0,dim=c(100,51,3))
beta_real=matrix(0,51,3)
eta_posterior=array(0,dim=c(100,3))
beta_posterior_upper=array(0,dim=c(100,51,3))
beta_posterior_lower=array(0,dim=c(100,51,3))
eta_posterior_upper=array(0,dim=c(100,3))
eta_posterior_lower=array(0,dim=c(100,3))
z_posterior=matrix(0,100,51)
rd1=rep(0,100)
LPML=matrix(0,100,11)
for(step in 0:10){
  SampleNumber=100
  
  for(Num in 1:SampleNumber){
    set.seed(Num)
    X1=matrix(0,51,3)
    for(i in 1:51){
      X1[i,]<-rdirichlet(1,c(1,3,6))
    }
    X3=matrix(0,51,3)
    for(i in 1:51){
      for(j in 1:3){
        X3[i,j]=log(X1[i,j])
      }
    }
    X2=matrix(0,51,3)
    data2=data.matrix(read.csv("data_others.csv",header=TRUE))
    for(i in 1:51){
      for(j in 1:3){
        X2[i,j]<-runif(1,-1,1)
      }
    }
    g1=data2[,8]
    #Simulate Beta
    beta=matrix(0,3,3)
    beta[1,]=c(1,-2,1)
    beta[2,]=c(-4,-3,7)
    beta[3,]=c(10,-9,-1)
    eta=c(1,2,1)
    #Simulate y
    y=matrix(0,51)
    for(i in 1:51){
      e=rnorm(1,0,1)
      y[i]=X3[i,]%*%beta[g1[i],]+X2[i,]%*%eta+e
    }
    N_grid <- length(y)
    time=1
    #dim
    H=helm(3)
    M=orthogonal_comp(H)
    
    X=X3%*%M
    Dim<-length(X[1,])
    Dim2<-length(X2[1,])
    # initialize z ##########################################
    ClusterNumber<- 1
    z_cur=rep(0,51+10)
    beta_0=matrix(0,51,2)
    eta_0=matrix(0,1,3)
    beta_cur=matrix(0,51+10,2)
    eta_cur=matrix(0,1,3)
    for(i in 1:51){
      beta_0[i,]=rep(0,2)
    }
    eta_0[1,]=c(0,0,0)
    c=matrix(0,51+10)
    c[1]=51
    c[2:61]=rep(0,60)
    # priors for sigma
    a <- 0.01#shape
    b <- 0.01#rate
    k_0<-1
    k_1<-1
    sigma=matrix(0,51+10)
    sigma<-rep(0.1,51+10)
    Sigma0=diag(2)
    Sigma1=diag(3)
    #other parameters 
    lambda<-step*0.5
    Lambda=1
    gamma<-1
    #VN
    N=N_grid ##  is the number of oberservations
    VN=matrix(0,100)
    VN<-0
    tmax = N+10
    for (t in 1:tmax)
    { 
      r=log(0)
      for (k in t:500)
      {
        l = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, Lambda, log = TRUE)
        m = max(l,r)
        r =log(exp(r-m) + exp(l-m)) + m
      }
      VN[t] = r
    }
    # Initialization ###########################
    eta_cur[1,]<-c(0,0,0)
    cur.sample<-NULL
    cur.sample$z_cur<-z_cur
    cur.sample$beta_cur<-beta_cur
    cur.sample$c<-c
    cur.sample$sigma<-sigma
    cur.sample$ClusterNumber<-ClusterNumber
    cur.sample$eta_cur<-eta_cur
    #Do MCMC
    n.MCMC.iter <-1500
    ##not burn.in
    nburn.in <- 500
    
    save.result<- NULL
    save.result$z <- array(0,dim = c(nburn.in, N_grid))
    save.result$beta<-array(0,dim=c(nburn.in,N_grid+10,Dim))
    save.result$c <- array(0,dim = c(nburn.in, N_grid))
    save.result$ClusterNumber<-array(0,dim = c(nburn.in))
    save.result$eta<-array(0,dim=c(nburn.in,Dim2))
    save.result$sigma<-array(0,dim=c(nburn.in,N_grid+10))
    
    for(i.iter in 1:n.MCMC.iter)
    {
      # Update z
      new.sample<-UpdateZ(cur.sample,y,X, a,b,lambda,gamma, N_grid,Dim,beta_0,k_0,time,X2,eta_0,Sigma0,VN,USState,Dim2)
      cur.sample$z_cur<-new.sample$z_cur
      cur.sample$beta_cur<-new.sample$beta_cur
      cur.sample$c<-new.sample$c
      cur.sample$sigma<-new.sample$sigma
      cur.sample$ClusterNumber<-new.sample$ClusterNumber
      cur.sample$eta_cur<-UpdateEta(cur.sample,y,X, a,b,lambda,gamma, N_grid,Dim,beta_0,k_0,time,X2,eta_0,Sigma0,Sigma1,Dim2)
      if(i.iter>1000){
        j=i.iter-1000
        save.result$z[j,]<-cur.sample$z_cur[1:51]
        save.result$c[j,]<-cur.sample$c[1:51]
        save.result$ClusterNumbe[j]<-cur.sample$ClusterNumber
        save.result$beta[j,,]<-cur.sample$beta_cur
        save.result$eta[j,]<-cur.sample$eta_cur
        save.result$sigma[j,]<-cur.sample$sigma
      }
      
    }
    
    CPO=0
    k=0
    for(j in 1:nburn.in){
      L=0
      for(i in 1:51){
        a=1/(exp(-(y[i]-X[i,]%*%save.result$beta[j,1+save.result$z[j,i],]-X2[i,]%*%save.result$eta[j,])^2/(2*(save.result$sigma[j,1+save.result$z[j,i]])^2))/(sqrt(2*pi)*save.result$sigma[j,1+save.result$z[j,i]]))
        L=L+a
      }
      ##remove the iterations with extreme values
      if(L>1e6){
        k=k+1
      }else{
        CPO=CPO+log(1/L)
      }
      
      
    }
    
    LPML[Num,step+1]=CPO/(nburn.in-k)
    
  }
  
}
best_lambda=rep(0,100)
for(i in 1:100){
  best_lambda[i]=(which.max(LPML[i,])-1)*0.5
}
lambda_best=data.frame("lambda"=best_lambda)
write.csv(x=lambda_best,file="best_lambda2.csv")




