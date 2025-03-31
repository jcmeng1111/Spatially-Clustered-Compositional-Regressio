library(Rcpp)
library(MASS)
library(ggplot2)
library(usmap)
library(fossil)
library(Matrix)
library(LaplacesDemon)
library(GGally)
library(mclust)
sourceCpp("MFM.cpp")
source("function.R")
##Load Data
data=read.csv("data_Industry.csv",header=TRUE,stringsAsFactors = FALSE)
data=data.matrix(data)
USgdist=readRDS("USgdist.rds")
USState=as.matrix(USgdist)
X1=matrix(0,51,16)
X=matrix(0,51,2)
X1[,1:16]=data[,3:18]
X3=matrix(0,51,3)
X3[,1]=X1[,1]
X3[,2]=rowSums(X1[,2:5])
X3[,3]=rowSums(X1[,6:16])
for(i in 1:51){
  for(j in 1:3){
    X3[i,j]<-log(X3[i,j])
  }
}
Group1<-rep(0,51)
Group2<-matrix(0,100,51)
for(sampleN in 1:1){
  for(step in 0:10){
    RD<-rep(0,100)
    set.seed(sampleN+100)
    H=helm(3)
    M=orthogonal_comp(H)
    X[,1:2]=X3%*%M
    data2=data.matrix(read.csv("data_others.csv",header=TRUE))
    X2<-matrix(0,51,1)
    ##Convert percentage to decimal
    X2[,1]=data2[,3]/1e2
    y=scale(as.matrix(data2[,2]))
    N_grid <- length(y)
    time=1
    #dim
    Dim<-length(X[1,])
    Dim2<-length(X2[1,])
    # initialize z ##########################################
    ClusterNumber_upper<- 1
    z_cur=rep(0,51+10)
    z_cur[1:51]<-sample(0:(ClusterNumber_upper-1), size=51, replace=TRUE)
    ClusterNumber<-length(unique(z_cur[1:51]))
    beta_0=matrix(0,51,Dim)
    eta_0=matrix(0,1,Dim2)
    beta_cur=matrix(0,51+10,Dim)
    eta_cur=matrix(0,1,Dim2)
    for(i in 1:51){
      beta_0[i,]=rep(0,Dim)
    }
    eta_0[1,]=rep(0,Dim2)
    c_ini=matrix(0,51+10)
    c_ini[1]=51
    # priors for sigma
    a <- 0.01#shape
    b <- 0.01#rate
    k_0<-1
    k_1<-1
    sigma=matrix(0.01,51+10)
    Sigma0=diag(1,Dim)
    Sigma1=diag(1,Dim2)
    #other parameters 
    Lambda<-1
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
    for(p in 1:Dim){
      beta_cur[1,p]<-0
    }
    
    
    lambda<-0.5*step
    eta_cur[1,]<-rep(0,Dim2)
    cur.sample<-NULL
    cur.sample$z_cur<-z_cur
    cur.sample$beta_cur<-beta_cur
    cur.sample$c<-c_ini
    cur.sample$sigma<-sigma
    cur.sample$ClusterNumber<-ClusterNumber
    cur.sample$eta_cur<-eta_cur
    #Do MCMC
    n.MCMC.iter <-40000
    nburn.in <-20000
    
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
      if(i.iter>20000){
        j=i.iter-20000
        save.result$z[j,]<-cur.sample$z_cur[1:51]
        save.result$c[j,]<-cur.sample$c[1:51]
        save.result$ClusterNumbe[j]<-cur.sample$ClusterNumber
        save.result$beta[j,,]<-cur.sample$beta_cur
        save.result$eta[j,]<-cur.sample$eta_cur
        save.result$sigma[j,]<-cur.sample$sigma
      }
      
    }
    ##Dahl method
    H=array(0,c(nburn.in,51,51))
    for(i in 1:nburn.in){
      for(j in 1:51){
        for(k in 1:51){
          if(save.result$z[i,j]==save.result$z[i,k]){
            H[i,j,k]=1
          }else{
            H[i,j,k]=0
          }
        }
      }
    }
    
    H_hat=matrix(0,51,51)
    for(i in 1:51){
      for(j in 1:51){
        H_hat[i,j]<-sum(H[1:nburn.in,i,j])/nburn.in
      }
    }
    min=1
    for(k in 1:nburn.in){
      L<-0
      S<-0
      for(i in 1:51){
        for(j in 1:51){
          L<-L+(H[min,i,j]-H_hat[i,j])^2
          S<-S+(H[k,i,j]-H_hat[i,j])^2
          
        }
      }
      if(S<L){
        min=k
      }
    }
    
    
    randindex_1=rep(0,nburn.in)
    for(i in 1:nburn.in){
      
      g1<-save.result$z[i,]
      g2<-save.result$z[min,]
      randindex_1[i]=rand.index(g1,g2)
      
    }
    randindex=subset(randindex_1,randindex_1>0)
    rand_index=as.data.frame(randindex)
    ggplot(rand_index,aes(y=randindex))+geom_boxplot()
    
    
    
    #LPML
    CPO=0
    k=0
    for(j in 1:nburn.in){
      L=0
      for(i in 1:51){
        a=1/(exp(-(y[i]-X[i,]%*%save.result$beta[j,1+save.result$z[j,i],]-X2[i,]%*%save.result$eta[j,])^2/(2*(save.result$sigma[j,1+save.result$z[j,i]])^2))/(sqrt(2*pi)*save.result$sigma[j,1+save.result$z[j,i]]))
        L=L+a
      }
      ##remove the iterations with extreme values
      if(L==Inf){
        k=k+1
      }else{
        CPO=CPO+log(1/L)
      }
      
      
    }
    LPML=CPO/(nburn.in-k)
    cat(LPML,",")
    
    
    if(sampleN==1){
      Group1<- save.result$z[min,]
    }
    Group2[sampleN,]<-save.result$z[min,]
    
  }
}





for(i in 1:1){
  for(sampleN in 1:100){
    RD[sampleN]<-rand.index(Group2[i,],Group2[sampleN,])
  }  
  cat(mean(RD[1:100]),"\n")
}
Group2[2,]
RD


cbPalette <- c( "#56B4E9", "#009E73", "#E69F00","#F0E442",
                "#0072B2", "#D55E00", "#CC79A7","#FF4040","#FFD39B","#00BFC4","#FC717F")
stateinfo <- readRDS("state_fips_name.rds")
c0=matrix(0,11,51)
sum1=1
for(k in 1:51){
  sum2=1
  index=0
  for(i in 1:51){
    if(save.result$z[min,i]==k-1){
      c0[sum1,sum2]<-stateinfo$abbr[i]
      sum2=sum2+1
      index=1
    }
  }
  if(index==1){
    sum1=sum1+1
  }
}
stateinfo$cluster <- as.factor(2*(stateinfo$abbr %in% c0[1,])+1*(stateinfo$abbr %in% c0[2,])+3*(stateinfo$abbr %in% c0[3,])+4*(stateinfo$abbr %in% c0[4,])+5*(stateinfo$abbr %in% c0[5,])+6*(stateinfo$abbr %in% c0[6,])+7*(stateinfo$abbr %in% c0[7,])+8*(stateinfo$abbr %in% c0[8,])+9*(stateinfo$abbr %in% c0[9,])+10*(stateinfo$abbr %in% c0[10,])+11*(stateinfo$abbr %in% c0[11,]))
plot_usmap(data = stateinfo, values = "cluster", labels = TRUE)+scale_fill_manual(values = cbPalette, name = "Cluster")+theme(legend.position = "right")



