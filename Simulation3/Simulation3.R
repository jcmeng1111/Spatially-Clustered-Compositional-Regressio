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
##Dimension of beta
Dim1=5
##Dimension of eta
Dim2=3
USgdist=readRDS("USgdist.rds")
USState=as.matrix(USgdist)
beta_posterior=array(0,dim=c(100,51,Dim1))
beta_real=matrix(0,51,Dim1)
eta_posterior=array(0,dim=c(100,Dim2))
beta_posterior_upper=array(0,dim=c(100,51,Dim1))
beta_posterior_lower=array(0,dim=c(100,51,Dim1))
eta_posterior_upper=array(0,dim=c(100,Dim2))
eta_posterior_lower=array(0,dim=c(100,Dim2))
lambda_best=matrix(0,100)
lambda_best=as.matrix(read.csv("best_lambda3.csv")$lambda)
rd1=matrix(0,100,6)
Number=matrix(0,100,6)
MMSE_Eta<-matrix(0,Dim2,3)
MSD_Eta<-matrix(0,Dim2,3)
MAB_Eta<-matrix(0,Dim2,3)
MAB<-matrix(0,Dim1,3)
MSD<-matrix(0,Dim1,3)
MMSE<-matrix(0,Dim1,3)
z_posterior=matrix(0,100,51)
Number=matrix(0,100,6)
for(step in 1:3){
  for(SampleNumber in 1:100){
    if(step==1){
      lambda=0
    }else if(step==2){
      lambda=1
    }else if(step==3){
      lambda=lambda_best[SampleNumber]
    }
    set.seed(SampleNumber)
    X1=matrix(0,51,Dim1)
    for(i in 1:51){
      X1[i,]<-rdirichlet(1,c(1,5,8,1,3))
    }
    X3=matrix(0,51,Dim1)
    for(i in 1:51){
      for(j in 1:Dim1){
        X3[i,j]=log(X1[i,j])
      }
    }
    X2=matrix(0,51,Dim2)
    data2=data.matrix(read.csv("data_others.csv",header=TRUE))
    for(i in 1:51){
      for(j in 1:Dim2){
        X2[i,j]<-runif(1,-0.5,0.5)
      }
    }
    g1=data2[,7]
    #Simulate Beta
    beta=matrix(0,3,Dim1)
    beta[1,]=c(2,-1,1,-1,-1)
    beta[2,]=c(3,5,3,-8,-3)
    beta[3,]=c(5,-6,-1,5,-3)
    beta[3,]%*%c(1,5,8,1,3)
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
    H=helm(Dim1)
    M=orthogonal_comp(H)
    
    X=X3%*%M
    Dim<-length(X[1,])
    Dim2<-length(X2[1,])
    # initialize z ##########################################
    ClusterNumber<- 1
    z_cur=rep(0,51+10)
    beta_0=matrix(0,51,Dim1-1)
    eta_0=matrix(0,1,Dim2)
    beta_cur=matrix(0,51+10,Dim1-1)
    eta_cur=matrix(0,1,Dim2)
    for(i in 1:51){
      beta_0[i,]=rep(0,Dim1-1)
    }
    eta_0[1,]=rep(0,Dim2)
    c=matrix(0,51+10)
    c[1]=51
    c[2:61]=rep(0,60)
    # priors for sigma
    a <- 0.01#shape
    b <- 0.01#scale
    
    k_0<-1
    k_1<-1
    sigma<-rep(0.1,51+10)
    Sigma0=diag(Dim1-1)
    Sigma1=diag(Dim2)
    #other parameters 
    gamma<-1
    Lambda=1
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
    eta_cur[1,]<-rep(0,Dim2)
    cur.sample<-NULL
    cur.sample$z_cur<-z_cur
    cur.sample$beta_cur<-beta_cur
    cur.sample$c<-c
    cur.sample$sigma<-sigma
    cur.sample$ClusterNumber<-ClusterNumber
    cur.sample$eta_cur<-eta_cur
    #Do MCMC
    n.MCMC.iter <-1500
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
    
    
    ##Dahl method
    H=array(0,c(nburn.in,51,51))
    for(i in 1:nburn.in){
      for(j in 1:51){
        for(k in 1:51){
          if(save.result$c[i,j]==save.result$c[i,k]){
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
        H_hat[i,j]<-sum(H[1:500,i,j])/500
      }
    }
    min=1
    for(k in 1:500){
      L<-sum((H[min,1:51,1:51]-H_hat[1:51,1:51])^2)
      S<-sum((H[k,1:51,1:51]-H_hat[1:51,1:51])^2)
      if(S<L){
        min=k
      }
    }
    
    ###store computation results
    for(i in 1:51){
      beta_posterior[SampleNumber,i,]<-M%*%save.result$beta[min,save.result$z[min,i]+1,]
    }
    for(i in 1:51){
      z_posterior[SampleNumber,i]<-save.result$z[min,i]
    }
    g2<-save.result$z[min,1:51]
    rd1[SampleNumber,step]<-rand.index(g1,g2)
    beta_result=array(0,c(nburn.in,51,Dim1))
    for(j in 1:nburn.in){
      for(i in 1:51){
        beta_result[j,i,]<-M%*%save.result$beta[j,save.result$z[j,i]+1,]
      }
    }
    eta_posterior[SampleNumber,]<-save.result$eta[min,]
    for(i in 1:51){
      beta_real[i,]<-beta[g1[i],]
    }
    for(i in 1:SampleNumber){
      Number[i,step]=as.character(length(unique(z_posterior[i,])))
    }
  
  }
  
  
  if(step<4){
    #BetaMAB
    AB<-rep(0,Dim1)
    for(m in 1:Dim1){
      for(l in 1:51){
        AB[m]=AB[m]+sum(abs(beta_posterior[1:SampleNumber,l,m]-beta_real[l,m]))/SampleNumber
      }
    }
    MAB[,step]=AB/51
    #BetaMSD
    SD<-rep(0,Dim1)
    for(m in 1:Dim1){
      for(l in 1:51){
        SD[m]=SD[m]+sqrt(sum((beta_posterior[1:SampleNumber,l,m]-sum(beta_posterior[1:SampleNumber,l,m])/SampleNumber)^2)/(SampleNumber-1))
      }
    }
    MSD[,step]=SD/51
    #BetaMMSE
    MSE<-rep(0,Dim1)
    for(m in 1:Dim1){
      for(l in 1:51){
        MSE[m]=MSE[m]+sum((beta_posterior[1:SampleNumber,l,m]-beta_real[l,m])^2)/SampleNumber
      }
    }
    MMSE[,step]<-MSE/51
    #BetaMCR
    CR<-rep(0,Dim1)
    MCR<-rep(0,Dim1)
    for(m in 1:Dim1){
      for(l in 1:51){
        sum=0
        for(p in 1:SampleNumber){
          if(beta_real[l,m]<beta_posterior_upper[p,l,m] && beta_real[l,m]>beta_posterior_lower[p,l,m]){
            sum=sum+1
          }
        }
        CR[m]=CR[m]+sum/SampleNumber
      }
    }
    MCR=CR/51
    for(m in 1:3){
      MAB_Eta[m,step]=sum(abs(eta_posterior[1:SampleNumber,m]-eta[m]))/SampleNumber
    }
    #EtaMSD
    for(m in 1:3){
      MSD_Eta[m,step]=sqrt(sum((eta_posterior[1:SampleNumber,m]-sum(eta_posterior[1:SampleNumber,m])/SampleNumber)^2)/(SampleNumber-1))
    }
    #EtaMMSE
    for(m in 1:3){
      MMSE_Eta[m,step]=sum((eta_posterior[1:SampleNumber,m]-eta[m])^2)/SampleNumber
    }
    #EtaMCR
    MCR_Eta<-rep(0,Dim2)
    for(m in 1:Dim2){
      sum=0
      for(p in 1:SampleNumber){
        if(eta[m]<eta_posterior_upper[p,m] && eta[m]>eta_posterior_lower[p,m]){
          sum=sum+1
        }
        MCR_Eta[m]=sum/SampleNumber
        
      }
    }
  }
  
}


#map
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7","#FF4040","#FFD39B","#00BFC4","#FC717F")
stateinfo <- readRDS("state_fips_name.rds")
c0=matrix(0,11,51)
sum1=1
for(k in 1:cur.sample$ClusterNumber){
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
stateinfo$cluster <- as.factor(0*(stateinfo$abbr %in% c0[1,])+1*(stateinfo$abbr %in% c0[2,])+2*(stateinfo$abbr %in% c0[3,])+3*(stateinfo$abbr %in% c0[4,])+4*(stateinfo$abbr %in% c0[5,])+5*(stateinfo$abbr %in% c0[6,])+6*(stateinfo$abbr %in% c0[7,])+7*(stateinfo$abbr %in% c0[8,])+8*(stateinfo$abbr %in% c0[9,])+9*(stateinfo$abbr %in% c0[10,])+10*(stateinfo$abbr %in% c0[11,]))
plot_usmap(data = stateinfo, values = "cluster", labels = TRUE) +
  scale_fill_manual(values = cbPalette, name = "Cluster")








MAB_1=rep(0,15)
MSD_1=rep(0,15)
MMSE_1=rep(0,15)
MAB_1[1:5]<-MAB[,1]
MAB_1[6:10]<-MAB[,2]
MAB_1[11:15]<-MAB[,3]
MSD_1[1:5]<-MSD[,1]
MSD_1[6:10]<-MSD[,2]
MSD_1[11:15]<-MSD[,3]
MMSE_1[1:5]<-MMSE[,1]
MMSE_1[6:10]<-MMSE[,2]
MMSE_1[11:15]<-MMSE[,3]

Beta=data.frame("beta"=rep(c("1","2","3","4","5"),9),"y"=c(MAB_1,MSD_1,MMSE_1),group1=c(rep("MAB",15),rep("MSD",15),rep("MMSE",15)),group2=rep("design1",45),"lambda"=rep(c(rep("0",5),rep("1",5),rep("LPML",5)),9))
Beta$beta<-factor(Beta$beta,level=c("1","2","3","4","5"))
p1<-ggplot(data=Beta,aes(x=beta,y=y,group=lambda),group=group1)+geom_point()+geom_line(aes(color=lambda))+facet_grid(group1~group2,scales="free")+labs(x=expression(hat(beta)),y="")+theme_bw()+scale_color_manual(values=c("green4","blue4","orange"),labels=c(expression(paste(lambda,"=0")),expression(paste(lambda,"=1")),expression(paste("LPML"))))+theme(legend.position = "top",legend.title = element_blank())
p1

MAB_Eta_1=rep(0,9)
MSD_Eta_1=rep(0,9)
MMSE_Eta_1=rep(0,9)

MAB_Eta_1[1:3]<-MAB_Eta[,1]
MAB_Eta_1[4:6]<-MAB_Eta[,2]
MAB_Eta_1[7:9]<-MAB_Eta[,3]
MSD_Eta_1[1:3]<-MSD_Eta[,1]
MSD_Eta_1[4:6]<-MSD_Eta[,2]
MSD_Eta_1[7:9]<-MSD_Eta[,3]
MMSE_Eta_1[1:3]<-MMSE_Eta[,1]
MMSE_Eta_1[4:6]<-MMSE_Eta[,2]
MMSE_Eta_1[7:9]<-MMSE_Eta[,3]
Eta=data.frame("eta"=rep(c("1","2","3"),9),"y"=c(MAB_Eta_1,MSD_Eta_1,MMSE_Eta_1),group1=c(rep("MAB",9),rep("MSD",9),rep("MMSE",9)),group2=rep("design1",27),"lambda"=rep(c(rep("0",3),rep("5",3),rep("from mDIC",3)),3))
p2<-ggplot(data=Eta,aes(x=eta,y=y,group=lambda),group=group1)+geom_point()+geom_line(aes(color=lambda))+facet_grid(group1~group2,scales="free")+labs(x=expression(eta),y="")+theme_bw()+scale_color_manual(values=c("green4","blue4","orange"),labels=c(expression(paste(lambda,"=0")),expression(paste(lambda,"=1")),expression(paste("LPML"))))+theme(legend.position = "top",legend.title = element_blank())
p2


randindex=rep(0,300)
randindex[1:100]<-rd1[,1]
randindex[101:200]<-rd1[,2]
randindex[201:300]<-rd1[,3]
m1<-round(median(rd1[,1]),4)
m2<-round(median(rd1[,2]),4)
m3<-round(median(rd1[,3]),4)
df1=data.frame("randindex"=randindex,"lambda"=c(rep("0",100),rep("1",100),rep("LPML",100)))
df1$lambda<-factor(df1$lambda,levels=c("0","1","LPML"))
means <- data.frame("lambda"=c(0,1,"LPML"),randindex=c(m1,m2,m3))
p3<-ggplot(df1,aes(y=randindex,x=lambda))+geom_boxplot()+theme_bw()+labs(x=expression(lambda),y="Rand Index")+geom_text(data = means, aes(label = randindex,y=1.01),color="red")
p3




Number_total=rep(0,300)
Number_total[1:100]=Number[,1]
Number_total[101:200]=Number[,2]
Number_total[201:300]=Number[,3]
Number_total=as.numeric(Number_total)
df=data.frame("Number"=Number_total,"lambda"=c(rep("lambda=0",100),rep("lambda=1",100),rep("lambda from DIC",100)))
df$lambda1<-factor(df$lambda,levels=c("lambda=0","lambda=1","lambda from DIC"),labels=c("lambda=0"=expression(paste(lambda,"=0")),"lambda=1"=expression(paste(lambda,"=1")),"lambda from DIC"="LPML"))
df$flag=df$Number==3
p4<-ggplot(df,aes(x=Number,fill=flag))+geom_bar(width=0.5,color="black")+facet_wrap(.~lambda1,labeller = label_parsed)+geom_text(stat='count',aes(label=..count..), color="black", size=3.5,position=position_dodge(0.5),vjust=-0.5)+labs(x="Inferred Number of Clusters",y="Number of Simulation Replicates")+theme_bw()+scale_y_continuous(breaks =c(seq(0,100,25)),limits=c(0,100))+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8),limits=c(0,9))+scale_fill_manual(values=c("dodgerblue4","darkorange1"))+theme(legend.position = "none")
p4
library(gridExtra)
grid.arrange(p4,p3,p1,p2,top="Cluster Setting 1")

