#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
using namespace Rcpp;
using namespace arma;


double likelihoodfunction(mat y,mat X,mat X2,mat beta_cur,mat eta_cur,mat sigma,int time,int j,double lambda,int k,IntegerVector z_cur ){
  double pi=3.14;
  double sum=0;
  for(int i=0;i<time;i++){
    sum+=pow((y(j,i)-det(X.row(j)*trans(beta_cur.row(k)))-det(X2.row(j)*trans(eta_cur))),2);
  }
  double likelihood=exp(-sum/(2*pow(sigma[k],2)))/pow(sqrt(2*pi*sigma[k]*sigma[k]),time);
  return likelihood;
}

double residue(mat y,mat X,mat X2,mat beta_cur,mat eta_cur,int time,int j){
  double sum=0;
  for(int i=0;i<time;i++){
    sum+=y(j,i)-det(X2.row(j)*trans(eta_cur));
  }
  return sum;
}

mat residue1(mat y_0,mat X_0,mat X2,mat beta_cur,mat eta_cur,int time,int j){
  mat A(time,1,fill::zeros);
  for(int i=0;i<time;i++){
    A.submat(i,0,i,0)=y_0(i,0)-det(X2.row(j)*trans(eta_cur));
  }
  return A;
}

mat residue2(mat y_hat,mat X_hat,mat X2_hat,mat beta_0_cur,mat eta_cur,int time,int N_grid){
  mat A(N_grid*time,1,fill::zeros);
  for(int i=0;i<N_grid*time;i++){
    A.submat(i,0,i,0)=y_hat(i,0)-det(X2_hat.row(i)*trans(eta_cur));
    
  }
  return A;
}

mat residue3(mat y_1,mat X4,mat X2_hat,mat beta_cur,mat eta_cur,int time,int N_grid,IntegerVector z_cur){
  mat A(N_grid*time,1,fill::zeros);
  for(int i=0;i<N_grid;i++){
    for(int p=0;p<time;p++){
      A.submat(i*time+p,0,i*time+p,0)=y_1(i*time+p,0)-det(X4.row(i*time+p)*trans(beta_cur.row(z_cur[i])));
    }
  }
  return A;
}

// [[Rcpp::export]]
List fn_updateZ(IntegerVector z_cur,mat beta_cur,mat beta_0,mat eta_cur,mat eta_0,int Dim,NumericVector c,NumericVector sigma,int N_grid,mat y,mat X,double k_0,double a,double b,int ClusterNumber,int time,mat X2,double lambda,mat Sigma0,NumericVector VN,mat USState,int Dim2) {
  double pi=3.14;
  double a_hat,b_hat;
  double gamma=1;
  mat A,B,C;
  NumericVector P(N_grid+10);
  NumericVector random,random1;
  for(int j=0;j<N_grid;j++){
    int Cluster=ClusterNumber;
    //Step 1:Compute the probability   
    for(int k=0;k<Cluster;k++){ 
      int Neighbor=0;
      for(int state=0;state<N_grid;state++){
        if(USState(state,j)==1&&z_cur[state]==z_cur[j]){
          Neighbor++;
        }
      }
      int Neighbor2=0;
      for(int state=0;state<N_grid;state++){
        if(USState(state,j)==2&&z_cur[state]==z_cur[j]){
          Neighbor2++;
        }
      }
      //For existing c
      if(z_cur[j]==k){
        //*exp(Neighbor*lambda)
        P[k]=(c[k]+gamma-1)*exp(Neighbor*lambda+Neighbor2*lambda)*likelihoodfunction(y,X,X2,beta_cur,eta_cur,sigma,time,j,lambda,k,z_cur);
      }else {
        P[k]=(c[k]+gamma)*exp(Neighbor*lambda+Neighbor2*lambda)*likelihoodfunction(y,X,X2,beta_cur,eta_cur,sigma,time,j,lambda,k,z_cur);
      }
    }
    //For new c
    mat X_0(time,Dim,fill::zeros),y_0(time,1,fill::zeros);
    for(int p=0;p<time;p++){
      X_0.submat(p,0,p,Dim-1)=X.submat(j,0,j,Dim-1);
      y_0.submat(p,0,p,0)=y.submat(j,p,j,p);
    }
    a_hat=a+(0.5*time);
    A=inv(trans(X_0)*X_0+k_0*inv(Sigma0))*(trans(trans(residue1(y_0,X_0,X2,beta_cur,eta_cur,time,j))*X_0+k_0*beta_0.row(j)*inv(Sigma0)));
    B=trans(X_0)*X_0+k_0*inv(Sigma0); 
    b_hat=(det(trans(residue1(y_0,X_0,X2,beta_cur,eta_cur,time,j))*(residue1(y_0,X_0,X2,beta_cur,eta_cur,time,j)))+2*b+k_0*(det(beta_0.row(j)*inv(Sigma0)*trans(beta_0.row(j))))-det(trans(A)*(B)*A))/2.0;
    random1=rgamma(1,a_hat,1/b_hat);
    P[Cluster]=gamma*exp(VN[Cluster]-VN[Cluster-1])*pow(b,a)*tgamma(a_hat)*sqrt(det(Sigma0))/(pow(sqrt(2*pi),time)*tgamma(a)*(pow(b_hat,a_hat))*sqrt(det(B)));
    double v=P[Cluster];
    for(int iter=0;iter<Cluster;iter++){
      v+=P[iter];
    }
    double v_=1/v;
    NumericVector s=Rcpp::runif(1,0,1);
    double t=v_*P[0];
    //Step 2:Draw a new value for c
    for(int l=0;l<Cluster;l++){
      //For an existing c
      if(s[0]<t){
        c[z_cur[j]]=c[z_cur[j]]-1;
        z_cur[j]=l;
        c[l]=c[l]+1;
        break;
      }else if(s[0]>t&&l==Cluster-1){
        int flag=0;
        //For a new c
        for(int m=0;m<Cluster;m++){
          if(c[m]==0){  
            int index=z_cur[j];
            c[z_cur[j]]=c[z_cur[j]]-1;
            z_cur[j]=m;
            mat beta_cur_0=beta_cur.row(index);
            random1=rgamma(1,a_hat,1/b_hat);
            sigma[m]=1/sqrt(random1[0]);
            for(int p=0;p<Dim;p++){
              mat V=inv(trans(X_0)*X_0+k_0*inv(Sigma0))*pow(sigma[m],2);
              mat V_rp,V_cp,A_0,V_rcp, beta_cur_1;
              V_rp=V;
              V_cp=V;
              A_0=A;
              V_rp.shed_row(p);
              V_cp.shed_col(p);
              V_rcp= V_rp;
              V_rcp.shed_col(p);
              beta_cur_1=beta_cur_0;
              A_0.shed_row(p);
              beta_cur_1.shed_col(p);
              double mu,Sigma;
              mu=det(A.row(p)+V_cp.row(p)*inv(V_rcp)*trans(beta_cur_1-trans(A_0)));
              Sigma=sqrt(det(V(p,p)-V_cp.row(p)*inv(V_rcp)*V_rp.col(p)));
              random=Rcpp::rnorm(1,mu,Sigma); 
              beta_cur.submat(m,p,m,p)=random[0];
              beta_cur_0.submat(0,p,0,p)=random[0];
            }
            c[m]=1;
            flag=1;
            break;
          }
        }
        if(flag!=1){
          int index=z_cur[j];
          c[z_cur[j]]=c[z_cur[j]]-1;
          z_cur[j]=Cluster;
          random1=rgamma(1,a_hat,1/b_hat);
          sigma[Cluster]=1/sqrt(random1[0]);
          mat beta_cur_0=beta_cur.row(index);
          for(int p=0;p<Dim;p++){
            mat V=inv((trans(X_0)*X_0+k_0*inv(Sigma0)))*pow(sigma[Cluster],2);
            mat V_rp,V_cp,A_0,V_rcp, beta_cur_1;
            V_rp=V;
            V_cp=V;
            A_0=A;
            V_rp.shed_row(p);
            V_cp.shed_col(p);
            V_rcp= V_rp;
            V_rcp.shed_col(p);
            beta_cur_1=beta_cur_0;
            A_0.shed_row(p);
            beta_cur_1.shed_col(p);
            double mu,Sigma;
            mu=det(A.row(p)+V_cp.row(p)*inv(V_rcp)*trans(beta_cur_1-trans(A_0)));
            Sigma=sqrt(det(V(p,p)-V_cp.row(p)*inv(V_rcp)*V_rp.col(p)));
            random=Rcpp::rnorm(1,mu,Sigma); 
            beta_cur.submat(Cluster,p,Cluster,p)=random[0];
            beta_cur_0.submat(0,p,0,p)=random[0];
          }
          c[Cluster]=1;
          ClusterNumber++;
        }
        break;
      }
      else{
        t=t+v_*P[l+1];
      }
      
    }
    
  }
  
  for(int k=0; k<ClusterNumber;k++){
    mat A_hat,B_hat,C_hat;
    mat X_hat(N_grid*time,Dim,fill::zeros),beta_0_cur=beta_0,X2_hat(N_grid*time,Dim2,fill::zeros);
    mat y_hat(N_grid*time,1,fill::zeros);
    mat V(Dim,Dim,fill::zeros);
    double sum_1=0;
    if(c[k]!=0){
      for(int s=0;s<N_grid;s++){
        if(z_cur[s]==k){
          for(int p=0;p<time;p++){
            X_hat.submat(sum_1*time+p,0,sum_1*time+p,Dim-1)=X.submat(s,0,s,Dim-1);
            y_hat.submat(sum_1*time+p,0,sum_1*time+p,0)=y.submat(s,p,s,p);
            X2_hat.submat(sum_1*time+p,0,sum_1*time+p,Dim2-1)=X2.submat(s,0,s,Dim2-1);
          }
          sum_1++;
        }
      }
      A_hat=(inv(trans(X_hat)*X_hat+k_0*inv(Sigma0)))*trans(trans(residue2(y_hat,X_hat,X2_hat,beta_0_cur,eta_cur,time,N_grid))*X_hat+(k_0*beta_0_cur.row(0)));
      B_hat=trans(X_hat)*X_hat+k_0*inv(Sigma0);
      a_hat=a+(sum_1*time/2.0);
      b_hat=(det(trans(residue2(y_hat,X_hat,X2_hat,beta_0_cur,eta_cur,time,N_grid))*residue2(y_hat,X_hat,X2_hat,beta_0_cur,eta_cur,time,N_grid))+2*b+(det(beta_0_cur.row(0)*inv(Sigma0)*trans(beta_0_cur.row(0))))-det(trans(A_hat)*(B_hat)*A_hat))/2.0;
      random=rgamma(1,a_hat,1/b_hat);
      sigma[k]=1/sqrt(random[0]);
      V=inv(trans(X_hat)*X_hat+k_0*inv(Sigma0))*pow(sigma[k],2);
      for(int p=0;p<Dim;p++){
        mat V_rp,V_cp,beta_cur_0,A_0,V_rcp;
        V_rp=V;
        V_cp=V;
        A_0=A_hat;
        V_rp.shed_row(p);
        V_cp.shed_col(p);
        V_rcp= V_rp;
        V_rcp.shed_col(p);
        A_0.shed_row(p);
        beta_cur_0= beta_cur.row(k);
        beta_cur_0.shed_col(p);
        double mu,Sigma;
        mu=det(A_hat.row(p)+V_cp.row(p)*(inv(V_rcp))*trans(beta_cur_0-trans(A_0)));
        Sigma=sqrt(det(V(p,p)-V_cp.row(p)*(inv(V_rcp))*V_rp.col(p)));
        random=Rcpp::rnorm(1,mu,Sigma);
        beta_cur.submat(k,p,k,p)=random[0];
      } 
      
    }
  } 
  
  return List::create(Named("z_cur")=z_cur,Named("c")=c,Named("beta_cur")=beta_cur,Named("sigma")=sigma,Named("ClusterNumber")=ClusterNumber);
}

// [[Rcpp::export]]
mat fn_updateEta(IntegerVector z_cur,mat beta_cur,mat beta_0,mat eta_cur,mat eta_0,int Dim,NumericVector c,NumericVector sigma,int N_grid,mat y,mat X,double k_0,double a,double b,int ClusterNumber,int time,mat X2,double lambda,mat Sigma0,int Dim2){
  double k_1=10;
  mat Sigma1(Dim2,Dim2,fill::eye);
  mat A,V;
  mat X3(N_grid*time,Dim2,fill::zeros),y1(N_grid*time,1,fill::zeros),X4(N_grid*time,Dim,fill::zeros);
  for(int i=0;i<N_grid;i++){
    for(int p=0;p<time;p++){
      X3.submat(i*time+p,0,i*time+p,Dim2-1)=X2.submat(i,0,i,Dim2-1);
      X4.submat(i*time+p,0,i*time+p,Dim-1)=X.submat(i,0,i,Dim-1);
      y1.submat(i*time+p,0,i*time+p,0)=y.submat(i,p,i,p);
    }
  }
  A=inv(trans(X3)*X3+k_1*inv(Sigma1))*(trans(trans(residue3(y1,X4,X3,beta_cur,eta_cur,time,N_grid,z_cur))*X3+k_1*eta_cur));
  V=inv(trans(X3)*X3+k_1*inv(Sigma1));
  for(int p=0;p<Dim2;p++){
    mat V_rp,V_cp,eta_cur_0,A_0,V_rcp;
    V_rp=V;
    V_cp=V;
    A_0=A;
    V_rp.shed_row(p);
    V_cp.shed_col(p);
    V_rcp= V_rp;
    V_rcp.shed_col(p);
    A_0.shed_row(p);
    eta_cur_0=eta_cur;
    eta_cur_0.shed_col(p);
    double mu,Sigma;
    mu=det(A(p)+V_cp.row(p)*inv(V_rcp)*trans(eta_cur_0-trans(A_0)));
    Sigma=sqrt(det(V(p,p)-V_cp.row(p)*(inv(V_rcp))*V_rp.col(p)));
    NumericVector random=Rcpp::rnorm(1,mu,Sigma);
    eta_cur.submat(0,p,0,p)=random[0];
  }
  return eta_cur;
}



