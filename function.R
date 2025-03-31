UpdateZ <- function(this.sample, y,X, a,b,lambda,gamma, N_grid,Dim,beta_0,k_0,time,X2,eta_0,Sigma0,VN,USState,Dim2){
  new.sample<-fn_updateZ(this.sample$z_cur, this.sample$beta_cur, beta_0,this.sample$eta_cur,eta_0,Dim,this.sample$c,this.sample$sigma,N_grid,y,X,k_0,a,b,this.sample$ClusterNumber,time,X2,lambda,Sigma0,VN,USState,Dim2)
  this.sample$z_cur<-new.sample$z_cur
  this.sample$beta_cur<-new.sample$beta_cur
  this.sample$c<-new.sample$c
  this.sample$sigma<-new.sample$sigma
  this.sample$ClusterNumber<-new.sample$ClusterNumber
  return(this.sample)
}
UpdateEta<-function(this.sample,y,X, a,b,lambda,gamma, N_grid,Dim,beta_0,k_0,time,X2,eta_0,Sigma0,Sigma1,Dim2){
  this.sample$eta_cur<-fn_updateEta(this.sample$z_cur, this.sample$beta_cur, beta_0,this.sample$eta_cur,eta_0,Dim,this.sample$c,this.sample$sigma,N_grid,y,X,k_0,a,b,this.sample$ClusterNumber,time,X2,lambda,Sigma0,Dim2)
  return(this.sample$eta_cur)
}

orthogonal_comp<-function(X){
  q=dim(X)[1]
  n=dim(X)[2]
  t=Matrix::rankMatrix(X)[1]
  temp=svd(X,n,n)
  Fmatrix=temp$u[,1:t]%*%diag(temp$d[1:t])
  Qmatrix=t(temp$v)[1:t,]
  temp=svd(Qmatrix,n,n)
  Q_final=rbind(Qmatrix,temp$v[,((t+1):n)])
  M=solve(Q_final)
  M_1=M[,1:t]
  return(M_1)
}

helm <- function(n) {
  mat <- matrix(0, n - 1, n) 
  i <- 2:n
  r <- 1 / sqrt( i * (i - 1) )
  for ( j in 1:(n - 1 ) )   mat[j, 1: c(j + 1) ] <- c( rep(r[j], j), - j * r[j]  ) 
  mat
}
