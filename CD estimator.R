
N=1000;n=100
#quantile_calculation
e=rnorm(10000,0,sqrt(0.1))
X=rgamma(10000,shape=1,scale = 2)
Y=X+e
p=ecdf(Y)
t1=Y[which(p(Y)==0.25)];t2=Y[which(p(Y)==0.5)]
t3=Y[which(p(Y)==0.75)];t4=Y[which(p(Y)==1)]
#Iteration process
ITR=function(r,t){
e=rnorm(N,0,sqrt(0.1))
X=rgamma(N,shape=1,scale = 2)
Y=X+e
#Finite population cdf
Q=ecdf(Y)
F_Nt=Q(t)
#CD_estimator & Naive estimator
  s=sample(c(1:1000),n,replace=FALSE)
  y=Y[s];x=X[s]
  x_ns=X[-s]
  beta_est=as.numeric((solve(t(x)%*%x))%*%t(x)%*%y)
  x_ij=matrix(NA,nrow = length(x),ncol=length(x_ns))
  cnt=0;
  for (i in 1:length(x)) {
    for(j in 1:length(x_ns)){
      x_ij[i,j]=x_ns[j]-x[i]
      if(y[i]<=t-(beta_est*x_ij[i,j])){
        cnt=cnt+1
       }
      }
    }
  CD=((length(y[y<=t]))/N)+((cnt)/(N*n))
  F_est=((length(y[y<=t]))/n)
  list(F_Nt=F_Nt,CD=CD,F_est=F_est)
}
R=500
k=c(1:R)
H=lapply((k),t=t1,ITR)
F_Nt=CD=F_est=array(NA,R)
for(i in 1:R){
  F_Nt[i]=H[[i]]$F_Nt
  CD[i]=H[[i]]$CD
  F_est[i]=H[[i]]$F_est
}


#The CD estimator of distribuition function at Quartile 1 is

H[[1]]$CD

#The Naive estimator of distribuition function at Quartile 1 is

H[[1]]$F_est

#The Relative Bias of CD estimator and Naive estimator respectively 

#Relative Bias
mean(CD-F_Nt)/mean(F_Nt);mean(F_est-F_Nt)/mean(F_Nt)
#The Relative MSE of CD estimator and Naive estimator respectively 

# Relative MSE
((sum((CD-F_Nt)^2))/R)/mean(F_Nt);((sum((F_est-F_Nt)^2))/R)/mean(F_Nt)
