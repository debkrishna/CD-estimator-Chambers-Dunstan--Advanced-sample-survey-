x1=iris$Sepal.Length[51:65];x2=iris$Sepal.Width[51:65]
x3=iris$Petal.Length[51:65];x4=iris$Petal.Width[51:65]
x0=c(rep(0.01,4))
sigma=cov(cbind(x1,x2,x3,x4))
# function gives largest eigen value and  corresponding eigen vector
eigen_fn=function(s){
  d=1;ic=0
  while( ic<10000 && d>0.000001){
    x_old=(s%*%x0)/((sum((x0)^2))^0.5)
    x_new=(s%*%x_old)/((sum((x_old)^2))^0.5)
    d=(sum((x_new-x_old)^2))^0.5
    x_old=x_new
    ic=ic+1
  }
  beta=x_old/(sum((x_old)^2)^0.5)
  lambda=(sum((x_old)^2)^0.5)
  list(beta=beta,lambda=lambda)
}
#function creates new sigma 
nsigma=function(ss){
  ss-((eigen_fn(ss)$lambda)*((eigen_fn(ss)$beta)%*%t((eigen_fn(ss)$beta))))
}
#calculating eigen values and eigen vectors
s1=sigma;iic=0
while(iic<5){
  s2=nsigma(s1)
  print(eigen_fn(s1))
  s1=s2
  iic=iic+1
}
