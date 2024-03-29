library(sn)
##' Set of functions required for two sample MDPDE based robust Wald-type test
##' in case of Skew Normal(SN) distributions
##' ##' @description  T^{\alpha}_{n1,n2} Statistic as given in the paper for both 
##' mle and DPD cases. For doing so, we first removed NA from data of each 
##' group.
##' 
##' @title T^{0}_{n1,n2} statistic calculation formula using mle estimate
##' @param g1: observations from group 1
##' @param g2: observations from group 2
##' @description Takes g1 and g2 as input and return the p-value of the two 
##' sample test
T_mle=function(g1,g2){
  g1=as.vector(g1)
  g2=as.vector(g2)
  g1=g1[!is.na(g1)]
  g2=g2[!is.na(g2)]
  m=length(g1)
  n=length(g2)
  dp.est1=sn_mle(g1)
  var_g1=est_var(dp.est1[1],dp.est1[2],dp.est1[3],0)
  dp.est2=sn_mle(g2)
  var_g2=est_var(dp.est2[1],dp.est2[2],dp.est2[3],0)
  w=m/(m+n)
  u=1-w
  sn_d=((m*n)/(m+n))*(t(dp.est1-dp.est2)%*%solve((w*var_g1)+(u*var_g2))%*%(dp.est1-dp.est2))
  p_value=pchisq(sn_d,df=3,lower.tail=FALSE)
  return(p_value)
}

##' @title T^{\alpha}_{n1,n2} statistic calculation formula using DPD estimate
##' @param g1: observations from group 1
##' @param g2: observations from group 2
##' @param alpha: the tuning parameter
##' @description Takes g1, g2 and $\alpha$ as input and return the p-value of 
##' the two sample test
T_dpd=function(g1,g2,alpha){
  g1=as.vector(g1)
  g2=as.vector(g2)
  g1=g1[!is.na(g1)]
  g2=g2[!is.na(g2)]
  alpha=as.numeric(alpha)
  m=length(g1)
  n=length(g2)
  dp.est1=sn_dpd(g1,alpha)
  var_g1=est_var(dp.est1[1],dp.est1[2],dp.est1[3],alpha)
  dp.est2=sn_dpd(g2,alpha)
  var_g2=est_var(dp.est2[1],dp.est2[2],dp.est2[3],alpha)
  w=m/(m+n)
  u=1-w
  sn_d=((m*n)/(m+n))*(t(dp.est1-dp.est2)%*%solve((w*var_g1)+(u*var_g2))%*%(dp.est1-dp.est2))
  p_value=pchisq(sn_d,df=3,lower.tail=FALSE)
  return(p_value)
}
##'
##' @decsription Functions for finding MLE of SN distribution
##' @title Likelihood function for estimating MLE
##' @param data: the set of observations
##' @param mu: the location parameter of SN distribution
##' @param sigma: the scale parameter of SN distribution
##' @param gamma: the shape parameter of SN distribution
obj_func_mle=function(data,mu,sigma,gamma){
  data=as.vector(data)
  mu=as.numeric(mu)
  sigma=as.numeric(sigma)
  gamma=as.numeric(gamma)
  n=length(data)
  l=(n*log(2/sigma))-((n/2)*log(2*pi))-(0.5*sum(((data-mu)/sigma)^2))+sum(log(pnorm((gamma*((data-mu)/sigma)),0,1)))
  return(-l)
}
##' @title Gradient function of the likelihood
##' @param data: the set of observations
##' @param mu: the location parameter of SN distribution
##' @param sigma: the scale parameter of SN distribution
##' @param gamma: the shape parameter of SN distribution
gradient_func_mle=function(data,mu,sigma,gamma){
  data=as.vector(data)
  mu=as.numeric(mu)
  sigma=as.numeric(sigma)
  gamma=as.numeric(gamma)
  n=length(data)
  g1=sum((data-mu)/(sigma^2))-((gamma/sigma)*sum(dnorm((gamma*((data-mu)/sigma)),0,1)/pnorm((gamma*((data-mu)/sigma)),0,1)))
  g2=-(n/sigma)+((1/(sigma^3))*sum((data-mu)^2))-((gamma/(sigma^2))*sum((data-mu)*(dnorm((gamma*((data-mu)/sigma)),0,1)/pnorm((gamma*((data-mu)/sigma)),0,1))))
  g3=sum(((data-mu)/sigma)*(dnorm((gamma*((data-mu)/sigma)),0,1)/pnorm((gamma*((data-mu)/sigma)),0,1)))
  g=c(g1,g2,g3)
  return(as.vector(-g))
}
##' @title Maximization of obj_func_mle() using gradient descent
##' @param data: the set of observations
##' @param mu_ini: the initial value of location parameter of SN distribution
##' @param sigma_ini: the initial value of scale parameter of SN distribution
##' @param gamma_ini: the initial value of shape parameter of SN distribution
##' @param lambda: the step size constant for gradient descent algorithm
grad_descent_mle=function(data,mu_ini,sigma_ini,gamma_ini,lambda){
  data=as.vector(data)
  mu_ini=as.numeric(mu_ini)
  sigma_ini=as.numeric(sigma_ini)
  gamma_ini=as.numeric(gamma_ini)
  lambda=as.numeric(lambda)
  theta_ini=c(mu_ini,sigma_ini,gamma_ini)
  diff=1
  while(diff>0)
  {
    func_ini=obj_func_mle(data,theta_ini[1],abs(theta_ini[2]),theta_ini[3])
    theta_new=theta_ini-(c(lambda)*c(gradient_func_mle(data,mu_ini,sigma_ini,gamma_ini)))
    func_new=obj_func_mle(data,theta_new[1],abs(theta_new[2]),theta_new[3])
    diff=func_ini-func_new
    theta_ini=theta_new
    mu_ini=theta_ini[1]
    sigma_ini=abs(theta_ini[2])
    gamma_ini=theta_ini[3]
  }
  return(as.vector(theta_new))
}
##'
##' @title Maximum likelihood estimates of SN parameters
##' @param  data_sn: The data in which SN distribution is to be fitted
##' @description It takes data as input and returns MLE of the parameters
##' of the SN distribution
sn_mle=function(data_sn)
{
  data=as.vector(data_sn)
  n=length(data_sn)
  dp.est=cp2dp(cp=sn.mple(y=data_sn,opt.method="SANN")$cp,family="SN")
  ans_mle=grad_descent_mle(data_sn,dp.est[1],dp.est[2],dp.est[3],0.04)
  return(as.vector(ans_mle))
}
##' @decsription Functions for finding DPD estimate of SN distribution
##' @title The DPD function which is to be minimized
##' @param data: the set of observations
##' @param mu: the location parameter of SN distribution
##' @param sigma: the scale parameter of SN distribution
##' @param gamma: the shape parameter of SN distribution
##' @param alpha: the tuning parameter
obj_func_dpd=function(data,mu,sigma,gamma,alpha){
  data=as.vector(data)
  mu=as.numeric(mu)
  sigma=as.numeric(sigma)
  gamma=as.numeric(gamma)
  alpha=as.numeric(alpha)
  y=rsn(100000,mu,sigma,gamma)
  t=mean((dsn(y,mu,sigma))^alpha)-((1+(1/alpha))*(mean((dsn(data,mu,sigma,gamma))^alpha)))
  return(t)
}
##' @title Gradient function of the DPD function
##' @param data: the set of observations
##' @param mu: the location parameter of SN distribution
##' @param sigma: the scale parameter of SN distribution
##' @param gamma: the shape parameter of SN distribution
##' @param alpha: the tuning parameter
gradient_func_dpd=function(data,mu,sigma,gamma,alpha){
  data=as.vector(data)
  mu=as.numeric(mu)
  sigma=as.numeric(sigma)
  gamma=as.numeric(gamma)
  alpha=as.numeric(alpha)
  x=rsn(100000,mu,sigma,gamma)
  g11=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(((x-mu)/sigma^2)-((gamma/sigma)*(dnorm((gamma*((x-mu)/sigma)),0,1)/pnorm((gamma*((x-mu)/sigma)),0,1)))))
  g12=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(-(1/sigma)+(((x-mu)^2)/(sigma^3))-((x-mu)*(gamma/sigma^2)*(dnorm((gamma*((x-mu)/sigma)),0,1)/pnorm((gamma*((x-mu)/sigma)),0,1)))))
  g13=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(dnorm((gamma*((x-mu)/sigma)),0,1)/pnorm((gamma*((x-mu)/sigma)),0,1))*((x-mu)/sigma))
  g1=c(g11,g12,g13)
  g21=mean(((dsn(data,mu,sigma,gamma))^(alpha))*(((data-mu)/sigma^2)-((gamma/sigma)*(dnorm((gamma*((data-mu)/sigma)),0,1)/pnorm((gamma*((data-mu)/sigma)),0,1)))))
  g22=mean(((dsn(data,mu,sigma,gamma))^(alpha))*(-(1/sigma)+(((data-mu)^2)/(sigma^3))-((data-mu)*(gamma/sigma^2)*(dnorm((gamma*((data-mu)/sigma)),0,1)/pnorm((gamma*((data-mu)/sigma)),0,1)))))
  g23=mean(((dsn(data,mu,sigma,gamma))^(alpha))*(dnorm((gamma*((data-mu)/sigma)),0,1)/pnorm((gamma*((data-mu)/sigma)),0,1))*((data-mu)/sigma))
  g2=c(g21,g22,g23)
  g=g1-g2
  return(as.vector(g))
}
##' @title Maximization of obj_func_dpd() using gradient descent
##' @param data: the set of observations
##' @param mu_ini: the initial value of location parameter of SN distribution
##' @param sigma_ini: the initial value of scale parameter of SN distribution
##' @param gamma_ini: the initial value of shape parameter of SN distribution
##' @param lambda: the step size constant for gradient descent algorithm
##' @param alpha: the tuning parameter
grad_descent_dpd=function(data,mu_ini,sigma_ini,gamma_ini,lambda,alpha){
  data=as.vector(data)
  mu_ini=as.numeric(mu_ini)
  sigma_ini=as.numeric(sigma_ini)
  gamma_ini=as.numeric(gamma_ini)
  alpha=as.numeric(alpha)
  lambda=as.numeric(lambda)
  theta_ini=c(mu_ini,sigma_ini,gamma_ini)
  diff=1
  while(diff>0)
  {
    obj_func_ini=obj_func_dpd(data,theta_ini[1],abs(theta_ini[2]),theta_ini[3],alpha)
    theta_new=theta_ini-(c(lambda)*c(gradient_func_dpd(data,mu_ini,sigma_ini,gamma_ini,alpha)))
    obj_func_new=obj_func_dpd(data,theta_new[1],abs(theta_new[2]),theta_new[3],alpha)
    diff=obj_func_ini-obj_func_new
    theta_ini=theta_new
    mu_ini=theta_ini[1]
    sigma_ini=abs(theta_ini[2])
    gamma_ini=theta_ini[3]
  }
  return(as.vector(theta_new))
}
##'
##' @title Maximum likelihood estimates of SN parameters
##' @param data_sn: The data in which SN distribution is to be fitted
##' @param alpha: the tuning parameter
##' @description It takes data and alpha as input and returns the DPD estimate
##'  of the parameters of SN distribution
sn_dpd=function(data_sn,alpha)
{
  data_sn=as.vector(data_sn)
  alpha=as.numeric(alpha)
  n=length(data_sn)
  dp.est=cp2dp(cp=sn.mple(y=data_sn,opt.method="SANN")$cp,family="SN")
  ans_dpd=grad_descent_dpd(data_sn,dp.est[1],dp.est[2],dp.est[3],0.04,alpha)
  return(as.vector(ans_dpd))
}
##'
##' @title Calculating the \eta vector described in paper
##' @param mu is the location parameter of SN distribution
##' @param sigma is the scale parameter
##' @param gamma is the shape parameter
##' @param alpha is the tuning parameter
eta=function(mu,sigma,gamma,alpha)
{
  alpha=as.numeric(alpha)
  mu=as.numeric(mu)
  gamma=as.numeric(gamma)
  sigma=as.numeric(sigma)
  x=rsn(100000,mu,sigma,gamma,alpha)
  z1=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(((x-mu)/(sigma^2))-((gamma/sigma)*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1)))))
  z2=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(-(1/sigma)+(((x-mu)^2)/(sigma^3))-(gamma*((x-mu)/(sigma^2))*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1)))))
  z3=mean(((dsn(x,mu,sigma,gamma))^(alpha))*((x-mu)/sigma)*(dnorm(((x-mu)*(gamma/sigma)),0,1)/pnorm(((x-mu)*(gamma/sigma)),0,1)))
  z=c(z1,z2,z3)
  return(z)
}

##' @title Calculating the J matrix described in paper
##' @param mu is the location parameter of SN distribution
##' @param sigma is the scale parameter
##' @param gamma is the skewness parameter
##' @param alpha is the tuning parameter
J=function(mu,sigma,gamma,alpha){
  mu=as.numeric(mu)
  gamma=as.numeric(gamma)
  sigma=as.numeric(sigma)
  alpha=as.numeric(alpha)
  x=rsn(100000,mu,sigma,gamma)
  j11=mean(((dsn(x,mu,sigma,gamma))^(alpha))*((((x-mu)/(sigma^2))-((gamma/sigma)*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))^2))
  j22=mean(((dsn(x,mu,sigma,gamma))^(alpha))*((-(1/sigma)+(((x-mu)^2)/(sigma^3))-(gamma*((x-mu)/(sigma^2))*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))^2))
  j33=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(((x-mu)/sigma)*(dnorm(((x-mu)*(gamma/sigma)),0,1)/pnorm(((x-mu)*(gamma/sigma)),0,1))^2))
  j12=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(((x-mu)/(sigma^2))-((gamma/sigma)*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))*(-(1/sigma)+(((x-mu)^2)/(sigma^3))-(gamma*((x-mu)/(sigma^2))*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1)))))
  j21=j12
  j13=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(((x-mu)/(sigma^2))-((gamma/sigma)*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))*((x-mu)/sigma)*(dnorm(((x-mu)*(gamma/sigma)),0,1)/pnorm(((x-mu)*(gamma/sigma)),0,1)))
  j31=j13
  j23=mean(((dsn(x,mu,sigma,gamma))^(alpha))*(-(1/sigma)+(((x-mu)^2)/(sigma^3))-(gamma*((x-mu)/(sigma^2))*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))*((x-mu)/sigma)*(dnorm(((x-mu)*(gamma/sigma)),0,1)/pnorm(((x-mu)*(gamma/sigma)),0,1)))
  j32=j23
  j=rbind(c(j11,j12,j13),c(j21,j22,j23),c(j31,j32,j33))
  return(j)
}
##'
##' @title Calculating the K matrix described in paper
##' @param mu is the location parameter of SN distribution
##' @param sigma is the scale parameter
##' @param gamma is the skewness parameter
##' @param alpha is the tuning parameter
K=function(mu,sigma,gamma,alpha){
  mu=as.numeric(mu)
  gamma=as.numeric(gamma)
  sigma=as.numeric(sigma)
  alpha=as.numeric(alpha)
  x=rsn(100000,mu,sigma,gamma)
  k11=mean(((dsn(x,mu,sigma,gamma))^((2*alpha)))*((((x-mu)/(sigma^2))-((gamma/sigma)*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))^2))
  k22=mean(((dsn(x,mu,sigma,gamma))^((2*alpha)))*((-(1/sigma)+(((x-mu)^2)/(sigma^3))-(gamma*((x-mu)/(sigma^2))*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))^2))
  k33=mean(((dsn(x,mu,sigma,gamma))^((2*alpha)))*(((x-mu)/sigma)*(dnorm(((x-mu)*(gamma/sigma)),0,1)/pnorm(((x-mu)*(gamma/sigma)),0,1))^2))
  k12=mean(((dsn(x,mu,sigma,gamma))^((2*alpha)))*(((x-mu)/(sigma^2))-((gamma/sigma)*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))*(-(1/sigma)+(((x-mu)^2)/(sigma^3))-(gamma*((x-mu)/(sigma^2))*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1)))))
  k21=k12
  k13=mean(((dsn(x,mu,sigma,gamma))^((2*alpha)))*(((x-mu)/(sigma^2))-((gamma/sigma)*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))*((x-mu)/sigma)*(dnorm(((x-mu)*(gamma/sigma)),0,1)/pnorm(((x-mu)*(gamma/sigma)),0,1)))
  k31=k13
  k23=mean(((dsn(x,mu,sigma,gamma))^((2*alpha)))*(-(1/sigma)+(((x-mu)^2)/(sigma^3))-(gamma*((x-mu)/(sigma^2))*(dnorm(((gamma/sigma)*(x-mu)),0,1)/pnorm(((gamma/sigma)*(x-mu)),0,1))))*((x-mu)/sigma)*(dnorm(((x-mu)*(gamma/sigma)),0,1)/pnorm(((x-mu)*(gamma/sigma)),0,1)))
  k32=k23
  k=rbind(c(k11,k12,k13),c(k21,k22,k23),c(k31,k32,k33))
  k=k-(eta(mu,sigma,gamma,alpha)%*%t(eta(mu,sigma,gamma,alpha)))
  return(k)
}

##' @title Calculating the asymptotic variance covariance matrix for 
##' particular \alpha
##' @description It will be $J^{-1}KJ^{-1}$. The diagonal elements are the 
##' asymptotic variances of the DPD estimates
##' @param mu is the location parameter of SN distribution
##' @param sigma is the scale parameter
##' @param gamma is the skewness parameter
##' @param alpha is the tuning parameter
est_var=function(mu,sigma,gamma,alpha){
  mu=as.numeric(mu)
  gamma=as.numeric(gamma)
  sigma=as.numeric(sigma)
  alpha=as.numeric(alpha)
  var_mat=solve(J(mu,sigma,gamma,alpha))%*%K(mu,sigma,gamma,alpha)%*%solve(J(mu,sigma,gamma,alpha))
  diag(var_mat)=abs(diag(var_mat))
  return(as.matrix(var_mat))
}
