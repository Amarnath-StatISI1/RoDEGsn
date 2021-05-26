library(sn)
#Objective function for estimating MLE of SN Distribution
obj_func_mle=function(data,mu,sigma,gamma){
  data=as.vector(data)
  mu=as.numeric(mu)
  sigma=as.numeric(sigma)
  gamma=as.numeric(gamma)
  n=length(data)
  l=(n*log(2/sigma))-((n/2)*log(2*pi))-(0.5*sum(((data-mu)/sigma)^2))+sum(log(pnorm((gamma*((data-mu)/sigma)),0,1)))
  return(-l)
}


#gradient function for estimating MLE

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

#Calculating the mle using gradient descent

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

#In the following function, we have implemented grad_descenet_mle() function to find MLE for particular data.

sn_mle=function(data_sn)
{
  data=as.vector(data_sn)
  n=length(data_sn)
  dp.est=cp2dp(cp=sn.mple(y=data_sn,opt.method="SANN")$cp,family="SN")
  ans_mle=grad_descent_mle(data_sn,dp.est[1],dp.est[2],dp.est[3],0.04)
  return(as.vector(ans_mle))
}

#objective function for finding DPD estimate at a particular \alpha

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

#Gradient function for dpd

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


#Finding the MDPD Estimator using gradient-descent

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

#In the following function, we have used grad_descent_dpd() for finding DPD estimate when the given data and
#\alpha will be the only output.

sn_dpd=function(data_sn,alpha)
{
  data_sn=as.vector(data_sn)
  alpha=as.numeric(alpha)
  n=length(data_sn)
  dp.est=cp2dp(cp=sn.mple(y=data_sn,opt.method="SANN")$cp,family="SN")
  ans_dpd=grad_descent_dpd(data_sn,dp.est[1],dp.est[2],dp.est[3],0.04,alpha)
  return(as.vector(ans_dpd))
}

##'Calculation of asymptotic variance in skew normal DPD
##'
##'Calculating the \eta vector described in paper
##' @param 
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

#Calculating the J Matrix
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

#Calculating the K Matrix
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

#Here we have estimated the SE of the DPD estimate of mu at a particular \alpha. It will be the [1,1] th 
#element of $J^{-1}KJ^{-1}$.
est_var=function(mu,sigma,gamma,alpha){
  mu=as.numeric(mu)
  gamma=as.numeric(gamma)
  sigma=as.numeric(sigma)
  alpha=as.numeric(alpha)
  var_mat=solve(J(mu,sigma,gamma,alpha))%*%K(mu,sigma,gamma,alpha)%*%solve(J(mu,sigma,gamma,alpha))
  diag(var_mat)=abs(diag(var_mat))
  return(as.matrix(var_mat))
}


# Here we have calculated dg Statistic as given in the paper for both mle and DPD cases.For doing so, we first removed NA
# from data of each group.

#dg statistic calculation formula using mle estimate

dg_mle=function(g1,g2){
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
  sn_d=((m*n)/(m+n))*(t(dp.est1-dp.est2)%*%solve(0.5*(var_g1+var_g2))%*%(dp.est1-dp.est2))
  p_value=pchisq(sn_d,df=3,lower.tail=FALSE)
  return(p_value)
}

#dg statistic calculation formula using dpd estimate

dg_dpd=function(g1,g2,alpha){
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
  sn_d=((m*n)/(m+n))*(t(dp.est1-dp.est2)%*%solve(0.5*(var_g1+var_g2))%*%(dp.est1-dp.est2))
  p_value=pchisq(sn_d,df=1,lower.tail=FALSE)
  return(p_value)
}
