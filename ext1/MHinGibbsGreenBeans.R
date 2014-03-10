beta<-c(1,1)
gamma<-c(1,1)

lambda<-function(beta,x){
  exp(beta[1]+ beta[2]*x)
}  
p<-function(gamma,x){
  exp(gamma[1]+gamma[2]*x)/(1+exp(gamma[1]+gamma[2]*x))
}

lhood<-function(y,x,beta,gamma){
  lambda1<-lambda(beta,x)
  p1<-p(gamma,x)
 sum(log(((1-p1+p1*exp(-lambda1))*(y==0)+(p1*(1/factorial(y))*exp(-lambda1)*(lambda1^y))*(y>0))))
}

test_dat<-rbind(head(beans),tail(beans))

lhood(y=test_dat$mvm,x=test_dat$price,gamma=gamma,beta=beta)
#good

prior.beta<-function(beta,Sigma.b,beta0=c(0,0)){
  det(Sigma.b)^(-1/2)*exp((-1/2)*t(beta-beta0)%*%solve(Sigma.b)%*%(beta-beta0))
}
#check this? 
prior.gamma<-function(gamma,Sigma.g,gamma0=c(0,0)){
  det(Sigma.g)^(-1/2)*exp((-1/2)*t(gamma-gamma0)%*%solve(Sigma.g)%*%(gamma-gamma0))
}

#Metropolis Hastings steps for beta, gamma

sample_beta <- function(y,x,betaj_1,sigma.b,Sigma.b,gammaj_1){
  beta_star<-c(0,0)
  beta_star[1] <- rnorm(1, betaj_1[1], sigma.b)
  beta_star[2] <- rnorm(1, betaj_1[2], sigma.b)
  
  f_bstar<-lhood(y,x,beta_star,gammaj_1)*prior.beta(beta_star,Sigma.b)
  f_b<-lhood(y,x,betaj_1,gammaj_1)*prior.beta(betaj_1,Sigma.b)
  
  varmat<-matrix(c(sigma.b,0,0,sigma.b),nrow=2)
  
  q_bstar<-det(varmat)^(-1/2)*exp((-1/2)*t(beta_star-betaj_1)%*%solve(varmat)%*%(beta_star-betaj_1))
  q_b<-det(varmat)^(-1/2)*exp((-1/2)*t(betaj_1-beta_star)%*%solve(varmat)%*%(betaj_1-beta_star))
  
  rho_b<-f_bstar*q_b/(f_b*q_bstar) 
  
  if (is.na(rho_b)==1){
    return(betaj_1)
  }
  
  if (rho_b<1){
    return(beta_star)
  }
  else {return(betaj_1)}
}
###CHECK THIS FUNCTION
sample_gamma <- function(y,x,gammaj_1, sigma.g,Sigma.g,beta_j){
  gamma_star<-c(0,0)
  gamma_star[1] <- rnorm(1, gammaj_1[1], sigma.g)
  gamma_star[2] <- rnorm(1, gammaj_1[2], sigma.g)
  f_gstar<-lhood(y,x,beta_j,gamma_star)*prior.gamma(gamma_star,Sigma.g)
  f_g<-lhood(y,x,beta_j,gammaj_1)*prior.gamma(gammaj_1,Sigma.g)
  
  varmat<-matrix(c(sigma.g,0,0,sigma.g),nrow=2)
  
  q_gstar<-det(varmat)^(-1/2)*exp((-1/2)*t(gamma_star-gammaj_1)%*%solve(varmat)%*%(gamma_star-gammaj_1))
  q_g<-det(varmat)^(-1/2)*exp((-1/2)*t(gammaj_1-gamma_star)%*%solve(varmat)%*%(gammaj_1-gamma_star))
  
  rho_g<-f_gstar*q_g/(f_g*q_gstar) 
  
  if (is.na(rho_g)){
    return(gammaj_1)
  }
  
  if (rho_g<1){
    return(gamma_star)
  }
  else {return(gammaj_1)}
}


#test sampling funtions
sigma.b<-1
Sigma.b<-matrix(c(100,0,0,100),nrow=2)

sample_beta(y=test_dat$mvm,x=test_dat$price,betaj_1=c(1,-1),sigma.b=sigma.b,Sigma.b=Sigma.B,gammaj_1=c(1,-1))

sample_gamma(y=test_dat$mvm,x=test_dat$price,gammaj_1=c(1,-1),sigma.g=sigma.b,Sigma.g=Sigma.B,beta_j=c(1,-1))

#GOOD!!!!!!

#MCMC!!!  

#beta0,gamma0 are initial values
#values for priors for gamma,beta (Sigma) and tuning parameter value (sigma)
sigma<-10
Sigma<-matrix(c(1000,0,0,1000),nrow=2)
beta0<-c(1,1)
gamma0<-c(1,1)

run_mcmc = function(y, x, beta0, gamma0, n.reps=1e3, tune=TRUE, sigma.b=10,sigma.g=10,Sigma=matrix(c(1000,0,0,1000),nrow=2)) {
  beta_keep = matrix(0,ncol=2,nrow=n.reps); beta_keep[1,] = beta = beta0
  gamma_keep =  matrix(0,ncol=2,nrow=n.reps); gamma_keep[1,] = gamma = gamma0
  
  for (i in 1:n.reps) {
    # Automatically tune beta
    beta_old <- beta
    beta <- sample_beta(y,x,beta_old,sigma.b,Sigma,gamma) #Metropolis step
    if (tune==TRUE) {
       if (sum(beta==beta_old)==2) {  #if next draw is rejected then make var smaller  (accept prob too low)
         sigma.b = sigma.b/1.1 
       } else {
         sigma.b = sigma.b*1.1  #increase var to get more of the density sampled  (accept prob too high)
       }
     }
    # Automatically tune gamma
    gamma_old <- gamma
    gamma <- sample_gamma(y,x,gamma,sigma.g,Sigma,beta) #Metropolis step
    if (tune==TRUE) {
      if (sum(gamma==gamma_old)==2) {  #if next draw is rejected then make var smaller  (accept prob too low)
        sigma.g = sigma.g/1.1 
      } else {
        sigma.g = sigma.g*1.1  #increase var to get more of the density sampled  (accept prob too high)
      }
    }
    
    beta_keep[i,]  = beta
    gamma_keep[i,]   = gamma
    }
  return(list(beta=beta_keep,gamma=gamma_keep,sigma.b,sigma.b))  #stores gamma and beta vectors in the same place
}

Sigma.b<-
set.seed(138818)
burnin1 = run_mcmc(y=store1$mvm,x=store1$price,beta0=c(-1,1), gamma0=c(-1,1),sigma.b=1,sigma.g=.5,Sigma=Sigma.b,tune=TRUE)  #3 chains 
burnin2 = run_mcmc(y=store1$mvm,x=store1$price,beta0=c(10,-10), gamma0=c(10,-10),sigma.b=10,sigma.g=10,Sigma=Sigma.b,tune=TRUE)
burnin3 = run_mcmc(y=store1$mvm,x=store1$price,beta0=c(-100,-100), gamma0=c(-100,-100),sigma.b=500,sigma.g=500)
#still not working 
chain1 <- run_mcmc(y=store1$mvm,x=store1$price,beta0=burnin1[[1]][1e3,], gamma0=burnin1[[2]][1e3,],Sigma=Sigma.b,sigma.b=burnin1[[3]],sigma.g=burnin1[[4]],n.reps=1e4,tune=FALSE)
chain2
chain3
