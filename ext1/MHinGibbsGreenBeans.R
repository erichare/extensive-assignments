beta<-c(1,1)
delta<-c(1,1)

lambda<-function(beta,x){
  exp(beta[1]+ beta[2]*x)
}  
p<-function(delta,x){
  exp(delta[1]+delta[2]*x)/(1+exp(delta[1]+delta[2]*x))
}

lhood<-function(y,x,beta,delta){
  lambda1<-lambda(beta,x)
  p1<-p(delta,x)
  one_obs<-log(( 1 - p1 + p1 * exp(-lambda1))) * (y==0) + (lgamma(y+1) + log(p1) - lambda1 +y*log(lambda1))*(y>0) 
  return(sum(one_obs))
}

beans <- read.csv("~/Desktop/601/Assignments/Extensive Assignments/greenbeans.csv")
test_dat<-rbind(head(beans),tail(beans))

lhood(y=test_dat$mvm,x=test_dat$price,delta=delta,beta=beta)
lhood(y=187,x=.4,delta=delta,beta=beta)
#good

prior.beta<-function(beta,Sigma.b,beta0=c(0,0)){
  log(det(Sigma.b)^(-1/2)*exp((-1/2)*t(beta-beta0)%*%solve(Sigma.b)%*%(beta-beta0)))
}
#check this? 
prior.delta<-function(delta,Sigma.d,delta0=c(0,0)){
  log(det(Sigma.d)^(-1/2)*exp((-1/2)*t(delta-delta0)%*%solve(Sigma.d)%*%(delta-delta0)))
}

#uniform prior on (-50,50)
prior.delta.u12<-function(delta){
  .01 * (delta > -50 && delta < 50)
}
  

#Metropolis Hastings steps for beta, delta
#fix these steps!!  See http://mcmcinirt.stat.cmu.edu/archives/320
sample_beta <- function(y,x,betaj_1,sigma.b,Sigma.b,deltaj_1){
  beta_star<-c(0,0)
  beta_star[1] <- rnorm(1, betaj_1[1], sigma.b)
  beta_star[2] <- rnorm(1, betaj_1[2], sigma.b)
  
  f_bstar<-lhood(y,x,beta_star,deltaj_1)+prior.beta(beta_star,Sigma.b)
  f_b<-lhood(y,x,betaj_1,deltaj_1)+prior.beta(betaj_1,Sigma.b)
  
  varmat<-matrix(c(sigma.b,0,0,sigma.b),nrow=2)
  
  q_bstar<-log(det(varmat)^(-1/2)*exp((-1/2)*t(beta_star-betaj_1)%*%solve(varmat)%*%(beta_star-betaj_1)))
  q_b<-log(det(varmat)^(-1/2)*exp((-1/2)*t(betaj_1-beta_star)%*%solve(varmat)%*%(betaj_1-beta_star)))
  
  log_rho_b<-min((f_bstar - f_b) + (q_b - q_bstar) , 0)
  
  if (log(runif(1)) <= log_rho_b){
    acc.new = 1
  }
   else acc.new = 0
  
  if (acc.new == 1){
    return(beta_star)
  }
  else return(betaj_1)
}
###CHECK THIS FUNCTION
sample_delta <- function(y,x,deltaj_1, sigma.d=1,Sigma.d,beta_j){
  delta_star<-c(0,0)
  delta_star[1] <- rnorm(1, deltaj_1[1], sigma.d)
  delta_star[2] <- rnorm(1, deltaj_1[2], sigma.d)
  f_dstar<-lhood(y,x,beta_j,delta_star)+prior.delta(delta_star,Sigma.d)
  f_d<-lhood(y,x,beta_j,deltaj_1)+prior.delta(deltaj_1,Sigma.d)
  
  varmat.d<-matrix(c(sigma.d,0,0,sigma.d),nrow=2)
  
  q_dstar<-log(det(varmat.d)^(-1/2)*exp((-1/2)*t(delta_star-deltaj_1)%*%solve(varmat.d)%*%(delta_star-deltaj_1)))
  q_d<-log(det(varmat.d)^(-1/2)*exp((-1/2)*t(deltaj_1-delta_star)%*%solve(varmat.d)%*%(deltaj_1-delta_star)))
  
  log_rho_d<-min((f_dstar - f_d) + (q_d - q_dstar),0) 
  
 if (log(runif(1)) <= log_rho_d){
    acc.new.d = 1
  }
  else acc.new.d = 0
  
  if (acc.new.d == 1){
    return(delta_star)
  }
  else return(deltaj_1)
}

#test sampling funtions
sigma.b<-1
Sigma.b<-matrix(c(100,0,0,100),nrow=2)

sample_beta(y=test_dat$mvm,x=test_dat$price,betaj_1=c(1,-1),sigma.b=10,Sigma.b=matrix(c(100,0,0,100),nrow=2),deltaj_1=c(1,-1))

sample_delta(y=test_dat$mvm,x=test_dat$price,deltaj_1=c(1,-1),sigma.d=10,Sigma.d=matrix(c(100,0,0,100),nrow=2),beta_j=c(1,-1))

#GOOD!!!!!!

#MCMC!!!  

#beta0,delta0 are initial values
#values for priors for delta,beta (Sigma) and tuning parameter value (sigma)
sigma<-10
Sigma<-matrix(c(1000,0,0,1000),nrow=2)
beta0<-c(1,1)
delta0<-c(1,1)

run_mcmc = function(y, x, beta0, delta0, n.reps=1e3, tune=TRUE, sigma.b=10,sigma.d=10,Sigma=matrix(c(100,0,0,100),nrow=2)) {
  beta_keep = matrix(0,ncol=2,nrow=n.reps); beta_keep[1,] = beta = beta0
  delta_keep =  matrix(0,ncol=2,nrow=n.reps); delta_keep[1,] = delta = delta0
  
  for (i in 1:n.reps) {
    # Automatically tune beta var
    beta_old <- beta
    beta <- sample_beta(y,x,beta_old,sigma.b,Sigma,delta) #Metropolis step
    if (tune==TRUE) {
       if (sum(beta==beta_old)==2) {  #if next draw is rejected then make var smaller  (accept prob too low)
         sigma.b = sigma.b/1.1 
       } else {
         sigma.b = sigma.b*1.1  #increase var to get more of the density sampled  (accept prob too high)
       }
     }
    # Automatically tune delta var
    delta_old <- delta
    delta <- sample_delta(y,x,delta,sigma.d,Sigma,beta) #Metropolis step
    #delta <- sample_delta_u(y,x,delta,sigma.d,beta) #Metropolis step new prior
    if (tune==TRUE) {
      if (sum(delta==delta_old)==2) {  #if next draw is rejected then make var smaller  (accept prob too low)
        sigma.d = sigma.d/1.1 
      } else {
        sigma.d = sigma.d*1.1  #increase var to get more of the density sampled  (accept prob too high)
      }
    }
    
    beta_keep[i,]  = beta
    delta_keep[i,]   = delta
    }
  return(list(beta=beta_keep,delta=delta_keep,sigma.b,sigma.d))  #stores delta and beta vectors in the same place
}

#store 1009 aka store 1
store1<-subset(beans,store==1009)
burnin1 = run_mcmc(y=store1$mvm,x=store1$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
burnin2 = run_mcmc(y=store1$mvm,x=store1$price,beta0=c(5,-5), delta0=c(5,-5),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
burnin3 = run_mcmc(y=store1$mvm,x=store1$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

chain1 <- run_mcmc(y=store1$mvm,x=store1$price,beta0=burnin1[[1]][1e3,], delta0=burnin1[[2]][1e3,],sigma.b=burnin1[[3]],sigma.d=burnin1[[4]],n.reps=1e4,tune=FALSE)
chain2 <- run_mcmc(y=store1$mvm,x=store1$price,beta0=burnin2[[1]][1e3,], delta0=burnin2[[2]][1e3,],sigma.b=burnin2[[3]],sigma.d=burnin2[[4]],n.reps=1e4,tune=FALSE)
chain3 <- run_mcmc(y=store1$mvm,x=store1$price,beta0=burnin3[[1]][1e3,], delta0=burnin3[[2]][1e3,],sigma.b=burnin3[[3]],sigma.d=burnin3[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store1<-NULL
mcmc_store1$delta0<-c(chain1[[1]][,1],chain2[[1]][,1],chain3[[1]][,1])
mcmc_store1$delta1<-c(chain1[[1]][,2],chain2[[1]][,2],chain3[[1]][,2])
mcmc_store1$beta0<-c(chain1[[2]][,1],chain2[[2]][,1],chain3[[2]][,1])
mcmc_store1$beta1<-c(chain1[[2]][,2],chain2[[2]][,2],chain3[[2]][,2])
mcmc_store1<-data.frame(mcmc_store1)
saveRDS(mcmc_store1,"Store1Samples.RDS")

mcmc_store1<-readRDS("Store1Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store1,main='Store 1, beta_0')+geom_vline(xintercept=mean(mcmc_store1$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store1,main='Store 1, beta_1')+geom_vline(xintercept=mean(mcmc_store1$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store1,main='Store 1, delta_0')+geom_vline(xintercept=mean(mcmc_store1$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store1,main='Store 1, delta_1')+geom_vline(xintercept=mean(mcmc_store1$delta1),color='red')

#store 1010 aka store 2
store2<-subset(beans,store==1010)
bin1_store2 = run_mcmc(y=store2$mvm,x=store2$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store2 = run_mcmc(y=store2$mvm,x=store2$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store2 = run_mcmc(y=store2$mvm,x=store2$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store2 <- run_mcmc(y=store2$mvm,x=store2$price,beta0=bin1_store2[[1]][1e3,], delta0=bin1_store2[[2]][1e3,],sigma.b=bin1_store2[[3]],sigma.d=bin1_store2[[4]],n.reps=1e4,tune=FALSE)
c2_store2 <- run_mcmc(y=store2$mvm,x=store2$price,beta0=bin2_store2[[1]][1e3,], delta0=bin2_store2[[2]][1e3,],sigma.b=bin2_store2[[3]],sigma.d=bin2_store2[[4]],n.reps=1e4,tune=FALSE)
c3_store2 <- run_mcmc(y=store2$mvm,x=store2$price,beta0=bin3_store2[[1]][1e3,], delta0=bin3_store2[[2]][1e3,],sigma.b=bin3_store2[[3]],sigma.d=bin3_store2[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store2<-NULL
mcmc_store2$delta0<-c(c1_store2[[1]][,1],c2_store2[[1]][,1],c3_store2[[1]][,1])
mcmc_store2$delta1<-c(c1_store2[[1]][,2],c2_store2[[1]][,2],c3_store2[[1]][,2])
mcmc_store2$beta0<-c(c1_store2[[2]][,1],c2_store2[[2]][,1],c3_store2[[2]][,1])
mcmc_store2$beta1<-c(c1_store2[[2]][,2],c2_store2[[2]][,2],c3_store2[[2]][,2])
mcmc_store2<-data.frame(mcmc_store2)
saveRDS(mcmc_store2,"store2Samples.RDS")

mcmc_store2<-readRDS("store2Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store2,main='Store 2, beta_0')+geom_vline(xintercept=mean(mcmc_store2$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store2,main='Store 2, beta_1')+geom_vline(xintercept=mean(mcmc_store2$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store2,main='Store 2, delta_0')+geom_vline(xintercept=mean(mcmc_store2$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store2,main='Store 2, delta_1')+geom_vline(xintercept=mean(mcmc_store2$delta1),color='red')

#store 1011 aka store 3
store3<-subset(beans,store==1011)
bin1_store3 = run_mcmc(y=store3$mvm,x=store3$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store3 = run_mcmc(y=store3$mvm,x=store3$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store3 = run_mcmc(y=store3$mvm,x=store3$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store3 <- run_mcmc(y=store3$mvm,x=store3$price,beta0=bin1_store3[[1]][1e3,], delta0=bin1_store3[[2]][1e3,],sigma.b=bin1_store3[[3]],sigma.d=bin1_store3[[4]],n.reps=1e4,tune=FALSE)
c2_store3 <- run_mcmc(y=store3$mvm,x=store3$price,beta0=bin2_store3[[1]][1e3,], delta0=bin2_store3[[2]][1e3,],sigma.b=bin2_store3[[3]],sigma.d=bin2_store3[[4]],n.reps=1e4,tune=FALSE)
c3_store3 <- run_mcmc(y=store3$mvm,x=store3$price,beta0=bin3_store3[[1]][1e3,], delta0=bin3_store3[[2]][1e3,],sigma.b=bin3_store3[[3]],sigma.d=bin3_store3[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store3<-NULL
mcmc_store3$delta0<-c(c1_store3[[1]][,1],c2_store3[[1]][,1],c3_store3[[1]][,1])
mcmc_store3$delta1<-c(c1_store3[[1]][,2],c2_store3[[1]][,2],c3_store3[[1]][,2])
mcmc_store3$beta0<-c(c1_store3[[2]][,1],c2_store3[[2]][,1],c3_store3[[2]][,1])
mcmc_store3$beta1<-c(c1_store3[[2]][,2],c2_store3[[2]][,2],c3_store3[[2]][,2])
mcmc_store3<-data.frame(mcmc_store3)
saveRDS(mcmc_store3,"store3Samples.RDS")

mcmc_store3<-readRDS("store3Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store3,main='Store 3, beta_0')+geom_vline(xintercept=mean(mcmc_store3$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store3,main='Store 3, beta_1')+geom_vline(xintercept=mean(mcmc_store3$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store3,main='Store 3, delta_0')+geom_vline(xintercept=mean(mcmc_store3$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store3,main='Store 3, delta_1')+geom_vline(xintercept=mean(mcmc_store3$delta1),color='red')

#store 1013 aka store 4
store4<-subset(beans,store==1013)
bin1_store4 = run_mcmc(y=store4$mvm,x=store4$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store4 = run_mcmc(y=store4$mvm,x=store4$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store4 = run_mcmc(y=store4$mvm,x=store4$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store4 <- run_mcmc(y=store4$mvm,x=store4$price,beta0=bin1_store4[[1]][1e3,], delta0=bin1_store4[[2]][1e3,],sigma.b=bin1_store4[[3]],sigma.d=bin1_store4[[4]],n.reps=1e4,tune=FALSE)
c2_store4 <- run_mcmc(y=store4$mvm,x=store4$price,beta0=bin2_store4[[1]][1e3,], delta0=bin2_store4[[2]][1e3,],sigma.b=bin2_store4[[3]],sigma.d=bin2_store4[[4]],n.reps=1e4,tune=FALSE)
c3_store4 <- run_mcmc(y=store4$mvm,x=store4$price,beta0=bin3_store4[[1]][1e3,], delta0=bin3_store4[[2]][1e3,],sigma.b=bin3_store4[[3]],sigma.d=bin3_store4[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store4<-NULL
mcmc_store4$delta0<-c(c1_store4[[1]][,1],c2_store4[[1]][,1],c3_store4[[1]][,1])
mcmc_store4$delta1<-c(c1_store4[[1]][,2],c2_store4[[1]][,2],c3_store4[[1]][,2])
mcmc_store4$beta0<-c(c1_store4[[2]][,1],c2_store4[[2]][,1],c3_store4[[2]][,1])
mcmc_store4$beta1<-c(c1_store4[[2]][,2],c2_store4[[2]][,2],c3_store4[[2]][,2])
mcmc_store4<-data.frame(mcmc_store4)
saveRDS(mcmc_store4,"store4Samples.RDS")

mcmc_store4<-readRDS("store4Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store4,main='Store 4, beta_0')+geom_vline(xintercept=mean(mcmc_store4$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store4,main='Store 4, beta_1')+geom_vline(xintercept=mean(mcmc_store4$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store4,main='Store 4, delta_0')+geom_vline(xintercept=mean(mcmc_store4$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store4,main='Store 4, delta_1')+geom_vline(xintercept=mean(mcmc_store4$delta1),color='red')

#store 1018 aka store 5
store5<-subset(beans,store==1018)
#store5<-store5[-83,]
bin1_store5 = run_mcmc(y=store5$mvm,x=store5$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store5 = run_mcmc(y=store5$mvm,x=store5$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store5 = run_mcmc(y=store5$mvm,x=store5$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store5 <- run_mcmc(y=store5$mvm,x=store5$price,beta0=bin1_store5[[1]][1e3,], delta0=bin1_store5[[2]][1e3,],sigma.b=bin1_store5[[3]],sigma.d=bin1_store5[[4]],n.reps=1e4,tune=FALSE)
c2_store5 <- run_mcmc(y=store5$mvm,x=store5$price,beta0=bin2_store5[[1]][1e3,], delta0=bin2_store5[[2]][1e3,],sigma.b=bin2_store5[[3]],sigma.d=bin2_store5[[4]],n.reps=1e4,tune=FALSE)
c3_store5 <- run_mcmc(y=store5$mvm,x=store5$price,beta0=bin3_store5[[1]][1e3,], delta0=bin3_store5[[2]][1e3,],sigma.b=bin3_store5[[3]],sigma.d=bin3_store5[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store5<-NULL
mcmc_store5$delta0<-c(c1_store5[[1]][,1],c2_store5[[1]][,1],c3_store5[[1]][,1])
mcmc_store5$delta1<-c(c1_store5[[1]][,2],c2_store5[[1]][,2],c3_store5[[1]][,2])
mcmc_store5$beta0<-c(c1_store5[[2]][,1],c2_store5[[2]][,1],c3_store5[[2]][,1])
mcmc_store5$beta1<-c(c1_store5[[2]][,2],c2_store5[[2]][,2],c3_store5[[2]][,2])
mcmc_store5<-data.frame(mcmc_store5)
saveRDS(mcmc_store5,"store5Samples.RDS")

mcmc_store5<-readRDS("store5Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store5,main='Store 5, beta_0')+geom_vline(xintercept=mean(mcmc_store5$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store5,main='Store 5, beta_1')+geom_vline(xintercept=mean(mcmc_store5$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store5,main='Store 5, delta_0')+geom_vline(xintercept=mean(mcmc_store5$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store5,main='Store 5, delta_1')+geom_vline(xintercept=mean(mcmc_store5$delta1),color='red')

#store 1019 aka store 6
store6<-subset(beans,store==1019)
bin1_store6 = run_mcmc(y=store6$mvm,x=store6$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store6 = run_mcmc(y=store6$mvm,x=store6$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store6 = run_mcmc(y=store6$mvm,x=store6$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store6 <- run_mcmc(y=store6$mvm,x=store6$price,beta0=bin1_store6[[1]][1e3,], delta0=bin1_store6[[2]][1e3,],sigma.b=bin1_store6[[3]],sigma.d=bin1_store6[[4]],n.reps=1e4,tune=FALSE)
c2_store6 <- run_mcmc(y=store6$mvm,x=store6$price,beta0=bin2_store6[[1]][1e3,], delta0=bin2_store6[[2]][1e3,],sigma.b=bin2_store6[[3]],sigma.d=bin2_store6[[4]],n.reps=1e4,tune=FALSE)
c3_store6 <- run_mcmc(y=store6$mvm,x=store6$price,beta0=bin3_store6[[1]][1e3,], delta0=bin3_store6[[2]][1e3,],sigma.b=bin3_store6[[3]],sigma.d=bin3_store6[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store6<-NULL
mcmc_store6$delta0<-c(c1_store6[[1]][,1],c2_store6[[1]][,1],c3_store6[[1]][,1])
mcmc_store6$delta1<-c(c1_store6[[1]][,2],c2_store6[[1]][,2],c3_store6[[1]][,2])
mcmc_store6$beta0<-c(c1_store6[[2]][,1],c2_store6[[2]][,1],c3_store6[[2]][,1])
mcmc_store6$beta1<-c(c1_store6[[2]][,2],c2_store6[[2]][,2],c3_store6[[2]][,2])
mcmc_store6<-data.frame(mcmc_store6)
saveRDS(mcmc_store6,"store6Samples.RDS")

mcmc_store6<-readRDS("store6Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store6,main='Store 6, beta_0')+geom_vline(xintercept=mean(mcmc_store6$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store6,main='Store 6, beta_1')+geom_vline(xintercept=mean(mcmc_store6$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store6,main='Store 6, delta_0')+geom_vline(xintercept=mean(mcmc_store6$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store6,main='Store 6, delta_1')+geom_vline(xintercept=mean(mcmc_store6$delta1),color='red')

#store 1022 aka store 7
store7<-subset(beans,store==1022)
bin1_store7 = run_mcmc(y=store7$mvm,x=store7$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store7 = run_mcmc(y=store7$mvm,x=store7$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store7 = run_mcmc(y=store7$mvm,x=store7$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store7 <- run_mcmc(y=store7$mvm,x=store7$price,beta0=bin1_store7[[1]][1e3,], delta0=bin1_store7[[2]][1e3,],sigma.b=bin1_store7[[3]],sigma.d=bin1_store7[[4]],n.reps=1e4,tune=FALSE)
c2_store7 <- run_mcmc(y=store7$mvm,x=store7$price,beta0=bin2_store7[[1]][1e3,], delta0=bin2_store7[[2]][1e3,],sigma.b=bin2_store7[[3]],sigma.d=bin2_store7[[4]],n.reps=1e4,tune=FALSE)
c3_store7 <- run_mcmc(y=store7$mvm,x=store7$price,beta0=bin3_store7[[1]][1e3,], delta0=bin3_store7[[2]][1e3,],sigma.b=bin3_store7[[3]],sigma.d=bin3_store7[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store7<-NULL
mcmc_store7$delta0<-c(c1_store7[[1]][,1],c2_store7[[1]][,1],c3_store7[[1]][,1])
mcmc_store7$delta1<-c(c1_store7[[1]][,2],c2_store7[[1]][,2],c3_store7[[1]][,2])
mcmc_store7$beta0<-c(c1_store7[[2]][,1],c2_store7[[2]][,1],c3_store7[[2]][,1])
mcmc_store7$beta1<-c(c1_store7[[2]][,2],c2_store7[[2]][,2],c3_store7[[2]][,2])
mcmc_store7<-data.frame(mcmc_store7)
saveRDS(mcmc_store7,"store7Samples.RDS")

mcmc_store7<-readRDS("store7Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store7,main='Store 7, beta_0')+geom_vline(xintercept=mean(mcmc_store7$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store7,main='Store 7, beta_1')+geom_vline(xintercept=mean(mcmc_store7$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store7,main='Store 7, delta_0')+geom_vline(xintercept=mean(mcmc_store7$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store7,main='Store 7, delta_1')+geom_vline(xintercept=mean(mcmc_store7$delta1),color='red')

#store 1026 aka store 8
store8<-subset(beans,store==1026)
bin1_store8 = run_mcmc(y=store8$mvm,x=store8$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store8 = run_mcmc(y=store8$mvm,x=store8$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store8 = run_mcmc(y=store8$mvm,x=store8$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store8 <- run_mcmc(y=store8$mvm,x=store8$price,beta0=bin1_store8[[1]][1e3,], delta0=bin1_store8[[2]][1e3,],sigma.b=bin1_store8[[3]],sigma.d=bin1_store8[[4]],n.reps=1e4,tune=FALSE)
c2_store8 <- run_mcmc(y=store8$mvm,x=store8$price,beta0=bin2_store8[[1]][1e3,], delta0=bin2_store8[[2]][1e3,],sigma.b=bin2_store8[[3]],sigma.d=bin2_store8[[4]],n.reps=1e4,tune=FALSE)
c3_store8 <- run_mcmc(y=store8$mvm,x=store8$price,beta0=bin3_store8[[1]][1e3,], delta0=bin3_store8[[2]][1e3,],sigma.b=bin3_store8[[3]],sigma.d=bin3_store8[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store8<-NULL
mcmc_store8$delta0<-c(c1_store8[[1]][,1],c2_store8[[1]][,1],c3_store8[[1]][,1])
mcmc_store8$delta1<-c(c1_store8[[1]][,2],c2_store8[[1]][,2],c3_store8[[1]][,2])
mcmc_store8$beta0<-c(c1_store8[[2]][,1],c2_store8[[2]][,1],c3_store8[[2]][,1])
mcmc_store8$beta1<-c(c1_store8[[2]][,2],c2_store8[[2]][,2],c3_store8[[2]][,2])
mcmc_store8<-data.frame(mcmc_store8)
saveRDS(mcmc_store8,"store8Samples.RDS")

mcmc_store8<-readRDS("store8Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store8,main='Store 8, beta_0')+geom_vline(xintercept=mean(mcmc_store8$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store8,main='Store 8, beta_1')+geom_vline(xintercept=mean(mcmc_store8$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store8,main='Store 8, delta_0')+geom_vline(xintercept=mean(mcmc_store8$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store8,main='Store 8, delta_1')+geom_vline(xintercept=mean(mcmc_store8$delta1),color='red')

#store 1027 aka store 9
store9<-subset(beans,store==1027)
bin1_store9 = run_mcmc(y=store9$mvm,x=store9$price,beta0=c(1,-1), delta0=c(1,-1),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)  #3 chains 
bin2_store9 = run_mcmc(y=store9$mvm,x=store9$price,beta0=c(2,-2), delta0=c(2,-2),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)
bin3_store9 = run_mcmc(y=store9$mvm,x=store9$price,beta0=c(0,0), delta0=c(0,0),sigma.b=2,sigma.d=2,Sigma=matrix(c(100,0,0,100),nrow=2),tune=T)

c1_store9 <- run_mcmc(y=store9$mvm,x=store9$price,beta0=bin1_store9[[1]][1e3,], delta0=bin1_store9[[2]][1e3,],sigma.b=bin1_store9[[3]],sigma.d=bin1_store9[[4]],n.reps=1e4,tune=FALSE)
c2_store9 <- run_mcmc(y=store9$mvm,x=store9$price,beta0=bin2_store9[[1]][1e3,], delta0=bin2_store9[[2]][1e3,],sigma.b=bin2_store9[[3]],sigma.d=bin2_store9[[4]],n.reps=1e4,tune=FALSE)
c3_store9 <- run_mcmc(y=store9$mvm,x=store9$price,beta0=bin3_store9[[1]][1e3,], delta0=bin3_store9[[2]][1e3,],sigma.b=bin3_store9[[3]],sigma.d=bin3_store9[[4]],n.reps=1e4,tune=FALSE)
#Change values to match ZIP section
mcmc_store9<-NULL
mcmc_store9$delta0<-c(c1_store9[[1]][,1],c2_store9[[1]][,1],c3_store9[[1]][,1])
mcmc_store9$delta1<-c(c1_store9[[1]][,2],c2_store9[[1]][,2],c3_store9[[1]][,2])
mcmc_store9$beta0<-c(c1_store9[[2]][,1],c2_store9[[2]][,1],c3_store9[[2]][,1])
mcmc_store9$beta1<-c(c1_store9[[2]][,2],c2_store9[[2]][,2],c3_store9[[2]][,2])
mcmc_store9<-data.frame(mcmc_store9)
saveRDS(mcmc_store9,"store9Samples.RDS")

mcmc_store9<-readRDS("store9Samples.RDS")
qplot(beta0,geom='histogram',binwidth=.1,data=mcmc_store9,main='Store 9, beta_0')+geom_vline(xintercept=mean(mcmc_store9$beta0),color='red')
qplot(beta1,geom='histogram',binwidth=.1,data=mcmc_store9,main='Store 9, beta_1')+geom_vline(xintercept=mean(mcmc_store9$beta1),color='red')
qplot(delta0,geom='histogram',binwidth=.01,data=mcmc_store9,main='Store 9, delta_0')+geom_vline(xintercept=mean(mcmc_store9$delta0),color='red')
qplot(delta1,geom='histogram',binwidth=.01,data=mcmc_store9,main='Store 9, delta_1')+geom_vline(xintercept=mean(mcmc_store9$delta1),color='red')
