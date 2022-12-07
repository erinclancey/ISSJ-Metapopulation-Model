library(vecsets)
library(tidyverse)
library(MASS)
library(Boom)
library(coda)
library(bayestestR)
library(lhs)
library(MCMCglmm)

#Simulate the life cycle, make simulated data and calculate means and SD's
process <- function(eta, theta_oak, theta_pine, gamma, SDa, SDe, SDs, m){
  N_P = 216
  N_O = 1587
  delta = 2.03
  alpha = 0.1
  n_pine = 62
  n_oak = 453
  sex = 1
  H_P = 0.12
  H_O = 0.88
  start = mean(theta_oak, theta_pine)
  
  BV_pine_f <- rnorm(N_P/2, mean = start, sd = SDa)
  BV_pine_m <- rnorm(N_P/2, mean = start, sd = SDa)
  BV_oak_f <- rnorm(N_O/2, mean = start, sd = SDa)
  BV_oak_m <- rnorm(N_O/2, mean = start, sd = SDa)
  
  Z_pine_f <- as.matrix(cbind(BV_pine_f, rnorm(N_P/2, mean = BV_pine_f-sex, sd = SDe)))
  colnames(Z_pine_f) = NULL
  Z_pine_m <- as.matrix(cbind(BV_pine_m, rnorm(N_P/2, mean = BV_pine_m+sex, sd = SDe)))
  colnames(Z_pine_m) = NULL
  Z_oak_f <- as.matrix(cbind(BV_oak_f, rnorm(N_O/2, mean = BV_oak_f-sex, sd = SDe)))
  colnames(Z_oak_f) = NULL
  Z_oak_m <- as.matrix(cbind(BV_oak_m, rnorm(N_O/2, mean = BV_oak_m+sex, sd = SDe)))
  colnames(Z_oak_m) = NULL
  
  R=100
  Zbar_pine_M <- rep(NA, R)
  Zbar_pine_F <- rep(NA, R)
  Zbar_oak_M <- rep(NA, R)
  Zbar_oak_F <- rep(NA, R)
  SDp_pine_M <-rep(NA, R)
  SDp_pine_F <-rep(NA, R)
  SDp_oak_M <-rep(NA, R)
  SDp_oak_F <-rep(NA, R)
  
  for (r in 1:R){
    Zbar_pine_M[r] <- mean(Z_pine_m[,2])
    Zbar_pine_F[r] <- mean(Z_pine_f[,2])
    Zbar_oak_M[r] <- mean(Z_oak_m[,2])
    Zbar_oak_F[r] <- mean(Z_oak_f[,2])
    SDp_pine_M[r]<-sd(Z_pine_m[,2])
    SDp_pine_F[r]<-sd(Z_pine_f[,2])
    SDp_oak_M[r]<-sd(Z_oak_m[,2])
    SDp_oak_F[r]<-sd(Z_oak_f[,2])
    
    #### Mate and make offspring in Pine
    Pairs_pine <- data.frame(FemBV = rep(NA, n_pine), FemZ = rep(NA,n_pine) , MalBV = rep(NA, n_pine), MalZ = rep(NA,n_pine))
    for (i in 1:n_pine){ # pairs breed in pine
      repeat{
        Pfem <- Z_pine_f[sample(nrow(Z_pine_f),size=1,replace=TRUE),, drop=FALSE]
        Pmal <- Z_pine_m[sample(nrow(Z_pine_m),size=1,replace=TRUE),, drop=FALSE]
        mate <- exp(-alpha*(Pmal[,2] - Pfem[,2] - delta)^2)
        if (mate > runif(1)){
          pair <- cbind(Pfem, Pmal)
          Pairs_pine[i,] <- pair
          break
        }
      }
    }
    
    BV <- cbind(Pairs_pine[,1], Pairs_pine[,3])
    BVmean <- apply(BV, 1, mean)
    BV_off_pine <- BVmean + rnorm(length(BVmean), 0, SDs)
    Z_off_pine <- BV_off_pine + rnorm(length(BV_off_pine), 0, SDe)
    sex_off_pine <- sample(c("M","F"), length(BV_off_pine), replace=TRUE)
    
    Zoff_pineMatrix <- as.matrix(cbind(BV_off_pine, Z_off_pine, sex_off_pine))
    colnames(Zoff_pineMatrix) = NULL
    
    Zoff_pine_f <- subset(Zoff_pineMatrix, Zoff_pineMatrix[,3]=="F")
    Zoff_pine_m <- subset(Zoff_pineMatrix, Zoff_pineMatrix[,3]=="M")
    
    Zoff_pine_f <- cbind(as.numeric(Zoff_pine_f[,1]), as.numeric(Zoff_pine_f[,2]) - sex)
    Zoff_pine_m <- cbind(as.numeric(Zoff_pine_m[,1]), as.numeric(Zoff_pine_m[,2]) + sex)
    
    Z_pine_m <- rbind(Z_pine_m, Zoff_pine_m)
    Z_pine_f <- rbind(Z_pine_f, Zoff_pine_f)
    
    #### Mate and make offspring in Oak
    
    Pairs_oak <- data.frame(FemBV = rep(NA, n_oak), FemZ = rep(NA,n_oak) , MalBV = rep(NA, n_oak), MalZ = rep(NA,n_oak))
    for (i in 1:n_oak){ # K pairs breed in oak
      repeat{
        Pfem <- Z_oak_f[sample(nrow(Z_oak_f),size=1,replace=TRUE),, drop=FALSE]
        Pmal <- Z_oak_m[sample(nrow(Z_oak_m),size=1,replace=TRUE),, drop=FALSE]
        mate <- exp(-alpha*(Pmal[,2] - Pfem[,2] - delta)^2)
        if (mate > runif(1)){
          pair <- cbind(Pfem, Pmal)
          Pairs_oak[i,] <- pair
          break
        }
      }
    }
    
    BV <- cbind(Pairs_oak[,1], Pairs_oak[,3])
    BVmean <- apply(BV, 1, mean)
    BV_off_oak <- BVmean + rnorm(length(BVmean), 0, SDs)
    Z_off_oak <- BV_off_oak + rnorm(length(BV_off_oak), 0, SDe)
    sex_off_oak <- sample(c("M","F"), length(BV_off_oak), replace=TRUE)
    
    Zoff_oakMatrix <- as.matrix(cbind(BV_off_oak, Z_off_oak, sex_off_oak))
    colnames(Zoff_oakMatrix) = NULL
    
    Zoff_oak_f <- subset(Zoff_oakMatrix, Zoff_oakMatrix[,3]=="F")
    Zoff_oak_m <- subset(Zoff_oakMatrix, Zoff_oakMatrix[,3]=="M")
    
    Zoff_oak_f <- cbind(as.numeric(Zoff_oak_f[,1]), as.numeric(Zoff_oak_f[,2]) - sex)
    Zoff_oak_m <- cbind(as.numeric(Zoff_oak_m[,1]), as.numeric(Zoff_oak_m[,2]) + sex)
    
    Z_oak_m <- rbind(Z_oak_m, Zoff_oak_m)
    Z_oak_f <- rbind(Z_oak_f, Zoff_oak_f)
    
    ##### Migration
    M_to_pine_f <- Z_oak_f[2*m*H_P*exp(-eta*(Z_oak_f[,2] - theta_pine)^2) / ( exp(-eta*(Z_oak_f[,2] - theta_pine)^2) + exp(-eta*(Z_oak_f[,2] - theta_oak)^2)) > runif(length(Z_oak_f[,2]), 0, 1),, drop=FALSE]
    M_to_pine_m <- Z_oak_m[2*m*H_P*exp(-eta*(Z_oak_m[,2] - theta_pine)^2) / ( exp(-eta*(Z_oak_m[,2] - theta_pine)^2) + exp(-eta*(Z_oak_m[,2] - theta_oak)^2)) > runif(length(Z_oak_m[,2]), 0, 1),, drop=FALSE]
    
    M_to_oak_f <- Z_pine_f[2*m*H_O*exp(-eta*(Z_pine_f[,2] - theta_oak)^2) / ( exp(-eta*(Z_pine_f[,2] - theta_oak)^2) + exp(-eta*(Z_pine_f[,2] - theta_pine)^2)) > runif(length(Z_pine_f[,2]), 0, 1),, drop=FALSE]
    M_to_oak_m <- Z_pine_m[2*m*H_O *exp(-eta*(Z_pine_m[,2] - theta_oak)^2) / ( exp(-eta*(Z_pine_m[,2] - theta_oak)^2) + exp(-eta*(Z_pine_m[,2] - theta_pine)^2)) > runif(length(Z_pine_m[,2]), 0, 1),, drop=FALSE]
    
    # Remove migrants from each population
    Z_pine_f <- cbind(vsetdiff(Z_pine_f[,1], M_to_oak_f[,1]), vsetdiff(Z_pine_f[,2], M_to_oak_f[,2]))
    Z_pine_m <- cbind(vsetdiff(Z_pine_m[,1], M_to_oak_m[,1]), vsetdiff(Z_pine_m[,2], M_to_oak_m[,2]))
    
    Z_oak_f <- cbind(vsetdiff(Z_oak_f[,1], M_to_pine_f[,1]), vsetdiff(Z_oak_f[,2], M_to_pine_f[,2]))
    Z_oak_m <- cbind(vsetdiff(Z_oak_m[,1], M_to_pine_m[,1]), vsetdiff(Z_oak_m[,2], M_to_pine_m[,2]))
    
    ## Pine after migration
    Z_pine_f <- rbind(Z_pine_f, M_to_pine_f)
    Z_pine_m <- rbind(Z_pine_m, M_to_pine_m)
    
    ## Oak after migration
    Z_oak_f <- rbind(Z_oak_f, M_to_oak_f)
    Z_oak_m <- rbind(Z_oak_m, M_to_oak_m)
    
    ##### Selection in Pine
    Z_pine_f <- Z_pine_f[exp(-gamma*(Z_pine_f[,2] - theta_pine)^2) > runif(length(Z_pine_f[,2]), 0, 1),]
    Z_pine_m <- Z_pine_m[exp(-gamma*(Z_pine_m[,2] - theta_pine)^2) > runif(length(Z_pine_m[,2]), 0, 1),]
    
    #### density dependent mortality
    if (length(Z_pine_f[,2]) > N_P/2){
      Z_pine_f <- Z_pine_f[sample(nrow(Z_pine_f), N_P/2, replace=FALSE),]
    }
    
    if (length(Z_pine_m[,2]) > N_P/2){
      Z_pine_m <- Z_pine_m[sample(nrow(Z_pine_m), N_P/2, replace=FALSE),]
    }
    
    ###### Selection in Oak
    Z_oak_f <- Z_oak_f[exp(-gamma*(Z_oak_f[,2] - theta_oak)^2) > runif(length(Z_oak_f[,2]), 0, 1),]
    Z_oak_m <- Z_oak_m[exp(-gamma*(Z_oak_m[,2] - theta_oak)^2) > runif(length(Z_oak_m[,2]), 0, 1),]
    
    ## density dependent mortality
    if (length(Z_oak_f[,2]) > N_O/2){
      Z_oak_f <- Z_oak_f[sample(nrow(Z_oak_f), N_O/2, replace=FALSE),]
    }
    
    if (length(Z_oak_m[,2]) > N_O/2){
      Z_oak_m <- Z_oak_m[sample(nrow(Z_oak_m), N_O/2, replace=FALSE),]
    }
  }
  
  return(list(Zbar_pine_M[R], SDp_pine_M[R],
              Zbar_pine_F[R],SDp_pine_F[R],
              Zbar_oak_M[R],SDp_oak_M[R],
              Zbar_oak_F[R], SDp_oak_F[R],
              Z_pine_m[,2], Z_pine_f[,2], Z_oak_m[,2], Z_oak_f[,2]
  ))
  
}

# Function to calculate the value of the log likelihood
calcLogLik <- function(eta, theta_oak, theta_pine, gamma,SDa, SDe, SDs,m, summary, Z_P_M, Z_P_F, Z_O_M, Z_O_F){
  if (eta<0|theta_oak<0|theta_pine<0|gamma<0|SDa<0|SDe<0|SDs<0|m<0){
    loglik <- log(0)
  } else{
    loglik <- sum(dnorm(Z_P_M, summary[[1]], summary[[2]], log = TRUE)) + 
      sum(dnorm(Z_P_F, summary[[3]], summary[[4]], log = TRUE)) +
      sum(dnorm(Z_O_M, summary[[5]], summary[[6]], log = TRUE)) +
      sum(dnorm(Z_O_F, summary[[7]], summary[[8]], log = TRUE)) + dgamma(m, shape=9.4, scale=0.01, log=TRUE) 
  }
  return(loglik)
}


Metro_Hast <- function (eta0, theta_oak0, theta_pine0, gamma0, SDa0, SDe0, SDs0, m0, iter, Z_P_M, Z_P_F, Z_O_M, Z_O_F) {
  
  eta=rep(NA, iter) 
  eta[1]=eta0
  theta_oak=rep(NA, iter)
  theta_oak[1]=theta_oak0
  theta_pine=rep(NA, iter)
  theta_pine[1]=theta_pine0
  gamma=rep(NA, iter)
  gamma[1]=gamma0
  SDa=rep(NA, iter)
  SDa[1]=SDa0
  SDe=rep(NA, iter)
  SDe[1]=SDe0
  SDs=rep(NA, iter)
  SDs[1]=SDs0
  m=rep(NA, iter) 
  m[1]=m0
  lik=rep(NA, iter)
  lik[1]=calcLogLik(eta0, theta_oak0, theta_pine0, gamma0,SDa0, SDe0, SDs0,m0, process(eta0, theta_oak0, theta_pine0, gamma0, SDa0, SDe0, SDs0, m0), Z_P_M, Z_P_F, Z_O_M, Z_O_F)
  
  prop.sd <- c(0.05, 0.5, 0.5, 0.01, 0.1, 0.1, 0.1, 0.05)
  
  for (i in 2:iter) {
    
    summary.old <- process(eta[i-1], theta_oak[i-1], theta_pine[i-1], gamma[i-1], SDa[i-1], SDe[i-1], SDs[i-1], m[i-1])
    
    eta.new <- rnorm(1, eta[i-1], sd=prop.sd[1])
    theta_oak.new <- rnorm(1, theta_oak[i-1], sd=prop.sd[2])
    theta_pine.new <- rnorm(1, theta_pine[i-1], sd=prop.sd[3])
    gamma.new <- rnorm(1, gamma[i-1], sd=prop.sd[4])
    SDa.new <- rnorm(1, SDa[i-1], sd=prop.sd[5])
    SDe.new <- rnorm(1, SDe[i-1], sd=prop.sd[6])
    SDs.new <- rnorm(1, SDs[i-1], sd=prop.sd[7])
    m.new <- rnorm(1, m[i-1], sd=prop.sd[7])
    
    theta.new <- c(eta.new, theta_oak.new, theta_pine.new, gamma.new, SDa.new, SDe.new, SDs.new, m.new)
    
    summary.new <- try(process(eta.new, theta_oak.new, theta_pine.new, gamma.new, SDa.new, SDe.new, SDs.new, m.new))
    if (class(summary.new) != "try-error") {
      loglik.new <- calcLogLik(eta.new, theta_oak.new, theta_pine.new, gamma.new,SDa.new,SDe.new,SDs.new,m.new, summary.new, Z_P_M, Z_P_F, Z_O_M, Z_O_F)
      loglik.old <- calcLogLik(eta[i-1], theta_oak[i-1], theta_pine[i-1], gamma[i-1],SDa[i-1], SDe[i-1], SDs[i-1], m[i-1], summary.old, Z_P_M, Z_P_F, Z_O_M, Z_O_F)
      log.ratio <- loglik.new - loglik.old
      log.accept.prob <- min(log.ratio, 0)
    } 
    else {log.accept.prob <- -Inf}
    
    if (log (runif (1)) < log.accept.prob && 
        eta.new<0.3 && 
        theta_oak.new>21.57 && theta_oak.new<23.79 &&
        theta_pine.new>24.67 && theta_pine.new<26.95 &&
        gamma.new<0.3 && 
        SDa.new>0.74 && SDa.new<2.0 && 
        SDe.new>0.5 && SDe.new<1.0 && 
        SDs.new>0.5 && SDs.new<1.0 &&
        m.new>0.01 && m.new<0.3) {
      eta[i]=eta.new
      theta_oak[i]=theta_oak.new
      theta_pine[i]=theta_pine.new
      gamma[i]=gamma.new
      SDa[i]=SDa.new
      SDe[i]=SDe.new
      SDs[i]=SDs.new
      m[i]=m.new
      lik[i]=loglik.new
    } else {
      eta[i]=eta[i-1]
      theta_oak[i]=theta_oak[i-1]
      theta_pine[i]=theta_pine[i-1]
      gamma[i]=gamma[i-1]
      SDa[i]=SDa[i-1]
      SDe[i]=SDe[i-1]
      SDs[i]=SDs[i-1]
      m[i]=m[i-1]
      lik[i]=loglik.old
    } 
    print(i)
  }
  return(list(eta, theta_oak, theta_pine, gamma, SDa, SDe, SDs, m, lik))
}

rows=20
names=c("eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m")
X <- randomLHS(rows, length(names))
X[,1] <- qunif(X[,1], 0, 0.3)
X[,2] <- qunif(X[,2], 21.57, 23.79)
X[,3] <- qunif(X[,3], 24.67, 26.95)
X[,4] <- qunif(X[,4], 0, 0.3)
X[,5] <- qunif(X[,5], 0.74, 2)
X[,6] <- qunif(X[,6], 0.5, 1)
X[,7] <- qunif(X[,7], 0.5, 1)
X[,8] <- qunif(X[,8], 0.01, 0.3)
params <- data.frame(X)
colnames(params)<- names

estimates <- data.frame(cbind(params, eta.mode=rep(NA, rows), theta_oak.mode=rep(NA, rows), theta_pine.mode=rep(NA, rows), 
                              gamma.mode=rep(NA, rows), SDa.mode=rep(NA, rows), SDe.mode=rep(NA, rows), SDs.mode=rep(NA, rows), m.mode=rep(NA, rows) ))

first=100
thin=50
iter=5000

for (i in 1:nrow(estimates)) {
  print(paste(round(i/nrow(estimates) * 100, 1), "%", sep = ""))
  
  sim_data <- try(with(estimates, process(eta[i], theta_oak[i], theta_pine[i], gamma[i], SDa[i], SDe[i], SDs[i], m[i])))
  
  if (class(sim_data) != "try-error") {
    Z_P_M <- sim_data[[9]]
    Z_P_F <- sim_data[[10]]
    Z_O_M <- sim_data[[11]]
    Z_O_F <- sim_data[[12]]
    results <- try(with(estimates, Metro_Hast(eta[i], theta_oak[i], theta_pine[i], gamma[i], SDa[i], SDe[i], SDs[i],m[i], 
                                         iter=iter, Z_P_M, Z_P_F, Z_O_M, Z_O_F)))
    
    if (class(results) != "try-error") {
    results.eta <- as.vector(window(mcmc(results[[1]]), start=first, end=iter, thin=thin))
    results.theta_oak <- as.vector(window(mcmc(results[[2]]), start=first, end=iter, thin=thin))
    results.theta_pine <- as.vector(window(mcmc(results[[3]]), start=first, end=iter, thin=thin))
    results.gamma <- as.vector(window(mcmc(results[[4]]), start=first, end=iter, thin=thin))
    results.SDa <- as.vector(window(mcmc(results[[5]]), start=first, end=iter, thin=thin))
    results.SDe <- as.vector(window(mcmc(results[[6]]), start=first, end=iter, thin=thin))
    results.SDs <- as.vector(window(mcmc(results[[7]]), start=first, end=iter, thin=thin))
    results.m <- as.vector(window(mcmc(results[[8]]), start=first, end=iter, thin=thin))
    
    estimates$eta.mode[i] <- posterior.mode(mcmc(results.eta),adjust=1)
    estimates$theta_oak.mode[i] <- posterior.mode(mcmc(results.theta_oak),adjust=1)
    estimates$theta_pine.mode[i] <- posterior.mode(mcmc(results.theta_pine),adjust=1)
    estimates$gamma.mode[i] <- posterior.mode(mcmc(results.gamma),adjust=1)
    estimates$SDa.mode[i] <- posterior.mode(mcmc(results.SDa),adjust=1)
    estimates$SDe.mode[i] <- posterior.mode(mcmc(results.SDe),adjust=1)
    estimates$SDs.mode[i] <- posterior.mode(mcmc(results.SDs),adjust=1)
    estimates$m.mode[i] <- posterior.mode(mcmc(results.m),adjust=1)
    }else {
      estimates$eta.mode[i] <- NA
      estimates$theta_oak.mode[i] <- NA
      estimates$theta_pine.mode[i] <- NA
      estimates$gamma.mode[i] <- NA
      estimates$SDa.mode[i] <- NA
      estimates$SDe.mode[i] <- NA
      estimates$SDs.mode[i] <- NA
      estimates$m.mode[i] <- NA
    }
    
  } else {
    estimates$eta.mode[i] <- NA
    estimates$theta_oak.mode[i] <- NA
    estimates$theta_pine.mode[i] <- NA
    estimates$gamma.mode[i] <- NA
    estimates$SDa.mode[i] <- NA
    estimates$SDe.mode[i] <- NA
    estimates$SDs.mode[i] <- NA
    estimates$m.mode[i] <- NA
  }
}

write.csv(estimates, "estimates_1.csv")



