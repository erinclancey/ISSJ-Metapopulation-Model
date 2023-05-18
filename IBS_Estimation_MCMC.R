library(vecsets)
library(tidyverse)
library(MASS)
library(Boom)
library(coda)
library(bayestestR)
library(lhs)
library(MCMCglmm)


# Empirical Data
ISSJ_phenotypes <- read.csv("Island_Scrub_Jay_Phenotypes.csv", header=TRUE)
Z_P_M <- subset(ISSJ_phenotypes, Habitat=="Pine" & Sex=="M")$Ave_Nares_mm
Z_P_F <- subset(ISSJ_phenotypes, Habitat=="Pine" & Sex=="F")$Ave_Nares_mm
Z_O_M <- subset(ISSJ_phenotypes, Habitat=="Oak" & Sex=="M")$Ave_Nares_mm
Z_O_F <- subset(ISSJ_phenotypes, Habitat=="Oak" & Sex=="F")$Ave_Nares_mm

#Simulate the life cycle, make simulated data and calculate means and SD's
process <- function(eta, theta_oak, theta_pine, gamma, SDe, SDs, m){
  N_P = 216
  N_O = 1587
  delta = 2.03
  alpha = 0.1
  n_pine = 62
  n_oak = 453
  sex = 1
  SDa=1.0
  H_P = 0.12
  H_O = 0.88
  start = theta_oak
  
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
  
  R=1000
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
    
    if(r>100 & abs(cor(seq(1:50),Zbar_pine_M[51:100])* sd(Zbar_pine_M[51:100])/sd(seq(1:50))) < 0.005 &
        abs(cor(seq(1:50),Zbar_pine_F[51:100])* sd(Zbar_pine_F[51:100])/sd(seq(1:50))) < 0.005 &
        abs(cor(seq(1:50),Zbar_oak_M[51:100])* sd(Zbar_oak_M[51:100])/sd(seq(1:50))) < 0.005 &
        abs(cor(seq(1:50),Zbar_oak_F[51:100])* sd(Zbar_oak_F[51:100])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>150 & abs(cor(seq(1:50),Zbar_pine_M[101:150])* sd(Zbar_pine_M[101:150])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[101:150])* sd(Zbar_pine_F[101:150])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[101:150])* sd(Zbar_oak_M[101:150])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[101:150])* sd(Zbar_oak_F[101:150])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>200 & abs(cor(seq(1:50),Zbar_pine_M[151:200])* sd(Zbar_pine_M[151:200])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[151:200])* sd(Zbar_pine_F[151:200])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[151:200])* sd(Zbar_oak_M[151:200])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[151:200])* sd(Zbar_oak_F[151:200])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>250 & abs(cor(seq(1:50),Zbar_pine_M[201:250])* sd(Zbar_pine_M[201:250])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[201:250])* sd(Zbar_pine_F[201:250])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[201:250])* sd(Zbar_oak_M[201:250])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[201:250])* sd(Zbar_oak_F[201:250])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>300 & abs(cor(seq(1:50),Zbar_pine_M[251:300])* sd(Zbar_pine_M[251:300])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[251:300])* sd(Zbar_pine_F[251:300])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[251:300])* sd(Zbar_oak_M[251:300])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[251:300])* sd(Zbar_oak_F[251:300])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>350 & abs(cor(seq(1:50),Zbar_pine_M[301:350])* sd(Zbar_pine_M[301:350])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[301:350])* sd(Zbar_pine_F[301:350])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[301:350])* sd(Zbar_oak_M[301:350])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[301:350])* sd(Zbar_oak_F[301:350])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>400 & abs(cor(seq(1:50),Zbar_pine_M[351:400])* sd(Zbar_pine_M[351:400])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[351:400])* sd(Zbar_pine_F[351:400])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[351:400])* sd(Zbar_oak_M[351:400])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[351:400])* sd(Zbar_oak_F[351:400])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>450 & abs(cor(seq(1:50),Zbar_pine_M[401:450])* sd(Zbar_pine_M[401:450])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[401:450])* sd(Zbar_pine_F[401:450])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[401:450])* sd(Zbar_oak_M[401:450])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[401:450])* sd(Zbar_oak_F[401:450])/sd(seq(1:50))) < 0.005) {break}
    
    if(r>500 & abs(cor(seq(1:50),Zbar_pine_M[451:500])* sd(Zbar_pine_M[451:500])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_pine_F[451:500])* sd(Zbar_pine_F[451:500])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_M[451:500])* sd(Zbar_oak_M[451:500])/sd(seq(1:50))) < 0.005 &
       abs(cor(seq(1:50),Zbar_oak_F[451:500])* sd(Zbar_oak_F[451:500])/sd(seq(1:50))) < 0.005) {break}
    
  }

  return(list(Zbar_pine_M[r], SDp_pine_M[r],
              Zbar_pine_F[r],SDp_pine_F[r],
              Zbar_oak_M[r],SDp_oak_M[r],
              Zbar_oak_F[r], SDp_oak_F[r]
              ))
  
}

# Function to calculate the value of the log likelihood
calcLogLik <- function(eta, theta_oak, theta_pine, gamma,SDe, SDs,m, summary, Z_P_M, Z_P_F, Z_O_M, Z_O_F){
  if (eta<0|theta_oak<0|theta_pine<0|gamma<0|SDe<0|SDs<0|m<0){
    loglik <- log(0)
  } else{
  loglik <- sum(dnorm(Z_P_M, summary[[1]], summary[[2]], log = TRUE)) + 
    sum(dnorm(Z_P_F, summary[[3]], summary[[4]], log = TRUE)) +
    sum(dnorm(Z_O_M, summary[[5]], summary[[6]], log = TRUE)) +
    sum(dnorm(Z_O_F, summary[[7]], summary[[8]], log = TRUE) + dgamma(m, shape=9.4, scale=0.01)) 
  }
  return(loglik)
}

Metro_Hast <- function (eta0, theta_oak0, theta_pine0, gamma0, SDe0, SDs0, m0, iter, Z_P_M, Z_P_F, Z_O_M, Z_O_F) {

eta=rep(NA, iter) 
eta[1]=eta0
theta_oak=rep(NA, iter)
theta_oak[1]=theta_oak0
theta_pine=rep(NA, iter)
theta_pine[1]=theta_pine0
gamma=rep(iter)
gamma[1]=gamma0
SDe=rep(NA, iter)
SDe[1]=SDe0
SDs=rep(NA, iter)
SDs[1]=SDs0
m=rep(NA, iter) 
m[1]=m0
lik=rep(NA, iter)
lik[1]=calcLogLik(eta0, theta_oak0, theta_pine0, gamma0, SDe0, SDs0, m0, process(eta0, theta_oak0, theta_pine0, gamma0, SDe0, SDs0, m0), Z_P_M, Z_P_F, Z_O_M, Z_O_F)
accept=rep(NA, iter)
accept[1]=0

prop.sd <- 0.193*matrix(c(0.000206881,-0.000479629,-0.000546419,-2.28E-05,-0.00043286,-0.000323722,-1.43E-06,
                         -0.000479629,0.036124693,-0.02545953,0.000522664,-0.001943777,-0.001929837,-7.78E-05,
                         -0.000546419,-0.02545953,0.087116475,-0.001559709,0.004462354,0.004449212,-0.000214598,
                         -2.28E-05,0.000522664,-0.001559709,0.000112756,-0.000161981,-7.21E-05,9.83E-06,
                         -0.00043286,-0.001943777,0.004462354,-0.000161981,0.008758784,-0.00098377,-3.35E-05,
                         -0.000323722,-0.001929837,0.004449212,-7.21E-05,-0.00098377,0.00577347,-1.33E-05,
                         -1.43E-06,-7.78E-05,-0.000214598,9.83E-06,-3.35E-05,-1.33E-05,9.11E-06), ncol = 7)


  for (i in 2:iter) {
    
    summary.old <- process(eta[i-1], theta_oak[i-1], theta_pine[i-1], gamma[i-1], SDe[i-1], SDs[i-1], m[i-1])
    
    mu <- c(eta[i-1], theta_oak[i-1], theta_pine[i-1], gamma[i-1], SDe[i-1], SDs[i-1], m[i-1])
    
    theta.new <- mvrnorm(1,mu, prop.sd)
    
    summary.new <- try(process(theta.new[1],theta.new[2],theta.new[3],theta.new[4],theta.new[5],theta.new[6],theta.new[7]))
    if (class(summary.new) != "try-error") {
      loglik.new <- calcLogLik(theta.new[1],theta.new[2],theta.new[3],theta.new[4],theta.new[5],theta.new[6],theta.new[7], summary.new, Z_P_M, Z_P_F, Z_O_M, Z_O_F)
      loglik.old <- calcLogLik(eta[i-1], theta_oak[i-1], theta_pine[i-1], gamma[i-1],SDe[i-1], SDs[i-1],m[i-1], summary.old, Z_P_M, Z_P_F, Z_O_M, Z_O_F)
      log.ratio <- loglik.new - loglik.old
      log.accept.prob <- min(log.ratio, 0)
    } 
    else {log.accept.prob <- -Inf}
    
    if (log (runif (1)) < log.accept.prob & 
                                  theta.new[1]<0.3 & 
                                  theta.new[2]>20 & theta.new[2]<30 &
                                  theta.new[3]>20 & theta.new[3]<30 &
                                  theta.new[4]<0.3 & 
                                  theta.new[5]>0.5 & theta.new[5]<1.5 & 
                                  theta.new[6]>0.5 & theta.new[6]<1.5) {
      eta[i]=theta.new[1]
      theta_oak[i]=theta.new[2]
      theta_pine[i]=theta.new[3]
      gamma[i]=theta.new[4]
      SDe[i]=theta.new[5]
      SDs[i]=theta.new[6]
      m[i]=theta.new[7]
      lik[i]=loglik.new
      accept[i]=1
    } else {
      eta[i]=eta[i-1]
      theta_oak[i]=theta_oak[i-1]
      theta_pine[i]=theta_pine[i-1]
      gamma[i]=gamma[i-1]
      SDe[i]=SDe[i-1]
      SDs[i]=SDs[i-1]
      m[i]=m[i-1]
      lik[i]=loglik.old
      accept[i]=0
    } 
    print(i)
  }
  return(list(eta, theta_oak, theta_pine, gamma, SDe, SDs, m, lik, accept))
}

#starting generated values from pilot runs
eta=0.1716631
theta_oak=23.6521
theta_pine=25.25092
gamma=0.006560123
SDe=0.7925147
SDs=0.6594796
m=0.084

results <- Metro_Hast(eta, theta_oak, theta_pine, gamma, SDe, SDs, m, iter=30000, Z_P_M, Z_P_F, Z_O_M, Z_O_F)
results.post <- do.call(cbind, results)

