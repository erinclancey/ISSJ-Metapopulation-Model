library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(latex2exp)

posterior <- read.csv("Joint Posterior Distribution.csv", header=TRUE)
rows <- posterior[sample(nrow(posterior), 375), ]

process <- function(eta, theta_oak, theta_pine, gamma, SDe, SDs, m){
  N_P = 216
  N_O = 1587
  delta = 2.03
  alpha = 0.1
  n_pine = 62
  n_oak = 453
  sex = 1
  SDa=1
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
       abs(cor(seq(1:50),Zbar_oak_F[451:500])* sd(Zbar_oak_F[451:500])/sd(seq(1:50))) < 0.005) {break}}
  
  delta_M <- Zbar_pine_M[r]-Zbar_oak_M[r]
  delta_F <- Zbar_pine_F[r]-Zbar_oak_F[r]
  return(data.frame(delta_M, delta_F))
  
}



library(foreach)
library(iterators)
library(parallel)
library(rngtools)
library(doParallel)
library(doRNG)
registerDoParallel()
registerDoRNG(2488820)
start_time <- Sys.time()
foreach(i=1:nrow(rows), .combine=rbind) %dopar%
  {library(vecsets)
    delta <- try(with(rows,process(eta[i], theta_oak[i], theta_pine[i], gamma[i], SDe[i], SDs[i], m[i])))
    if (class(delta) != "try-error") {delta} 
    else{
      delta$delta_M <- NA
      delta$delta_F <- NA
    }
  }-> deltazbar
end_time <- Sys.time()
end_time - start_time

#write.csv(deltazbar, file="M_F_deltazbar.csv")
deltazbar <- read.csv(file="M_F_deltazbar.csv", header=TRUE)

mean(deltazbar$delta_M)
mean(deltazbar$delta_F)
MOE_M <- sd(deltazbar$delta_M)*1.96
MOE_F <- sd(deltazbar$delta_F)*1.96



M_plot <- ggplot(deltazbar, aes(x=delta_M, color=delta_M, fill=delta_M))+
  scale_x_continuous(limits = c(-0.5,2))+
  geom_histogram(aes(y=..density..), color="grey30",fill="grey30",alpha=0.5, bins=30)+
  geom_density(alpha=0.2, color="grey30", fill="grey30", adjust=1.5)+
  labs(title="Male Divergence",x=TeX("$\\delta_{\\bar{z}_M}$"), y = "Density")+theme_minimal()+
  geom_vline(xintercept=0.794, color="blue", size=1, linetype=1)+
geom_vline(xintercept=mean(deltazbar$delta_M), color="black", size=1, linetype=2)+
  #theme(axis.title.x=element_text(size=18))+
  geom_segment(aes(x=mean(deltazbar$delta_M)-MOE_M,xend=mean(deltazbar$delta_M)+MOE_M),y=0,yend=0,color="black",size=2,lineend="round")
plot(M_plot)



F_plot <- ggplot(deltazbar, aes(x=delta_F, color=delta_F, fill=delta_F))+
  scale_x_continuous(limits = c(-0.5,2))+
  geom_histogram(aes(y=..density..), color="grey60",fill="grey60",alpha=0.5, bins=30)+
  geom_density(alpha=0.2, color="grey60", fill="grey60", adjust=1.5)+
  labs(title="Female Divergence",x=TeX("$\\delta_{\\bar{z}_F}$"), y = "Density")+theme_minimal()+
  geom_vline(xintercept=0.797, color="blue", size=1, linetype=1)+
  geom_vline(xintercept=mean(deltazbar$delta_F), color="black", size=1, linetype=2)+
  #theme(axis.title.x=element_text(size=18))+
geom_segment(aes(x=mean(deltazbar$delta_F)-MOE_F,xend=mean(deltazbar$delta_F)+MOE_F),y=0,yend=0,color="black",size=2,lineend="round")

plot(F_plot)


plot <- plot_grid(M_plot, F_plot, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot(plot)

