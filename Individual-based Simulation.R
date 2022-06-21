library(vecsets)

N_P = 216
N_O = 1587
a = 0.1 # strength of mate preference
d = 2.06 # match offset
gamma = 1 # max. prob. of mating
n_pine = 62 # number of nesting sites in pine
n_oak = 453 # number of nesting sites in oak
SDa = SDa # additive genetic standard deviation
SDe = SDe # environmental standard deviation
SDs = SDs # segregational standard deviation
sex = 1 # fixed sex effect
eta = eta # strength of habitat preference
m = 0.25 # max prob of migration
H_P = 0.12 # proportion of pine habitat
H_O = 0.88 # proportion of oak habitat
theta_pine = theta_pine # ecological optimum in pine
theta_oak = theta_oak # ecological optimum in oak
S = S # strength of selection

############# Initialize population
start=mean(c(theta_oak, theta_pine))
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

R=300 # number of generations
Zbar_pine <- rep(NA, R)
Zbar_oak <- rep(NA, R)
deltaZbar <- rep(NA, R)
Vp_pine <-rep(NA, R)
Vp_pine_f <-rep(NA, R)
Vp_pine_m <-rep(NA, R)
Vp_oak <-rep(NA, R)
Vp_oak_f <-rep(NA, R)
Vp_oak_m <-rep(NA, R)
Va_pine <-rep(NA, R)
Va_oak <-rep(NA, R)
h2 <- rep(NA, R)
N_pine <- rep(NA, R)
N_oak<- rep(NA, R)
BVbar_pine <- rep(NA, R)
BVbar_oak <- rep(NA, R)
Vp_POP<- rep(NA, R)
Va_POP<- rep(NA, R)

############ Run Simulations
for (r in 1:R){
   
   N_pine[r] <- length(c(Z_pine_f[,1], Z_pine_m[,1]))
   N_oak[r] <- length(c(Z_oak_f[,1], Z_oak_m[,1]))

   Zbar_pine[r] <- mean(c(Z_pine_f[,2], Z_pine_m[,2]))
   Zbar_oak[r] <- mean(c(Z_oak_f[,2], Z_oak_m[,2]))
   
   BVbar_pine[r] <- mean(c(Z_pine_f[,1], Z_pine_m[,1]))
   BVbar_oak[r] <- mean(c(Z_oak_f[,1], Z_oak_m[,1]))

   Va_pine[r]<-var(c(Z_pine_f[,1], Z_pine_m[,1]))
   Va_oak[r]<-var(c(Z_oak_f[,1], Z_oak_m[,1]))
   Va_POP[r] <- var(c(Z_pine_f[,1], Z_pine_m[,1], Z_oak_f[,1], Z_oak_m[,1]))
   
   Vp_pine[r]<-var(c(Z_pine_f[,2], Z_pine_m[,2]))
   Vp_pine_f[r]<-var(Z_pine_f[,2])
   Vp_pine_m[r]<-var(Z_pine_m[,2])
   Vp_oak[r]<-var(c(Z_oak_f[,2], Z_oak_m[,2]))
   Vp_oak_f[r]<-var(Z_oak_f[,2])
   Vp_oak_m[r]<-var(Z_oak_m[,2])
   Vp_POP[r] <- var(c(Z_pine_f[,2], Z_pine_m[,2], Z_oak_f[,2], Z_oak_m[,2]))
   
   h2[r] <- var(c(Z_pine_f[,1], Z_pine_m[,1], Z_oak_f[,1], Z_oak_m[,1])) / var(c(Z_pine_f[,2], Z_pine_m[,2], Z_oak_f[,2], Z_oak_m[,2]))



#### Mate and make offspring in Pine
   Pairs_pine <- data.frame(FemBV = rep(NA, n_pine), FemZ = rep(NA,n_pine) , MalBV = rep(NA, n_pine), MalZ = rep(NA,n_pine))
   for (i in 1:n_pine){ # pairs breed in pine
      repeat{
         Pfem <- Z_pine_f[sample(nrow(Z_pine_f),size=1,replace=TRUE),, drop=FALSE]
         Pmal <- Z_pine_m[sample(nrow(Z_pine_m),size=1,replace=TRUE),, drop=FALSE]
         mate <- gamma*exp(-a*(Pmal[,2] - Pfem[,2] - d)^2)
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
         mate <- gamma*exp(-a*(Pmal[,2] - Pfem[,2] - d)^2)
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
   

##### Dispersal
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
   Z_pine_f <- Z_pine_f[exp(-S*(Z_pine_f[,2] - theta_pine)^2) > runif(length(Z_pine_f[,2]), 0, 1),]
   Z_pine_m <- Z_pine_m[exp(-S*(Z_pine_m[,2] - theta_pine)^2) > runif(length(Z_pine_m[,2]), 0, 1),]

#### Density-dependent mortality in pine
   if (length(Z_pine_f[,2]) > N_P/2){
      Z_pine_f <- Z_pine_f[sample(nrow(Z_pine_f), N_P/2, replace=FALSE),]
   }

   if (length(Z_pine_m[,2]) > N_P/2){
      Z_pine_m <- Z_pine_m[sample(nrow(Z_pine_m), N_P/2, replace=FALSE),]
   }

###### Selection in Oak
   Z_oak_f <- Z_oak_f[exp(-S*(Z_oak_f[,2] - theta_oak)^2) > runif(length(Z_oak_f[,2]), 0, 1),]
   Z_oak_m <- Z_oak_m[exp(-S*(Z_oak_m[,2] - theta_oak)^2) > runif(length(Z_oak_m[,2]), 0, 1),]

#### Density-dependent mortality in oak
if (length(Z_oak_f[,2]) > N_O/2){
   Z_oak_f <- Z_oak_f[sample(nrow(Z_oak_f), N_O/2, replace=FALSE),]
}

if (length(Z_oak_m[,2]) > N_O/2){
   Z_oak_m <- Z_oak_m[sample(nrow(Z_oak_m), N_O/2, replace=FALSE),]
}
   print(r)
}

