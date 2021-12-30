#### SETUP
library(vecsets)
library(tidyverse)
set.seed(420)

#### SET PARAMETER VALUES
N_P = 216 # Pop. size in Pine
N_O = 1587 #Pop. size in Oak
a = 0.1 # Strength of mate preference
d = 2 # Offset match
n_pine = 62 # Number of nesting sites in Pine
n_oak = 453 # Number of nesting sites in Oak
SDa = 1.2 # Additive genetic variance to initialize pop.
SDe = 0.6 # Environmental variance 
SDs = 0.9 # Segregational variance
sex = 1 # Sexual dimorphism deviate
eta = 0.1 # Strength of habitat preference
M = 0.25*2 # Max prob of migration
M_P = M*0.12 # Probability of moving to pine given we are sampling oak birds
M_O = M*0.88 # Probability of moving to oak given we are sampling pine birds
theta_pine = 25.5 # Ecological optimum in Pine
theta_oak = 23.7 # Ecological optimum in Oak
S = 0.01 # Strength of selection

#### INITIALIZE THE POPULATION
BV_pine_f <- rnorm(N_P/2, mean = 23.5, sd = SDa)
BV_pine_m <- rnorm(N_P/2, mean = 23.5, sd = SDa)
BV_oak_f <- rnorm(N_O/2, mean = 23.5, sd = SDa)
BV_oak_m <- rnorm(N_O/2, mean = 23.5, sd = SDa)

Z_pine_f <- as.matrix(cbind(BV_pine_f, rnorm(N_P/2, mean = BV_pine_f-sex, sd = SDe)))
colnames(Z_pine_f) = NULL
Z_pine_m <- as.matrix(cbind(BV_pine_m, rnorm(N_P/2, mean = BV_pine_m+sex, sd = SDe)))
colnames(Z_pine_m) = NULL
Z_oak_f <- as.matrix(cbind(BV_oak_f, rnorm(N_O/2, mean = BV_oak_f-sex, sd = SDe)))
colnames(Z_oak_f) = NULL
Z_oak_m <- as.matrix(cbind(BV_oak_m, rnorm(N_O/2, mean = BV_oak_m+sex, sd = SDe)))
colnames(Z_oak_m) = NULL

#### RUN FOR R GENERATIONS; COLLECT DATA IN FOLLOWING VECTORS
R=500
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

#### BEGIN ITERATIONS

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



### MATE AND MAKE OFFSPRING IN PINE

   Pairs_pine <- data.frame(FemBV = rep(NA, n_pine), FemZ = rep(NA,n_pine) , MalBV = rep(NA, n_pine), MalZ = rep(NA,n_pine))
   for (i in 1:n_pine){ # pairs breed in pine
      repeat{
         Pfem <- Z_pine_f[sample(nrow(Z_pine_f),size=1,replace=TRUE),, drop=FALSE]
         Pmal <- Z_pine_m[sample(nrow(Z_pine_m),size=1,replace=TRUE),, drop=FALSE]
         mate <- exp(-a*(Pmal[,2] - Pfem[,2] - d)^2)
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

### MATE AND MAKE OFFSPRING IN OAK
   
   Pairs_oak <- data.frame(FemBV = rep(NA, n_oak), FemZ = rep(NA,n_oak) , MalBV = rep(NA, n_oak), MalZ = rep(NA,n_oak))
   for (i in 1:n_oak){
      repeat{
         Pfem <- Z_oak_f[sample(nrow(Z_oak_f),size=1,replace=TRUE),, drop=FALSE]
         Pmal <- Z_oak_m[sample(nrow(Z_oak_m),size=1,replace=TRUE),, drop=FALSE]
         mate <- exp(-a*(Pmal[,2] - Pfem[,2] - d)^2)
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
   

### MIGRATION - MIGRATION IS RANDOM IF ETA=0

   M_to_pine_f <- Z_oak_f[M_P*exp(-eta*(Z_oak_f[,2] - theta_pine)^2) / ( exp(-eta*(Z_oak_f[,2] - theta_pine)^2) + exp(-eta*(Z_oak_f[,2] - theta_oak)^2)) > runif(length(Z_oak_f[,2]), 0, 1),, drop=FALSE]


   M_to_pine_m <- Z_oak_m[M_P*exp(-eta*(Z_oak_m[,2] - theta_pine)^2) / ( exp(-eta*(Z_oak_m[,2] - theta_pine)^2) + exp(-eta*(Z_oak_m[,2] - theta_oak)^2)) > runif(length(Z_oak_m[,2]), 0, 1),, drop=FALSE]


   M_to_oak_f <- Z_pine_f[M_O*exp(-eta*(Z_pine_f[,2] - theta_oak)^2) / ( exp(-eta*(Z_pine_f[,2] - theta_oak)^2) + exp(-eta*(Z_pine_f[,2] - theta_pine)^2)) > runif(length(Z_pine_f[,2]), 0, 1),, drop=FALSE]

   M_to_oak_m <- Z_pine_m[M_O *exp(-eta*(Z_pine_m[,2] - theta_oak)^2) / ( exp(-eta*(Z_pine_m[,2] - theta_oak)^2) + exp(-eta*(Z_pine_m[,2] - theta_pine)^2)) > runif(length(Z_pine_m[,2]), 0, 1),, drop=FALSE]

### REMOVE MIGRANT FROM THEIR NATAL SUBPOPULATION
   
   Z_pine_f <- cbind(vsetdiff(Z_pine_f[,1], M_to_oak_f[,1]), vsetdiff(Z_pine_f[,2], M_to_oak_f[,2]))
   Z_pine_m <- cbind(vsetdiff(Z_pine_m[,1], M_to_oak_m[,1]), vsetdiff(Z_pine_m[,2], M_to_oak_m[,2]))
   Z_oak_f <- cbind(vsetdiff(Z_oak_f[,1], M_to_pine_f[,1]), vsetdiff(Z_oak_f[,2], M_to_pine_f[,2]))
   Z_oak_m <- cbind(vsetdiff(Z_oak_m[,1], M_to_pine_m[,1]), vsetdiff(Z_oak_m[,2], M_to_pine_m[,2]))

### PINE MALES AND FEMALES AFTER MIGRATION
   
   Z_pine_f <- rbind(Z_pine_f, M_to_pine_f)
   Z_pine_m <- rbind(Z_pine_m, M_to_pine_m)

### OAK MALES AND FEMALES AFTER MIGRATION
   
   Z_oak_f <- rbind(Z_oak_f, M_to_oak_f)
   Z_oak_m <- rbind(Z_oak_m, M_to_oak_m)

### ENVOKE SELECTION IN PINE
   
   Z_pine_f <- Z_pine_f[exp(-S*(Z_pine_f[,2] - theta_pine)^2) > runif(length(Z_pine_f[,2]), 0, 1),]
   Z_pine_m <- Z_pine_m[exp(-S*(Z_pine_m[,2] - theta_pine)^2) > runif(length(Z_pine_m[,2]), 0, 1),]

   ## Do not let population exceed carrying capacity
   if (length(Z_pine_f[,2]) > N_P/2){
      Z_pine_f <- Z_pine_f[sample(nrow(Z_pine_f), N_P/2, replace=FALSE),]
   }
   if (length(Z_pine_m[,2]) > N_P/2){
      Z_pine_m <- Z_pine_m[sample(nrow(Z_pine_m), N_P/2, replace=FALSE),]
   }

### ENVOKE SELECTION IN PINE
   
   Z_oak_f <- Z_oak_f[exp(-S*(Z_oak_f[,2] - theta_oak)^2) > runif(length(Z_oak_f[,2]), 0, 1),]
   Z_oak_m <- Z_oak_m[exp(-S*(Z_oak_m[,2] - theta_oak)^2) > runif(length(Z_oak_m[,2]), 0, 1),]

   ## Do not let population exceed carrying capacity
   if (length(Z_oak_f[,2]) > N_O/2){
   Z_oak_f <- Z_oak_f[sample(nrow(Z_oak_f), N_O/2, replace=FALSE),]
   }

   if (length(Z_oak_m[,2]) > N_O/2){
   Z_oak_m <- Z_oak_m[sample(nrow(Z_oak_m), N_O/2, replace=FALSE),]
   }
   
   print(r)
}

