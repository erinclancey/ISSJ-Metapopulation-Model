# Initial Setup
library(vecsets)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(latex2exp)
set.seed(420)


# Set Parameter Values
H_P = 0.12 # percent pine habitat
H_O = 0.88 # percent oak habitat
N_P = 216 # carrying capacity in pine  
N_O = 1587 # carrying capacity in oak
alpha = 0.1 # strength of mate preference 
delta = 2.06 # optimal mate offset match 
n_pine = 62 # number of nesting sites in pine
n_oak = 453 # number of nesting sites in oak
SDa = 1 # additive genetic standard deviation
SDe = 0.83 # environmental standard deviation
SDs = 0.6 # segregational standard deviation
sex = 1 # genetic factor causing sexual dimorphism
eta = 0.21 # strength of habitat preference
m = 0.084 # max prob of migration
theta_pine = 25.21 # ecological optimum in pine
theta_oak = 23.68 # ecological optimum in oak
gamma = 0.0057 # strength of selection



# Initialize population Breeding Values and Phenotypic Trait Values
start=theta_oak
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

R=1000 # Number of generations
Zbar_pine_F <- rep(NA, R)
Zbar_pine_M <- rep(NA, R)
Zbar_oak_F <- rep(NA, R)
Zbar_oak_M <- rep(NA, R)
Vp_pine_F <-rep(NA, R)
Vp_pine_M <-rep(NA, R)
Vp_oak_F <-rep(NA, R)
Vp_oak_M <-rep(NA, R)
Va_pine_F <-rep(NA, R)
Va_pine_M <-rep(NA, R)
Va_oak_F <-rep(NA, R)
Va_oak_M <-rep(NA, R)
h2 <- rep(NA, R)
N_pine <- rep(NA, R)
N_oak<- rep(NA, R)
BVbar_pine_F <- rep(NA, R)
BVbar_pine_M <- rep(NA, R)
BVbar_oak_F <- rep(NA, R)
BVbar_oak_M <- rep(NA, R)
Vp_POP<- rep(NA, R)
Va_POP<- rep(NA, R)

for (r in 1:R){
   
   N_pine[r] <- length(c(Z_pine_f[,1], Z_pine_m[,1]))
   N_oak[r] <- length(c(Z_oak_f[,1], Z_oak_m[,1]))

   Zbar_pine_F[r] <- mean(Z_pine_f[,2])
   Zbar_pine_M[r] <- mean(Z_pine_m[,2])
   Zbar_oak_M[r] <- mean(Z_oak_m[,2])
   Zbar_oak_F[r] <- mean(Z_oak_f[,2])
   
   BVbar_pine_M[r] <- mean(Z_pine_m[,1])
   BVbar_pine_F[r] <- mean(Z_pine_f[,1])
   BVbar_oak_M[r] <- mean(Z_oak_m[,1])
   BVbar_oak_F[r] <- mean(Z_oak_f[,1])

   Va_pine_M[r]<-var(Z_pine_m[,1])
   Va_pine_F[r]<-var(Z_pine_f[,1])
   Va_oak_M[r]<-var(Z_oak_m[,1])
   Va_oak_F[r]<-var(Z_oak_f[,1])
   Va_POP[r] <- var(c(Z_pine_f[,1], Z_pine_m[,1], Z_oak_f[,1], Z_oak_m[,1]))
   
   Vp_pine_F[r]<-var(Z_pine_f[,2])
   Vp_pine_M[r]<-var(Z_pine_m[,2])
   Vp_oak_F[r]<-var(Z_oak_f[,2])
   Vp_oak_M[r]<-var(Z_oak_m[,2])
   Vp_POP[r] <- var(c(Z_pine_f[,2], Z_pine_m[,2], Z_oak_f[,2], Z_oak_m[,2]))
   
   h2[r] <- var(c(Z_pine_f[,1], Z_pine_m[,1], Z_oak_f[,1], Z_oak_m[,1])) / var(c(Z_pine_f[,2], Z_pine_m[,2], Z_oak_f[,2], Z_oak_m[,2]))

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

   if (length(Z_pine_f[,2]) > N_P/2){
      Z_pine_f <- Z_pine_f[sample(nrow(Z_pine_f), N_P/2, replace=FALSE),]
   }

   if (length(Z_pine_m[,2]) > N_P/2){
      Z_pine_m <- Z_pine_m[sample(nrow(Z_pine_m), N_P/2, replace=FALSE),]
   }

   ###### Selection in Oak
   Z_oak_f <- Z_oak_f[exp(-gamma*(Z_oak_f[,2] - theta_oak)^2) > runif(length(Z_oak_f[,2]), 0, 1),]
   Z_oak_m <- Z_oak_m[exp(-gamma*(Z_oak_m[,2] - theta_oak)^2) > runif(length(Z_oak_m[,2]), 0, 1),]

   if (length(Z_oak_f[,2]) > N_O/2){
    Z_oak_f <- Z_oak_f[sample(nrow(Z_oak_f), N_O/2, replace=FALSE),]
    }

   if (length(Z_oak_m[,2]) > N_O/2){
    Z_oak_m <- Z_oak_m[sample(nrow(Z_oak_m), N_O/2, replace=FALSE),]
    }
   print(r)
}


########## Plot Figure B1

Gen <- seq(1,R, 1)
Zbar <- data.frame(cbind(Zbar_oak_F,Zbar_oak_M, Zbar_pine_F,Zbar_pine_M, Gen ))
var <- data.frame(cbind(Vp_pine_F, Vp_pine_M, Vp_oak_F, Vp_oak_M, Gen))

colors <- c("Zbar_pine_M" = "skyblue4","Zbar_pine_F" = "skyblue2", "Zbar_oak_M" = "tomato4","Zbar_oak_F" = "tomato2")
#plot means
p <- ggplot(Zbar, aes(x = Gen)) + theme_minimal()+ ggtitle("A")+ theme(axis.title.x = element_blank())+ theme(text = element_text(size = 10)) +theme(legend.text = element_text(size=6))
p <- p + geom_line(aes(y = Zbar_pine_F, color="Zbar_pine_F"))+ geom_line(aes(y = Zbar_pine_M, color="Zbar_pine_M"))
p <- p + geom_line(aes(y = Zbar_oak_F, color="Zbar_oak_F")) + labs(color = "") + 
  scale_color_manual(values = colors, labels = c("Oak Female","Oak Male","Pine Female","Pine Male"))+ 
  geom_line(aes(y = Zbar_oak_M, color="Zbar_oak_M"))
p <- p + geom_hline(yintercept = 25.57424, color="skyblue4", linetype=2)
p <- p + geom_hline(yintercept = 24.78057, color="tomato4", linetype=2)
p <- p + geom_hline(yintercept = 23.49823, color="skyblue2", linetype=2)
p <- p + geom_hline(yintercept = 22.70103, color="tomato2", linetype=2)
p <- p + scale_x_continuous() + scale_y_continuous(name = '$\\bar{z}$', breaks = seq(18,32,0.5), limits = c(22, 26))
legend <- cowplot::get_legend(p)
grid.newpage()
grid.draw(legend) 
p <- p+theme(legend.position="none")
A <- p
#plot variances
p <- ggplot(Zbar, aes(x = Gen)) + theme_minimal()+ ggtitle("B")+ theme(axis.title.x = element_blank())+ theme(text = element_text(size = 10)) +theme(legend.text = element_text(size=6))
p <- p + geom_line(aes(y = Vp_pine_F, color="Zbar_pine_F"))+ geom_line(aes(y = Vp_pine_M, color="Zbar_pine_M"))
p <- p + geom_line(aes(y = Vp_oak_F, color="Zbar_oak_F")) + labs(color = "") + geom_line(aes(y = Vp_oak_M, color="Zbar_oak_M"))+
  scale_color_manual(values = colors, labels = c("Oak Female","Oak Male","Pine Female","Pine Male"))
p <- p + geom_hline(yintercept = 1.290559, color="skyblue4", linetype=2)
p <- p + geom_hline(yintercept = 1.318115, color="tomato4", linetype=2)
p <- p + geom_hline(yintercept = 1.53861, color="skyblue2", linetype=2)
p <- p + geom_hline(yintercept = 1.177218, color="tomato2", linetype=2)
p <- p + scale_x_continuous() + scale_y_continuous(name = '$\\hat{\\mathcal{V}}(z)$', limits = c(0.5, 2))
p <- p+theme(legend.position="none")
B <- p

x.grob <- textGrob("Generation", gp=gpar(fontsize=13))
plot <- plot_grid(A, B, ncol = 1, nrow = 2, rel_heights=c(1.25,0.75))
plot2 <- grid.arrange(arrangeGrob(plot, bottom = x.grob))
plot3 <- plot_grid(plot2, legend, ncol = 2, rel_widths = c(0.98, 0.2))
plot3
# library(tikzDevice)
# tikz("equil.tex", width = 6, height = 3)
# plot3
# graphics.off()



