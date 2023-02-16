library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)

ISSJ_phenotypes <- read.csv("ISSJ_phenotypes.csv", header=TRUE)
Oak <- subset(ISSJ_phenotypes, Habitat.Type =="Oak")
Pine <- subset(ISSJ_phenotypes, Habitat.Type =="Pine")
mean(Pine$Ave_Nares..mm.)
mean(Oak$Ave_Nares..mm.)
var(Oak$Ave_Nares..mm.)
var(Pine$Ave_Nares..mm.)
Oak_M <- subset(Oak, Sex =="M")
Oak_F <- subset(Oak, Sex =="F")
Pine_M <- subset(Pine, Sex =="M")
Pine_F <- subset(Pine, Sex =="F")

mean(Oak_M$Ave_Nares..mm.)
mean(Oak_F$Ave_Nares..mm.)
mean(Pine_M$Ave_Nares..mm.)
mean(Pine_F$Ave_Nares..mm.)

var(Oak_M$Ave_Nares..mm.)
var(Oak_F$Ave_Nares..mm.)
var(Pine_M$Ave_Nares..mm.)
var(Pine_F$Ave_Nares..mm.)

var(c(Oak_M$Ave_Nares..mm.,Pine_M$Ave_Nares..mm. ))
var(c(Oak_F$Ave_Nares..mm.,Pine_F$Ave_Nares..mm. ))

length(Oak_M$Ave_Nares..mm.)
length(Oak_F$Ave_Nares..mm.)
length(Pine_M$Ave_Nares..mm.)
length(Pine_F$Ave_Nares..mm.)

bartlett.test(list(Oak_M$Ave_Nares..mm., Oak_F$Ave_Nares..mm., Pine_M$Ave_Nares..mm., Pine_F$Ave_Nares..mm.))

pooled_var = 1/(211+175+92+79 - 4) * (210*var(Oak_M$Ave_Nares..mm.) + 174*var(Oak_F$Ave_Nares..mm.) +
                                     91*var(Pine_M$Ave_Nares..mm.) + 78*var(Pine_F$Ave_Nares..mm.))

t.test(Oak_F$Ave_Nares..mm., Oak_M$Ave_Nares..mm.)
t.test(Pine_F$Ave_Nares..mm., Pine_M$Ave_Nares..mm.)
t.test(Pine_F$Ave_Nares..mm., Oak_F$Ave_Nares..mm.)
t.test(Pine_M$Ave_Nares..mm., Oak_M$Ave_Nares..mm.)
t.test(Pine_F$Ave_Nares..mm., Oak_M$Ave_Nares..mm.)

n <- length(Oak_F$Ave_Nares..mm.)
xbar <- mean(Oak_F$Ave_Nares..mm.)
s <- sd(Oak_F$Ave_Nares..mm.)

n <- length(Oak_M$Ave_Nares..mm.)
xbar <- mean(Oak_M$Ave_Nares..mm.)
s <- sd(Oak_M$Ave_Nares..mm.)

n <- length(Pine_F$Ave_Nares..mm.)
xbar <- mean(Pine_F$Ave_Nares..mm.)
s <- sd(Pine_F$Ave_Nares..mm.)

n <- length(Pine_M$Ave_Nares..mm.)
xbar <- mean(Pine_M$Ave_Nares..mm.)
s <- sd(Pine_M$Ave_Nares..mm.)


#calculate margin of error
margin <- qt(0.975,df=n-1)*s/sqrt(n)
#calculate lower and upper bounds of confidence interval
low <- xbar - margin
low
high <- xbar + margin
high



chi_upper <- qchisq(1-0.975,df=(211+175+92+79 - 4))
chi_lower <- qchisq(0.975,df=(211+175+92+79 - 4))
#calculate lower and upper bounds of confidence interval
low <- (211+175+92+79 - 4)*1.3/chi_lower
low
high <- (211+175+92+79 - 4)*1.3/chi_upper
high

chi_upper <- qchisq(1-0.975,df=(79 - 1))
chi_lower <- qchisq(0.975,df=(79 - 1))
#calculate lower and upper bounds of confidence interval
low <- (79 - 1)*1.53861/chi_lower
low
high <- (79 - 1)*1.53861/chi_upper
high

oak_plot <- ggplot(Oak, aes(x=Ave_Nares..mm., color=Sex, fill=Sex))+ ylim(0,0.5)+
  scale_x_continuous(limits = c(19, 30), breaks = seq(19,30))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=30)+
  geom_density(alpha=0.2, adjust=1)+
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Status),linetype="solid", size=1)+
  scale_color_manual(values=c("tomato3", "tomato4"))+
  scale_fill_manual(values=c("tomato3", "tomato4")) +
  labs(title="Oak Habitat",x="Bill Length (mm)", y = "Density")+theme_minimal()+
  geom_vline(xintercept=mean(Oak_M$Ave_Nares..mm.), color="tomato4", size=1, linetype=2)+
  geom_vline(xintercept=mean(Oak_F$Ave_Nares..mm.), color="tomato3", size=1, linetype=2)
plot(oak_plot)

pine_plot <- ggplot(Pine, aes(x=Ave_Nares..mm., color=Sex, fill=Sex)) +  ylim(0,0.5)+
  scale_x_continuous(limits = c(19, 30), breaks = seq(19,30))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=30)+
  geom_density(alpha=0.2, adjust=1)+
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Status),linetype="solid", size=1)+
  scale_color_manual(values=c("skyblue3", "skyblue4"))+
  scale_fill_manual(values=c("skyblue3", "skyblue4"))+
  labs(title="Pine Habitat",x="Bill Length (mm)", y = "Density")+theme_minimal()+
  geom_vline(xintercept=mean(Pine_M$Ave_Nares..mm.), color="skyblue4", size=1, linetype=2)+
  geom_vline(xintercept=mean(Pine_F$Ave_Nares..mm.), color="skyblue3", size=1, linetype=2)
plot(pine_plot)


plot <- plot_grid(oak_plot, pine_plot, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot(plot)

library(tikzDevice)
# tikz("empirical_data.tex", width = 5, height = 5)
# plot(plot)
# graphics.off()




ggplot(ISSJ_phenotypes, aes(x=Ave_Nares..mm., color=Habitat.Type, fill=Habitat.Type))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=20)+
  geom_density(alpha=0.2)+
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Status),linetype="solid", size=1)+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  scale_fill_manual(values=c("skyblue4", "tomato3"))+
  labs(title="Pine Habitat",x="Bill Length (mm)", y = "Frequency")+theme_minimal()

M <- c(0.03663594,
       0.061827148,
       0.469829663,
       0.026593346,
       0.040104373,
       0.15256591,
       0.030012967,
       0.164821222,
       0.070126309,
       0.047853452,
       0.068546425,
       0.111085012,
       0.023791194,
       0.03327019,
       0.044090991)
mean(M)


M <- c(0.469829663,
       0.164821222,
       0.111085012,
       0.061827148,
       0.040104373,
       0.15256591,
       0.070126309,
       0.047853452,
       0.068546425,
       0.03327019,
       0.044090991)
mean(M)

M <- c(0.032679,
0.045247384,
0.020142233,
0.01973583,
0.01380469,
0.066683876,
0.030644569,
0.051444297,
0.040853588,
0.043242299,
0.042801777,
0.052918821,
0.018488543,
0.011743302,
0.022859345,
0.034219304)
mean(M)



n <- length(M)
xbar <- mean(M)
s <- sd(M)

#calculate margin of error
margin <- qt(0.975,df=n-1)*s/sqrt(n)

#calculate lower and upper bounds of confidence interval
low <- xbar - margin
low


high <- xbar + margin
high

shape=9.4
scale=0.01

shape*scale
(shape-1)*scale

a <- rgamma(100000, shape=shape, scale=scale)
mean(a)
posterior.mode(mcmc(a),adjust=1)
hist(a)

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
eta = 0.01 # strength of habitat preference
m = 0.084 # max prob of migration
theta_pine =25.2 # ecological optimum in pine
theta_oak = 23.7	 # ecological optimum in oak
gamma = 0.3

exp(-gamma*(25.6 - theta_pine)^2)
exp(-gamma*(23.5 - theta_pine)^2)
exp(-gamma*(24.8 - theta_oak)^2)
exp(-gamma*(22.7 - theta_oak)^2)

2*m*H_P*exp(-eta*(24.8 - theta_pine)^2) / ( exp(-eta*(24.8 - theta_pine)^2) + exp(-eta*(24.8 - theta_oak)^2)) 
2*m*H_P*exp(-eta*(22.7 - theta_pine)^2) / ( exp(-eta*(22.7 - theta_pine)^2) + exp(-eta*(22.7 - theta_oak)^2)) 
2*m*H_O*exp(-eta*( 25.6 - theta_oak)^2) / ( exp(-eta*(25.6 - theta_oak)^2) + exp(-eta*(25.6 - theta_pine)^2)) 
2*m*H_O*exp(-eta*(23.5 - theta_oak)^2) / ( exp(-eta*(23.5 - theta_oak)^2) + exp(-eta*(23.5 - theta_pine)^2)) 

