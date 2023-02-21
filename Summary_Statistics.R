library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)

Phenotypes <- read.csv("Island_Scrub_Jay_Phenotypes.csv", header=TRUE)
Oak <- subset(Phenotypes, Habitat =="Oak")
Pine <- subset(Phenotypes, Habitat =="Pine")
Oak_M <- subset(Oak, Sex =="M")
Oak_F <- subset(Oak, Sex =="F")
Pine_M <- subset(Pine, Sex =="M")
Pine_F <- subset(Pine, Sex =="F")

# Calculate the means for males and female in each subpopulation
mean(Oak_M$Ave_Nares_mm)
mean(Oak_F$Ave_Nares_mm)
mean(Pine_M$Ave_Nares_mm)
mean(Pine_F$Ave_Nares_mm)

# Calculate a 95% CI for male and female means in each subpopulation by finding the margin of error (MOE)
n <- length(Oak_F$Ave_Nares_mm)
xbar <- mean(Oak_F$Ave_Nares_mm)
s <- sd(Oak_F$Ave_Nares_mm)
MOE_Oak_F <- qt(0.975,df=n-1)*s/sqrt(n)

n <- length(Oak_M$Ave_Nares_mm)
xbar <- mean(Oak_M$Ave_Nares_mm)
s <- sd(Oak_M$Ave_Nares_mm)
MOE_Oak_M  <- qt(0.975,df=n-1)*s/sqrt(n)

n <- length(Pine_F$Ave_Nares_mm)
xbar <- mean(Pine_F$Ave_Nares_mm)
s <- sd(Pine_F$Ave_Nares_mm)
MOE_Pine_F <- qt(0.975,df=n-1)*s/sqrt(n)

n <- length(Pine_M$Ave_Nares_mm)
xbar <- mean(Pine_M$Ave_Nares_mm)
s <- sd(Pine_M$Ave_Nares_mm)
MOE_Pine_M<- qt(0.975,df=n-1)*s/sqrt(n)

# Test for differences between means in males and female in each subpopulation
mean(Pine_F$Ave_Nares_mm)-mean(Oak_F$Ave_Nares_mm)
t.test(Pine_F$Ave_Nares_mm, Oak_F$Ave_Nares_mm)
mean(Pine_M$Ave_Nares_mm) - mean(Oak_M$Ave_Nares_mm)
t.test(Pine_M$Ave_Nares_mm, Oak_M$Ave_Nares_mm)

# Calculate the variances for males and female in each subpopulation
var(Oak_M$Ave_Nares_mm)
var(Oak_F$Ave_Nares_mm)
var(Pine_M$Ave_Nares_mm)
var(Pine_F$Ave_Nares_mm)

# Calculate the sample size for males and female in each subpopulation
length(Oak_M$Ave_Nares_mm)
length(Oak_F$Ave_Nares_mm)
length(Pine_M$Ave_Nares_mm)
length(Pine_F$Ave_Nares_mm)

# Test for differences between variances in males and female in each subpopulation
bartlett.test(list(Oak_M$Ave_Nares_mm, Oak_F$Ave_Nares_mm, Pine_M$Ave_Nares_mm, Pine_F$Ave_Nares_mm))

# Calculate the composite variance for the population
pooled_var = 1/(211+175+92+79 - 4) * (210*var(Oak_M$Ave_Nares_mm) + 174*var(Oak_F$Ave_Nares_mm) +
                                        91*var(Pine_M$Ave_Nares_mm) + 78*var(Pine_F$Ave_Nares_mm))

# calculate lower and upper bounds of confidence interval for the composite variance
chi_upper <- qchisq(1-0.975,df=(211+175+92+79 - 4))
chi_lower <- qchisq(0.975,df=(211+175+92+79 - 4))
low <- (211+175+92+79 - 4)*1.3/chi_lower
high <- (211+175+92+79 - 4)*1.3/chi_upper

# Plot the male and female distributions in each habitat
oak_plot <- ggplot(Oak, aes(x=Ave_Nares_mm, color=Sex, fill=Sex))+ ylim(0,0.5)+
  scale_x_continuous(limits = c(19, 30), breaks = seq(19,30))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=30)+
  geom_density(alpha=0.2, adjust=1)+
  scale_color_manual(values=c("tomato3", "tomato4"))+
  scale_fill_manual(values=c("tomato3", "tomato4")) +
  labs(title="Oak Habitat",x="Bill Length (mm)", y = "Density")+theme_minimal()+
  geom_vline(xintercept=mean(Oak_M$Ave_Nares_mm), color="tomato4", size=1, linetype=2)+
  geom_vline(xintercept=mean(Oak_F$Ave_Nares_mm), color="tomato3", size=1, linetype=2)
plot(oak_plot)

pine_plot <- ggplot(Pine, aes(x=Ave_Nares_mm, color=Sex, fill=Sex)) +  ylim(0,0.5)+
  scale_x_continuous(limits = c(19, 30), breaks = seq(19,30))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=30)+
  geom_density(alpha=0.2, adjust=1)+
  scale_color_manual(values=c("skyblue3", "skyblue4"))+
  scale_fill_manual(values=c("skyblue3", "skyblue4"))+
  labs(title="Pine Habitat",x="Bill Length (mm)", y = "Density")+theme_minimal()+
  geom_vline(xintercept=mean(Pine_M$Ave_Nares_mm), color="skyblue4", size=1, linetype=2)+
  geom_vline(xintercept=mean(Pine_F$Ave_Nares_mm), color="skyblue3", size=1, linetype=2)
plot(pine_plot)


plot <- plot_grid(oak_plot, pine_plot, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot(plot)
