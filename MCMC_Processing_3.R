library(reshape2)
library(vecsets)
library(tidyverse)
library(MASS)
library(Boom)
library(coda)
library(bayestestR)
library(lhs)
library(MCMCglmm)
library(ggExtra)
library(ggplot2)

chain1 <- read.csv("Results_New_0.csv", header=TRUE)
colnames(chain1) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik")
chain1$chain <- rep("chain_1", length(chain1$iter))

chain2 <- read.csv("Results_New_1.csv", header=TRUE)
colnames(chain2) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik")
chain2$chain <- rep("chain_2", length(chain2$iter))

chain3 <- read.csv("Results_New_2.csv", header=TRUE)
colnames(chain3) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik")
chain3$chain <- rep("chain_3", length(chain3$iter))

chain4 <- read.csv("Results_New_3.csv", header=TRUE)
colnames(chain4) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik")
chain4$chain <- rep("chain_4", length(chain4$iter))

chain5 <- read.csv("Results_New_4.csv", header=TRUE)
colnames(chain5) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik")
chain5$chain <- rep("chain_5", length(chain5$iter))

chain6 <- read.csv("Results_New_5.csv", header=TRUE)
colnames(chain6) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik")
chain6$chain <- rep("chain_6", length(chain5$iter))



chains <- rbind(chain1,chain2,chain3,chain4,chain5,chain6)
chains.long <- melt(chains, id=c("chain", "iter"))

#plot the three unprocessed pmcmc chains
p <- ggplot(chains.long, aes(x = iter, y = value, group=chain)) + 
  geom_line(aes(color=chain)) +  theme_minimal()+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#999999", "red", "blue", "pink"))+ 
  facet_wrap(vars(variable), scales='free',   ncol = 3)
plot(p)


#post-process the chains
mcmc.chain1 <- mcmc(chain1[,-11])
mcmc.chain2 <- mcmc(chain2[,-11])
mcmc.chain3 <- mcmc(chain3[,-11])
mcmc.chain4 <- mcmc(chain4[,-11])
mcmc.chain5 <- mcmc(chain5[,-11])
mcmc.chain6 <- mcmc(chain6[,-11])
post.list <- list(mcmc.chain1,mcmc.chain2,mcmc.chain3,mcmc.chain4,mcmc.chain5,mcmc.chain6)
mcmc.list <- mcmc.list(post.list)
raftery.diag(mcmc.list, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
autocorr(mcmc.list, lags = c(0, 1, 5, 10, 50), relative=TRUE)
autocorr.plot(mcmc.list, lag.max=500, auto.layout = TRUE)

gelman.diag(mcmc.list, confidence = 0.95, transform=FALSE, autoburnin=TRUE,multivariate=TRUE)

processed <- window(mcmc.list, start=500, end=60000, thin=500)
processed <- data.frame(do.call(rbind, processed))

ci(processed$eta, ci=0.95, method = "ETI")

ci(processed$theta_oak, ci=0.95, method = "ETI")

ci(processed$theta_pine, ci=0.95, method = "ETI")

ci(processed$gamma, ci=0.95, method = "ETI")

ci(processed$SDa, ci=0.95, method = "ETI")

ci(processed$SDe, ci=0.95, method = "ETI")

ci(processed$SDs, ci=0.95, method = "ETI")

ci(processed$m, ci=0.95, method = "ETI")





processed.long <- melt(processed, id=c("iter"))
modes <- processed.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
modes <- modes[-9,]

p <- ggplot(processed.long, aes(x = value)) + theme_minimal()+
  geom_histogram(aes(y=..density..), colour="#00AFBB", fill="white")+
  geom_density(alpha=.1, fill="#00AFBB", color="#00AFBB")+
  geom_vline(data=modes, aes(xintercept=mode), color="red", size=0.75)+
  facet_wrap(vars(variable), scales='free',   ncol = 2)
plot(p)

p <- ggplot(processed, aes(x=gamma, y=eta)) + 
  geom_hex(bins = 10, aes(fill=..count..)) +
  geom_point(alpha = 0, show.legend=FALSE) +
  scale_fill_continuous(type = "viridis") +
  xlab(expression(paste("Ecological Selection ", (gamma)))) +  
  ylab(expression(paste("Habitat Preference",  (eta)))) +
  labs(fill='Count')+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.95, 0.80),
        legend.background = element_rect(fill = "gainsboro",size=0.1, linetype="solid", colour ="white"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"))+
  scale_x_continuous(limits=c(0, 0.1)) + scale_y_continuous(breaks=seq(0,0.3,0.05))
ggMarginal(p,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", size = 5, alpha=1)

p <- ggplot(processed, aes(x=gamma, y=eta))+
  geom_density_2d(size=1,aes(color = ..level..), show.legend=TRUE)+scale_color_viridis_c(name = "Density")+
  geom_point(alpha = 1, size=0.5,show.legend=FALSE) +
  xlab(expression(paste("Ecological Selection ", (gamma)))) +  
  ylab(expression(paste("Habitat Preference",  (eta)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.95, 0.80),
        legend.background = element_rect(fill = "gainsboro",size=0.1, linetype="solid", colour ="white"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"))+
  scale_x_continuous(limits=c(0, 0.1)) + scale_y_continuous(breaks=seq(0,0.3,0.05))
ggMarginal(p,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", size = 5, alpha=1)


p <- ggplot(processed, aes(x=gamma, y=eta))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75)+ 
  geom_density_2d(colour="black")+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste("Ecological Selection ", (gamma)))) +  
  ylab(expression(paste("Habitat Preference",  (eta)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.85, 0.50),
        legend.background = element_rect(fill = "gainsboro",size=0.1, linetype="solid", colour ="white"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"))+
  guides(fill = guide_legend(title = "Density"))+
  scale_x_continuous(limits=c(0, 0.1)) + scale_y_continuous(breaks=seq(0,0.3,0.05))
ggMarginal(p,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", size = 5, alpha=1)
