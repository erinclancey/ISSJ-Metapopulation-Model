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
library(cowplot)
library(grid)
library(latex2exp)
library(ggh4x)

chain1 <- read.csv("Final_COV_0.csv", header=TRUE)
colnames(chain1) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain1$chain <- rep("chain_1", length(chain1$iter))

chain2 <- read.csv("Final_COV_1.csv", header=TRUE)
colnames(chain2) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain2$chain <- rep("chain_2", length(chain2$iter))

chain3 <- read.csv("Final_COV_2.csv", header=TRUE)
colnames(chain3) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain3$chain <- rep("chain_3", length(chain3$iter))

chain4 <- read.csv("Final_COV_3.csv", header=TRUE)
colnames(chain4) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain4$chain <- rep("chain_4", length(chain4$iter))

chain5 <- read.csv("Final_COV_4.csv", header=TRUE)
colnames(chain5) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain5$chain <- rep("chain_5", length(chain5$iter))

chain6 <- read.csv("Final_COV_5.csv", header=TRUE)
colnames(chain6) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain6$chain <- rep("chain_6", length(chain5$iter))


chains <- rbind(chain1,chain2,chain3,chain4,chain5,chain6)
chains.long <- melt(chains, id=c("chain", "iter"))
chains.long <- subset(chains.long, variable!="accept")

plot_names1 <- as_labeller(c('eta' = "paste(eta)",
                             'theta_oak' = "paste(theta)[O]",
                             'theta_pine' = "paste(theta)[P]",
                             'gamma' = "paste(gamma)",
                             'SDa' = "paste(sigma)[a]",
                             'SDe' = "paste(sigma)[e]",
                             'SDs' = "paste(sigma)[s]",
                             'm' = "m",
                             'loglik' = "Log-Likelihood"),label_parsed)

#plot the three unprocessed pmcmc chains
C <- ggplot(chains.long, aes(x = iter, y = value, group=chain)) + 
  geom_line(aes(color=chain)) +  theme_minimal()+
  scale_color_manual(values = c("#F0E442","#000000", "#009E73","#999999", "#CC79A7","#0072B2"))+ 
  facet_wrap(vars(variable),labeller = plot_names1, scales='free',   ncol = 3)+
  theme(legend.position="none")+ labs(x="Iteration",y="Parameter Value")+
  theme(panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))
C <- C + facetted_pos_scales(
  y = list(scale_y_continuous(limits=c(0, 0.3)),
           scale_y_continuous(limits=c(20, 30)),
           scale_y_continuous(limits=c(20, 30)),
           scale_y_continuous(limits=c(0, 0.3)),
           scale_y_continuous(limits=c(0.5, 2.0)),
           scale_y_continuous(limits=c(0.5, 1.5)),
           scale_y_continuous(limits=c(0.5, 1.5)),
           scale_y_continuous(limits=c(0.05, 0.15)),
           scale_y_continuous(limits=c(0, 1700))))
plot(C)




#post-process the chains
mcmc.chain1 <- mcmc(chain1[,-11:-12])
mcmc.chain2 <- mcmc(chain2[,-11:-12])
mcmc.chain3 <- mcmc(chain3[,-11:-12])
mcmc.chain4 <- mcmc(chain4[,-11:-12])
mcmc.chain5 <- mcmc(chain5[,-11:-12])
mcmc.chain6 <- mcmc(chain6[,-11:-12])
post.list <- list(mcmc.chain1,mcmc.chain2,mcmc.chain3,mcmc.chain4,mcmc.chain5,mcmc.chain6)
mcmc.list <- mcmc.list(post.list)
raftery.diag(mcmc.list, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
autocorr(mcmc.list, lags = c(0, 1, 5, 10, 50), relative=TRUE)
autocorr.plot(mcmc.list, lag.max=500, auto.layout = TRUE)

acf(chain2$eta, lag.max = 500)
acf(chain2$theta_oak, lag.max = 500)
acf(chain2$theta_pine, lag.max = 500)
acf(chain2$gamma, lag.max = 500)
acf(chain2$SDa, lag.max = 500)
acf(chain2$SDe, lag.max = 500)
acf(chain2$SDs, lag.max = 500)
acf(chain2$m, lag.max = 500)

gelman.diag(mcmc.list, confidence = 0.95, transform=TRUE, autoburnin=TRUE, multivariate=TRUE)

processed_final <- window(mcmc.list, start=300, end=30000, thin=300)
processed_final <- data.frame(do.call(rbind, processed_final))

chain1 <- read.csv("Results_COV_0.csv", header=TRUE)
colnames(chain1) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain1$chain <- rep("chain_1", length(chain1$iter))
chain2 <- read.csv("Results_COV_1.csv", header=TRUE)
colnames(chain2) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain2$chain <- rep("chain_2", length(chain2$iter))
chain3 <- read.csv("Results_COV_2.csv", header=TRUE)
colnames(chain3) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain3$chain <- rep("chain_3", length(chain3$iter))
chain4 <- read.csv("Results_COV_3.csv", header=TRUE)
colnames(chain4) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain4$chain <- rep("chain_4", length(chain4$iter))
chain5 <- read.csv("Results_COV_4.csv", header=TRUE)
colnames(chain5) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain5$chain <- rep("chain_5", length(chain5$iter))
chain6 <- read.csv("Results_COV_5.csv", header=TRUE)
colnames(chain6) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDa", "SDe", "SDs", "m", "loglik", "accept")
chain6$chain <- rep("chain_6", length(chain5$iter))
mcmc.chain1 <- mcmc(chain1[,-11:-12])
mcmc.chain2 <- mcmc(chain2[,-11:-12])
mcmc.chain3 <- mcmc(chain3[,-11:-12])
mcmc.chain4 <- mcmc(chain4[,-11:-12])
mcmc.chain5 <- mcmc(chain5[,-11:-12])
mcmc.chain6 <- mcmc(chain6[,-11:-12])
post.list <- list(mcmc.chain1,mcmc.chain2,mcmc.chain3,mcmc.chain4,mcmc.chain5,mcmc.chain6)
mcmc.list <- mcmc.list(post.list)

processed_short <- window(mcmc.list, start=300, end=10000, thin=300)
processed_short <- data.frame(do.call(rbind, processed_short))

processed <- rbind(processed_final, processed_short)
processed <- processed_final

ci(processed$eta, ci=0.95, method = "ETI")

ci(processed$theta_oak, ci=0.95, method = "ETI")

ci(processed$theta_pine, ci=0.95, method = "ETI")

ci(processed$gamma, ci=0.95, method = "ETI")

ci(processed$SDa, ci=0.95, method = "ETI")

ci(processed$SDe, ci=0.95, method = "ETI")

ci(processed$SDs, ci=0.95, method = "ETI")

ci(processed$m, ci=0.95, method = "ETI")


processed.long <- melt(processed, id=c("iter"))
processed.long <- subset(processed.long, variable!="loglik")
modes <- processed.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode <- as.vector(modes$mode)

plot_names2 <- c('eta' = "$\\eta$",
                'theta_oak' = "$\\theta_{O}$",
                'theta_pine' = "$\\theta_{P}$",
                'gamma' = "$\\gamma$",
                'SDa' = "$\\sigma_a$",
                'SDe' = "$\\sigma_e$",
                'SDs' = "$\\sigma_s$",
                'm' = "$m$")

plot_names2 <- as_labeller(c('eta' = "paste(eta)",
                             'theta_oak' = "paste(theta)[O]",
                             'theta_pine' = "paste(theta)[P]",
                             'gamma' = "paste(gamma)",
                             'SDa' = "paste(sigma)[a]",
                             'SDe' = "paste(sigma)[e]",
                             'SDs' = "paste(sigma)[s]",
                             'm' = "m"),label_parsed)

P <- ggplot(processed.long, aes(x = value)) + theme_minimal()+
  geom_histogram(aes(y=..density..),position="identity", alpha=0.2, bins=30, colour="#333366", fill="#333366")+
  geom_density(alpha=.2, fill="#333366", color="#333366")+
  geom_vline(data=modes, aes(xintercept=mode), color="black", linewidth=1, linetype=2)+
  facet_wrap(vars(variable),labeller = as_labeller(plot_names2), scales='free', ncol = 2)+ labs(x="Parameter Value",y="Density")+
  theme(panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))
plot(P)



J <- ggplot(processed, aes(x=gamma, y=eta))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75)+ 
  geom_density_2d(colour="black")+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste("Ecological Selection ", (gamma)))) +  
  ylab(expression(paste("Habitat Preference",  (eta)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.85, 0.50),
        legend.background = element_rect(fill = "gainsboro",linewidth=0.1, linetype="solid", colour ="gainsboro"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"))+
  guides(fill = guide_legend(title = "Density"))+
  scale_x_continuous(limits=c(0, 0.1)) + scale_y_continuous(breaks=seq(0,0.3,0.05))
J <- ggMarginal(J,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", bins=35, size = 5, alpha=1)
J
K <- ggplot(processed, aes(x=eta, y=theta_pine-theta_oak))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75)+ 
  geom_density_2d(colour="black")+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  ylab(expression(paste("Distance between Ecological Optima", (delta[theta])))) +  
  xlab(expression(paste("Habitat Preference",  (eta)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.12, 0.32),
        legend.background = element_rect(fill = "gainsboro",linewidth=0.1, linetype="solid", colour ="gainsboro"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"))+
  guides(fill = guide_legend(title = "Density"))+
  scale_y_continuous(limits=c(0, 5)) + scale_x_continuous(breaks=seq(0,0.3,0.05))
K <- ggMarginal(K,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", bins=35, size = 5, alpha=1)
K

L <- ggplot(processed, aes(x=theta_pine-theta_oak, y=gamma))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75)+ 
  geom_density_2d(colour="black")+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste("Distance between Ecological Optima", (delta[theta])))) +  
  ylab(expression(paste("Ecological Selection ", (gamma)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.12, 0.56),
        legend.background = element_rect(fill = "gainsboro",linewidth=0.1, linetype="solid", colour ="gainsboro"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"),
        legend.key.height= unit(0.5,'cm'),legend.key.width= unit(0.5,'cm'),
        legend.title = element_text(size=10),legend.text = element_text(size=8))+
  guides(fill = guide_legend(title = "Density"))+
  scale_x_continuous(limits=c(0, 5)) + scale_y_continuous(breaks=seq(0,0.1, 0.01))
L <- ggMarginal(L,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", bins=35, size = 5, alpha=1)
L

# tikz("bi_post.tex", width = 5, height = 6)
# ggMarginal(J,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", bins=35, size = 5, alpha=1)
# graphics.off()


plot <- plot_grid(J,K,L, ncol = 3, nrow = 1, rel_heights=c(1,1,1))
plot(plot)
