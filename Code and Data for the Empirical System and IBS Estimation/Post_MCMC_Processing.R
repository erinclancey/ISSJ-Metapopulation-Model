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
library(latex2exp)
library(ggh4x)
library(cowplot)

chain1 <- read.csv("Chain1.csv", header=TRUE)
colnames(chain1) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDe", "SDs", "m", "loglik", "accept")
chain1<- chain1[, c("iter","eta","gamma", "theta_oak", "theta_pine", "SDe", "SDs", "m", "loglik", "accept")]
chain1$Chain <- rep("1", length(chain1$iter))

chain2 <- read.csv("Chain2.csv", header=TRUE)
colnames(chain2) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDe", "SDs", "m", "loglik", "accept")
chain2<- chain2[, c("iter","eta","gamma", "theta_oak", "theta_pine", "SDe", "SDs", "m", "loglik", "accept")]
chain2$Chain <- rep("2", length(chain2$iter))

chain3 <- read.csv("Chain3.csv", header=TRUE)
colnames(chain3) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDe", "SDs", "m", "loglik", "accept")
chain3<- chain3[, c("iter","eta","gamma", "theta_oak", "theta_pine", "SDe", "SDs", "m", "loglik", "accept")]
chain3$Chain <- rep("3", length(chain3$iter))

chain4 <- read.csv("Chain4.csv", header=TRUE)
colnames(chain4) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDe", "SDs", "m", "loglik", "accept")
chain4<- chain4[, c("iter","eta","gamma", "theta_oak", "theta_pine", "SDe", "SDs", "m", "loglik", "accept")]
chain4$Chain <- rep("4", length(chain4$iter))

chain5 <- read.csv("Chain5.csv", header=TRUE)
colnames(chain5) <- c("iter","eta", "theta_oak", "theta_pine", "gamma", "SDe", "SDs", "m", "loglik", "accept")
chain5<- chain5[, c("iter","eta","gamma", "theta_oak", "theta_pine", "SDe", "SDs", "m", "loglik", "accept")]
chain5$Chain <- rep("5", length(chain5$iter))


chains <- rbind(chain1,chain2,chain3,chain4,chain5)
chains.long <- melt(chains, id=c("Chain", "iter"))
chains.long <- subset(chains.long, variable!="accept")
chains.long <- subset(chains.long, variable!="loglik")


plot_names1 <- as_labeller(c('eta' = "paste(eta)",
                             'gamma' = "paste(gamma)",
                             'theta_oak' = "paste(theta)[O]",
                             'theta_pine' = "paste(theta)[P]",
                             'SDe' = "paste(sigma)[e]",
                             'SDs' = "paste(sigma)[s]",
                             'm' = "m"),label_parsed)

#plot the three unprocessed pmcmc chains
C <- ggplot(chains.long, aes(x = iter, y = value, group=Chain)) + 
  geom_line(aes(color=Chain), size=0.25) +  theme_minimal()+
  scale_color_manual(values = c("#F0E442","#000000", "#009E73","#999999","#0072B2"))+ 
  facet_wrap(vars(variable),labeller = plot_names1, scales='free',   ncol = 2)+
  labs(x="Iteration",y="Parameter Value")+
  theme(panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))
C <- C + facetted_pos_scales(
  y = list(scale_y_continuous(limits=c(0, 0.3)),
           scale_y_continuous(limits=c(0, 0.3)),
           scale_y_continuous(limits=c(20, 30)),
           scale_y_continuous(limits=c(20, 30)),
           scale_y_continuous(limits=c(0.5, 1.5)),
           scale_y_continuous(limits=c(0.5, 1.5)),
           scale_y_continuous(limits=c(0.05, 0.15))))
plot(C)



#post-process the chains
mcmc.chain1 <- mcmc(chain1[,-10:-11])
mcmc.chain2 <- mcmc(chain2[,-10:-11])
mcmc.chain3 <- mcmc(chain3[,-10:-11])
mcmc.chain4 <- mcmc(chain4[,-10:-11])
mcmc.chain5 <- mcmc(chain5[,-10:-11])
post.list <- list(mcmc.chain1,mcmc.chain2, mcmc.chain3,mcmc.chain4,mcmc.chain5)
mcmc.list <- mcmc.list(post.list)
raftery.diag(mcmc.list, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
geweke.diag(mcmc.list, frac1=0.25, frac2=0.5)
gelman.diag(mcmc.list, confidence = 0.95, transform=TRUE, autoburnin=FALSE, multivariate=TRUE)
autocorr.diag(mcmc.list, lags = c(50, 100, 200, 500),relative=FALSE)
# Compile the posterior distributions
processed <- window(mcmc.list, start=200, end=30000, thin=200)
processed <- data.frame(do.call(rbind, processed))
length(processed$iter)

### delta_theta
D_theta <- processed$theta_pine-processed$theta_oak
D_theta_mode <- posterior.mode(mcmc(D_theta), adjust=1)
ci(D_theta, ci=0.95, method = "HDI")


# Generate 95% High Density Intervals for each parameter
ci(processed$eta, ci=0.95, method = "HDI")
ci(processed$theta_oak, ci=0.95, method = "HDI")
ci(processed$theta_pine, ci=0.95, method = "HDI")
ci(processed$gamma, ci=0.95, method = "HDI")
ci(processed$SDe, ci=0.95, method = "HDI")
ci(processed$SDs, ci=0.95, method = "HDI")
ci(processed$m, ci=0.95, method = "HDI")

#Calculate posterior modes
processed.long <- melt(processed, id=c("iter"))
processed.long <- subset(processed.long, variable!="loglik")
modes <- processed.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=2))
modes$hdi_low <- c(0.05,0,23.08,24.23,0.5,0.5,0.08)
modes$hdi_high <- c(0.3,0.02,24.08,26.78,1.09,0.89,0.09)
mode <- as.vector(modes$mode)

#Plot the marginal posterior distributions
plot_names2 <- as_labeller(c('eta' = "paste(eta)",
                             'gamma' = "paste(gamma)",
                             'theta_oak' = "paste(theta)[O]",
                             'theta_pine' = "paste(theta)[P]",
                             'SDe' = "paste(sigma)[e]",
                             'SDs' = "paste(sigma)[s]",
                             'm' = "m"),label_parsed)

P <- ggplot(processed.long, aes(x = value)) + theme_minimal()+
  geom_histogram(aes(y=..density..),position="identity", alpha=0.2, bins=30, colour="#333366", fill="#333366")+
  geom_density(alpha=.2, fill="#333366", color="#333366", adjust=2)+
  geom_vline(data=modes, aes(xintercept=mode), color="black", linewidth=1, linetype=2)+
  facet_wrap(vars(variable),labeller = as_labeller(plot_names2), scales='free', ncol = 2)+ labs(x="Parameter Value",y="Density")+
  theme(panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))+
  geom_segment(data=modes, aes(x=hdi_low,xend=hdi_high),y=0,yend=0,color="#333366",size=2,lineend="round")
plot(P)

# Plot the bivariate posterior distributions
J <- ggplot(processed, aes(x=gamma, y=eta))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75, adjust=2)+ 
  geom_density_2d(colour="black", adjust=2)+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste("Natural Selection ", (gamma)))) +  
  ylab(expression(paste("Habitat Preference",  (eta)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.85, 0.83),
        legend.background = element_rect(fill = "gainsboro",linewidth=0.1, linetype="solid", colour ="gainsboro"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"),
        legend.key.height= unit(0.25,'cm'),legend.key.width= unit(0.25,'cm'),
        legend.title = element_text(size=8),legend.text = element_text(size=6))+
  guides(fill = guide_legend(title = "Density"))+
  scale_x_continuous(limits=c(0, 0.03)) + scale_y_continuous(limits = c(0, 0.3))
J <- ggMarginal(J,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", bins=30, size = 5, alpha=1)


K <- ggplot(processed, aes(x=eta, y=theta_pine-theta_oak))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75, adjust=2)+ 
  geom_density_2d(colour="black", adjust=2)+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  ylab(expression(paste("Distance between Ecological Optima", (delta[theta])))) +  
  xlab(expression(paste("Habitat Preference",  (eta)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.856, 0.835),
        legend.background = element_rect(fill = "gainsboro",linewidth=0.1, linetype="solid", colour ="gainsboro"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"),
        legend.key.height= unit(0.25,'cm'),legend.key.width= unit(0.25,'cm'),
        legend.title = element_text(size=8),legend.text = element_text(size=6))+
  guides(fill = guide_legend(title = "Density"))+
  scale_y_continuous(limits=c(0, 4.0)) + scale_x_continuous(limits=c(0,0.3))
K <- ggMarginal(K,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", bins=30, size = 5, alpha=1)


L <- ggplot(processed, aes(x=theta_pine-theta_oak, y=gamma))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75, adjust=2)+ 
  geom_density_2d(colour="black", adjust=2)+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste("Distance between Ecological Optima", (delta[theta])))) +  
  ylab(expression(paste("Natural Selection ", (gamma)))) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.856, 0.82),
        legend.background = element_rect(fill = "gainsboro",linewidth=0.1, linetype="solid", colour ="gainsboro"),
        panel.background = element_rect(size=0.1, linetype="solid", color="white"),
        legend.key.height= unit(0.25,'cm'),legend.key.width= unit(0.25,'cm'),
        legend.title = element_text(size=8),legend.text = element_text(size=6))+
  guides(fill = guide_legend(title = "Density"))+
  scale_x_continuous(limits=c(0, 4)) + scale_y_continuous(limits = c(0,0.03))
L <- ggMarginal(L,data = processed, x=gamma, y=eta, type =  "histogram", color='black',margins = "both", bins=30, size = 5, alpha=1)



plot <- plot_grid(J,K,L, ncol = 3, nrow = 1, rel_heights=c(1,1,1), labels = c("A","B","C"))
plot(plot)

