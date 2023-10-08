#this file generates plots for a simulation result. 
library(ggplot2)
library(ggrepel)
library(qqman)
library(plotrix)


#load("SimulationStudy/sparse_smallN/N=100_P=100/simulation_results98056840_emptyGraph_N100.RDdata")

load("SimulationStudy/sparse_smallN/N=100_P=200/simulation_results3612351_true_N100_P200.RDdata")
mcmc_iterations = simulation_result$mcmcIterations

#MCMC diagnostics: model size
mcmc_iterations = simulation_result$mcmcIterations
model_size = rep(0, mcmc_iterations+1)


for (i in 1:(mcmc_iterations+1)) {
  model_size[i] = sum( as.numeric( simulation_result$result$gamma.p[i,] != 0))# / length(simulation_result$result$gamma.p[i,])
  
}
#which. 
df = data.frame(model_size)
true_model_size = sum( as.numeric(simulation_result$truePara$beta != 0) )# / length(simulation_result$truePara$beta)

ggplot(data=df, aes(x=0:mcmc_iterations, y=model_size)) + 
  geom_line() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line(aes(x=0:mcmc_iterations, y=true_model_size), linetype="dashed", color="red" ) + 
  ggtitle("Model size") +
  labs(x="MCMC iterations", y = "Model size") #+ 
  #geom_text(aes(label=c("", "True model size")))
  #geom_label_repel(aes(label),
  #                 nudge_x = 1,
  #                 na.rm = TRUE)
  
#add text for "true model size"



#MCMC diagnostics: log-likelihood
loglik = simulation_result$result$log.like
ggplot(data = data.frame(loglik), aes(x=1:mcmc_iterations, y=loglik)) + 
         geom_line() + 
         ggtitle("Log likelihood") + 
         labs(x="MCMC iterations", y = "Log likelihood")



#check frequency of each inclusion indicator to be 1 for each variable
meanPIF = rep(0, length(simulation_result$truePara$beta)-7)
for (i in 8:length(simulation_result$truePara$beta)) {
  meanPIF[i - 7] = mean(simulation_result$result$post.gamma[-(1:warmup),i])
}

plot(meanPIF)
  
  



# manhattan plot: mean posterior inclusion probability
warmup = 300
post_gamma= simulation_result$result$post.gamma[-(1:warmup),] #removing warmup period


lower_quantile=0.025
upper_quantile=0.975
cred_int = t(apply(post_gamma, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))
mPIP = colMeans(post_gamma) #mean posterior inclusion probability

gfg = data.frame(x=1:20, y=mPIP[1:20], low=cred_int[1:20,1], up=cred_int[1:20,2])
gfg = data.frame(x=1:100, y=mPIP, low=cred_int[,1], up=cred_int[,2])

ggplot(data=gfg, aes(x,y)) + geom_point() +# + geom_errorbar(aes(ymin=low, ymax=up))# +
  theme(aspect.ratio=1/2, plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limit = c(0, 1), breaks = c(0, 0.5, 1)) + 
 # coord_cartesian(xlim = c(0, 10), ylim = c(0, 5)) +
labs(x="Covariate index", y="mPIP") + 
  ggtitle("Mean Posterior Inclusion Probabilities (mPIP)")
  #geom_point(aes(x, simulation_result$truePara$beta), color="red")#, linetype="dashed")




#beta credible interval plots


lower_quantile=0.025
upper_quantile=0.975

credible_intervals <- t(apply(simulation_result$result$beta.p, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))


beta_after_warmup = simulation_result$result$beta.p[-(1:warmup),]

cred_int = t(apply(beta_after_warmup, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))
posterior_means = colMeans(beta_after_warmup)

gfg = data.frame(x=1:50, y=posterior_means[1:50], low=cred_int[1:50,1], up=cred_int[1:50,2])
gfg = data.frame(x=1:100, y=posterior_means, low=cred_int[,1], up=cred_int[,2])

ggplot(data=gfg, aes(x,y)) + geom_point() + geom_errorbar(aes(ymin=low, ymax=up)) +
  geom_point(aes(x, simulation_result$truePara$beta[1:50]), color="red") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Covariate index", y = "Beta posterior") + 
  ggtitle("Beta posterior means and 95% credible intervals")




