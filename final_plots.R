library(ggplot2)
library(cowplot)

#load("SimStudy/full_sim/plot_info.RData")
load("./SimStudy/real_sim_thin6/plot_info.RData")

### mcmc diagnostics

#model size
load("SimStudy/real_sim_thin6/true1.RData")
mcmc_iterations = simulation_result$mcmcIterations / simulation_result$thinning

#MCMC diagnostics: model size

model_size = rep(0, mcmc_iterations)


for (i in 1:(mcmc_iterations)) {
  model_size[i] = sum( as.numeric( simulation_result$result$gamma.p[i,] != 0))# / length(simulation_result$result$gamma.p[i,])
  
}

df = data.frame(model_size)
true_model_size = sum( as.numeric(simulation_result$truePara$beta != 0) )# / length(simulation_result$truePara$beta)

ggplot(data=df, aes(x=1:mcmc_iterations, y=model_size)) + 
  geom_line() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line(aes(x=1:mcmc_iterations, y=true_model_size), linetype="dashed", color="red" ) + 
  ggtitle("Model size") +
  labs(x="MCMC iterations", y = "Model size") #+ 
#geom_text(aes(label=c("", "True model size")))
#geom_label_repel(aes(label),
#                 nudge_x = 1,
#                 na.rm = TRUE)

#add text for "true model size"






# model size, all 5 scenarios in one figure

scenarios = c("true", "empty", "noise", "partial_uniform", "partial_non_uniform")
plt = list()

#layout(t(matrix(c(0,1,1,2,2,0,3,3,4,4,5,5), nrow=6, ncol=2)))

for (scen in scenarios) {
  
  
  load(sprintf("SimStudy/real_sim_thin6/%s1.RData", scen))
  mcmc_iterations = simulation_result$mcmcIterations / simulation_result$thinning
  
  #MCMC diagnostics: model size
  
  model_size = rep(0, mcmc_iterations)
  
  
  for (i in 1:(mcmc_iterations)) {
    model_size[i] = sum( as.numeric( simulation_result$result$gamma.p[i,] != 0))# / length(simulation_result$result$gamma.p[i,])
    
  }
  
  df = data.frame(model_size)
  true_model_size = sum( as.numeric(simulation_result$truePara$beta != 0) )# / length(simulation_result$truePara$beta)
  
  plt[[scen]] = ggplot(data=df, aes(x=1:mcmc_iterations, y=model_size)) + 
    geom_line() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_line(aes(x=1:mcmc_iterations, y=true_model_size), linetype="dashed", color="red" ) + 
    ggtitle(sprintf("%s", scen)) +
    labs(x="MCMC iterations", y = "Model size") #+ 
  #geom_text(aes(label=c("", "True model size")))
  #geom_label_repel(aes(label),
  #                 nudge_x = 1,
  #                 na.rm = TRUE)
  
  #add text for "true model size"
  plt[[scen]]
  
  
}

plt$true = plt$true + ggtitle("True")# + scale_x_discrete(labels=c(15, ))
plt$empty = plt$empty + ggtitle("Empty")
plt$noise = plt$noise + ggtitle("Noise")
plt$partial_uniform = plt$partial_uniform + ggtitle("Partial uniform")
plt$partial_non_uniform = plt$partial_non_uniform + ggtitle("Partial non-uniform")

for (scen in scenarios) {
  plt[[scen]]
}

top_row = plot_grid(plt$true, plt$empty, plt$noise, ncol=3, nrow=1)

bottom_row = plot_grid(NULL,plt$partial_uniform, plt$partial_non_uniform, NULL, ncol=4, rel_widths=c(0.125, 3/8, 3/8, 0.125 )  )



plot_grid(top_row,bottom_row, ncol=1, nrow=2)#, rel_heights=c(0.1,1,1))

top_row = plot_grid(plt$true, plt$empty, ncol=1, nrow=2)
mid_row = plot_grid(plt$partial_uniform, plt$partial_non_uniform, ncol=1, nrow=2)
bottom_row = plot_grid(plt$noise, NULL)

plot_grid(top_row, mid_row, bottom_row, ncol=1, nrow=3)


plot_grid(plt$true, plt$empty, plt$partial_uniform, plt$partial_non_uniform, plt$noise, ncol=1, nrow=5, 
          scale=1.5)

## different way: two + two + 1


ggsave(plt$partial_non_uniform + ylim(15,25), 
       filename="./SimStudy/real_sim_thin6/plots/mcmc_model_size_p_non_uni.png", 
       width=6, height=2.5)


#log likelihood

plt_ll = list()
for (scen in scenarios) {
  load(sprintf("./SimStudy/real_sim_thin6/%s1.RData", scen))
  loglik = simulation_result$result$log.like
  plt_ll[[scen]] = ggplot(data = data.frame(loglik), aes(x=1:mcmc_iterations, y=loglik)) + 
    geom_line() + 
    ggtitle("Log likelihood") + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="MCMC iterations", y = "Log likelihood")
}



plt_ll$true = plt_ll$true + ggtitle("True")# + scale_x_discrete(labels=c(15, ))
plt_ll$empty = plt_ll$empty + ggtitle("Empty")
plt_ll$noise = plt_ll$noise + ggtitle("Noise")
plt_ll$partial_uniform = plt_ll$partial_uniform + ggtitle("Partial uniform")
plt_ll$partial_non_uniform = plt_ll$partial_non_uniform + ggtitle("Partial non-uniform")

ggsave(plt_ll$partial_non_uniform, filename="./SimStudy/real_sim_thin6/plots/mcmc_loglike_p_non_uni.png", width=6, height=3)

top_row = plot_grid(plt_ll$true, plt_ll$empty, plt_ll$noise, ncol=3, nrow=1)

bottom_row = plot_grid(NULL,plt_ll$partial_uniform, plt_ll$partial_non_uniform, NULL, ncol=4, rel_widths=c(0.125, 3/8, 3/8, 0.125 )  )



plot_grid(top_row,bottom_row, ncol=1, nrow=2)#, rel_heights=c(0.1,1,1))










load("./SimStudy/real_sim_thin6/plot_info_feb5_initial.RData")
ordering = c("empty", "true", "partial_uniform", "partial_non_uniform", "noise")
labels = c("Empty", "True", "Partial uniform", "Partial non-uniform", "Noise")

load("./SimStudy/real_sim_thin6/plot_info_feb4_sensitivity_non_uni.RData")
ordering = c("empty", "non_uniform_max10_", "non_uniform_max", "partial_non_uniform10_", "partial_non_uniform", "true")
labels = c("Empty", "Non-uniform max 10", "Non-uniform max 5", "Partial non-uniform 10", "Partial non-uniform 5",  "True")

load("./SimStudy/real_sim_thin6/plot_info_feb3_sensitivity_uni.RData")
ordering = c("empty", "line", "line6_", "line3_", "partial_uniform", "true")
labels = c("Empty", "94% removal", "86% removal", "75% removal", "50% removal", "True")

load("./SimStudy/real_sim_thin6/plot_info_feb5_sensitivity_noise.RData")
ordering = c("empty", "true", "noise50_", "noise", "noise200_")
labels = c("Empty", "True", "50% noise", "100% noise", "200% noise")


#model size plot

plot_info$model_size$sim_scenario = factor(plot_info$model_size$sim_scenario, 
                                           levels = ordering)


ggplot(plot_info$model_size, aes(x = sim_scenario, y = model_size)) +
  geom_boxplot() +
  labs(#title = "Model size",
       x = "",#"Simulation scenarios",
       y = "Model size") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels=labels) + 
  geom_line(aes(group=dataset), linetype=2, linewidth=0.3, alpha=0.5) +
  geom_point(aes(group=dataset), alpha=0.5, size=1.2) + 
  ylim(0,23)




#beta mpm plot


lower_quantile=0.025
upper_quantile=0.975

#3credible_intervals <- t(apply(simulation_result$result$beta.p, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))

CI_empty = t(apply(plot_info$beta_mpm$empty, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))

warmup = 2500
#beta_after_warmup = simulation_result$result$beta.p[-(1:warmup),]

#cred_int = t(apply(beta_after_warmup, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))
#posterior_means = colMeans(beta_after_warmup)

CI_true = t(apply(plot_info$beta_mpm$true, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))

CI_noise = t(apply(plot_info$beta_mpm$noise, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))
CI_partial_uniform = t(apply(plot_info$beta_mpm$partial_uniform, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))
CI_partial_non_uniform = t(apply(plot_info$beta_mpm$partial_non_uniform, 2, function(x) quantile(x, c(lower_quantile, upper_quantile))))
x_offset = 0.1

num_vars = 20
num_runs = 20

xlabels = c()
for (i in num_vars:1) {
  xlabels = c(xlabels, paste("X", toString(i), sep=""))
}

gfg = data.frame(x=num_vars:1, y=colMeans(plot_info$beta_mpm$true)[num_vars:1], low=cred_int[num_vars:1,1], up=cred_int[num_vars:1,2])
#gfg = data.frame(x=1:200, y=posterior_means, z= low=cred_int[,1], up=cred_int[,2])

trueBeta = simulation_result$truePara$beta

beta_plt = ggplot(data=gfg, aes(x,y)) + 
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +

  geom_linerange(aes(x - x_offset, ymin=CI_empty[1:num_vars,1], ymax= CI_empty[1:num_vars,2]), color="grey" ) + 
  geom_linerange(aes(x-2*x_offset, ymin=CI_true[1:num_vars,1], ymax=CI_true[1:num_vars,2]), color="grey") + 
  geom_linerange(aes(x + x_offset, ymin=CI_partial_uniform[1:num_vars,1], ymax=CI_partial_uniform[1:num_vars,2]), color="grey" ) + 
  geom_linerange(aes(x + 2*x_offset, ymin=CI_partial_non_uniform[1:num_vars,1], ymax=CI_partial_non_uniform[1:num_vars,2]), color="grey" ) + 
  geom_linerange(aes(x, ymin=CI_noise[1:num_vars,1], ymax=CI_noise[1:num_vars,2]), color="grey" ) + 
  geom_segment(aes(x=num_vars:1 - 4*x_offset,xend=num_vars:1 + 4*x_offset, y=trueBeta[1:num_vars], yend=trueBeta[1:num_vars], color="True coefficient value"), linewidth=1.2, lineend="round") + 
  

  
  
  geom_point(aes(x - 2*x_offset, y=colMeans(plot_info$beta_mpm$true)[1:num_vars], color="True")) + 
  geom_point(aes(x- x_offset, colMeans(plot_info$beta_mpm$empty)[1:num_vars], color="Empty")) +
  geom_point(aes(x, colMeans(plot_info$beta_mpm$noise)[1:num_vars], color="Noise")) +
  geom_point(aes(x+x_offset, colMeans(plot_info$beta_mpm$partial_uniform)[1:num_vars], color="Partial uniform")) +
  geom_point(aes(x+x_offset*2, colMeans(plot_info$beta_mpm$partial_non_uniform)[1:num_vars], color="Partial non-uniform")) +
  
  #scale_x_continuous(labels = xlabels, breaks=seq(1,num_vars), minor_breaks = seq(1,num_vars)) + 
  scale_x_continuous(labels=c(rev(c("X1", "X50", "X100", "X150", "X200")))) + 
  scale_y_continuous(minor_breaks=c(-1,-0.5, 0, 0.5, 1)) + 
  
  #  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title=element_blank(),
        legend.position="bottom") + 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x="", y = "β conditional posterior mean +/- SD of β") + 
  #ggtitle("β posterior means and 95% credible intervals") + 
  scale_colour_manual("", 
                      breaks = c("Empty", "True", "Partial uniform", "Partial non-uniform", "Noise", "True coefficient value"),
                      values = c("darkgreen", "darkred", "orange", "purple", "blue", "red" ),
                      guide = guide_legend(override.aes = list(
                        linetype = c(rep("blank", 5), "solid"),
                        shape = c(rep(16, 5), NA)))) + 
  
  coord_flip()





beta_plt



#manhattan plot
num_runs = 20
num_vars = 200

xlabels = c()
for (i in num_vars:1) {
  xlabels = c(xlabels, paste("X", toString(i), sep=""))
}

CI = lapply(setNames(vector("list", length = length(plot_info$post.gamma)), 
                                  names(plot_info$post.gamma)), function(x) matrix(0, ncol=2, nrow=num_vars))

post.gamma_mean = lapply(setNames(vector("list", length = length(plot_info$post.gamma)), 
                                  names(plot_info$post.gamma)), function(x) rep(0,num_runs))

for (scen in names(plot_info$post.gamma)) {
  for (i in 1:num_vars) {
    post.gamma_mean[[scen]][i] = mean(plot_info$post.gamma[[scen]][,i])
    CI[[scen]][i,] = post.gamma_mean[[scen]][i] + sd(plot_info$post.gamma[[scen]][,i])*c(-1,1)
  }
  
}
gfg = data.frame(x=num_vars:1, y=post.gamma_mean$true)


manhattan = ggplot(data=gfg, aes(x,y)) + 
  geom_hline(yintercept = 0.5, linewidth = 0.5, color = "black") +
  geom_vline(xintercept = 200-20 + 0.5, linewidth=0.5, color="black", linetype=5) + 
  
  #geom_linerange(aes(x - 2*x_offset, ymin=CI$empty[1:num_vars,1], ymax= CI$empty[1:num_vars,2]), color="grey" ) + 
  #geom_linerange(aes(x-x_offset, ymin=CI$true[1:num_vars,1], ymax=CI$true[1:num_vars,2]), color="grey") + 
  #geom_linerange(aes(x, ymin=CI$partial_uniform[1:num_vars,1], ymax=CI$partial_uniform[1:num_vars,2]), color="grey" ) + 
  #geom_linerange(aes(x + x_offset, ymin=CI$partial_non_uniform[1:num_vars,1], ymax=CI$partial_non_uniform[1:num_vars,2]), color="grey" ) + 
  #geom_linerange(aes(x + 2*x_offset, ymin=CI$noise[1:num_vars,1], ymax=CI$noise[1:num_vars,2]), color="grey" ) + 

  
  
  
  geom_point(aes(x - 2*x_offset, y=post.gamma_mean$empty, color="Empty")) + 
  geom_point(aes(x- x_offset, y=post.gamma_mean$true[1:num_vars], color="True")) +
  geom_point(aes(x, post.gamma_mean$noise50_[1:num_vars], color="Partial uniform")) +
  geom_point(aes(x+x_offset, post.gamma_mean$noise[1:num_vars], color="Partial non-uniform")) +
  geom_point(aes(x+x_offset*2, post.gamma_mean$noise200_[1:num_vars], color="Noise")) +
  
  scale_x_continuous(labels = rev(c("X1", "X50", "X100", "X150", "X200")))+# breaks=seq(1,num_vars), minor_breaks = seq(1,num_vars)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0,1)) + 

theme(plot.title = element_text(hjust = 0.5), 
      legend.title=element_blank(),
      legend.position="bottom") + 
  
  #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x="", y = "Posterior mean of the inclusion indicators γ") + #"γ posterior mean +/- SD of γ") + 
  #ggtitle("γ posterior means and 95% credible intervals") + 
  scale_colour_manual("", 
                      breaks = c("Empty", "True", "Partial uniform", "Partial non-uniform", "Noise"),# "True coefficient value"),
                      values = c("darkgreen", "darkred", "orange", "purple", "blue"),#, "red" ),
                      guide = guide_legend(override.aes = list(
                        linetype = rep("blank", 5),#, "solid"),
                        shape = rep(16, 5))), 
                      labels=labels)+#, NA)))) + 
  
  coord_flip()

manhattan





#ibs plot
plot_info$ibs$sim_scenario = factor(plot_info$ibs$sim_scenario, 
                                    levels = ordering)

ggplot(plot_info$ibs, aes(x = sim_scenario, y = ibs)) +
  geom_boxplot() +
  labs(#title = "Integrated Brier score",
       x = "Simulation scenarios",
       y = "IBS") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels=labels)


#png(filename="./SimStudy/real_sim_thin6/plots/model_size_p_uni.png" ,width=600, height=300)

## ibs plot w lines
source("./plot_helper.R")

km_models = kaplan_meier(20, "./SimStudy/datasets")

km_ibs = kaplan_meier_ibs(km_models, 9, "./SimStudy/datasets/test_dataset.RData")
# 
# load("./SimStudy/datasets/test_dataset.RData")
# km_test = survfit(Surv(time, status)~ 1, data=test_dataset)
# data_test = data.frame(time = test_dataset$time, status = test_dataset$status, lp = 1)
# 
# Brier_train <- riskRegression::Score(list("Brier_test" = km_test), 
#                                      formula = Surv(time, status) ~ 1, 
#                                      data = data_test, conf.int = FALSE, 
#                                      metrics = "brier", summary = "ibs", 
#                                      times = seq(0, time_point, 0.1))$Brier$score
# Brier_score <- Brier_train[Brier_train$model != "Null model", ]
# 
# km_ibs <- Brier_score$IBS[time_point]#which.max(Brier_score$times)]



km_df = data.frame(dataset=1:20, sim_scenario=rep("kaplan_meier",20), ibs=km_ibs)
km_df$dataset = rep(NA,20)
#km_df$ibs = rep(NA,20)

ibs_df = rbind(plot_info$ibs, km_df)
ibs_df$sim_scenario = factor(ibs_df$sim_scenario, levels = c("kaplan_meier", ordering))
data_group = ibs_df[!is.na(ibs_df$dataset),]

plot_ibs = ggplot(ibs_df, aes(x = sim_scenario, y = ibs, label=sim_scenario)) +
  geom_boxplot() +
  geom_line(data=data_group, aes(group=dataset), linetype=2, linewidth=0.3, alpha=0.5) +
  geom_point(data=data_group, aes(group=dataset), alpha=0.5, size=1.2) + 
  labs(#title = "Integrated Brier score",
    x = "", #"Simulation scenarios",
    y = "IBS") +
  theme(plot.title = element_text(hjust = 0.5))+ 
  #geom_segment(aes(x=0.75, xend=1.25, y=km_ibs, yend=km_ibs), color="red") + 
  scale_x_discrete(labels=c("Kaplan-Meier", labels)) +
  ylim(0,max(km_ibs))
#geom_text_repel()
#geom_jitter(color="black", alpha = 0.5) + 



plot_ibs




## time-dependent brier score plot
time_dep_plot = time_dep_brier()
time_dep_plot


## metric tables (accuracy, sensitivity, specificity)



num_scenarios = dim(plot_info$metric_table$sensitivity)[1]

all = rep("", num_scenarios)
sens = cbind(plot_info$metric_table$sensitivity[ordering,], "Sensitivity"=all)
spec = cbind( plot_info$metric_table$specificity[ordering,], "Specificity"=all )
acc = cbind( plot_info$metric_table$accuracy[ordering,],"Accuracy"=all )
rec = cbind( plot_info$metric_table$recall[ordering,],"Recall"=all )
dig = 3

for (i in 1:num_scenarios) {
  sens$Sensitivity[i] = paste0( round(sens$mean[i], dig), " (", round(sens$mean[i] - sens$lower[i],dig), ")")#, ", ", round(sens$upper[i], dig), ")")
  spec$Specificity[i] = paste0( round(spec$mean[i], dig), " (", round(spec$mean[i] - spec$lower[i], dig), ")")# ", ", round(spec$upper[i], dig), ")")
  acc$Accuracy[i] = paste0( round(acc$mean[i],dig), " (", round(acc$mean[i] - acc$lower[i], dig), ")")# ", ", round(acc$upper[i], dig), ")")
  rec$Recall[i] = paste0( round(rec$mean[i],dig), " (", round(rec$mean[i] - rec$lower[i], dig),")")# ", ", round(rec$upper[i], dig), ")")
}


library(knitr)
library(dplyr)
full_table = cbind(sens$Sensitivity, spec$Specificity, acc$Accuracy)#, rec$Recall)
rownames(full_table) = labels#rownames(sens)
colnames(full_table) = c("Sensitivity", "Specificity", "Accuracy")#, "Recall")

knitr::kable(full_table, format="latex", digits=3)
knitr::kable(select(spec,all), format="latex", digits=3)
knitr::kable(select(acc,all), format="latex", digits=3)

knitr::kable(plot_info$metric_table$sensitivity, format="latex", digits=3)
knitr::kable(plot_info$metric_table$specificity, format="latex", digits=3)
knitr::kable(plot_info$metric_table$accuracy, format="latex", digits=3)

merged_table = cbind(plot_info$metric_table$sensitivity, plot_info$metric_table$specificity, 
                     plot_info$metric_table$accuracy, plot_info$metric_table$recall)


knitr::kable(merged_table, format="latex", digits=3)



