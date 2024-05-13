library(coda)

source("plot_helper.R")

load("./RealData/mcmc1000/result_basal_kegg.RData")

load("./RealData/march14clin/output_empty18.RData")
load("./RealData/march29test/output_full1_2.RData")
load("./RealData/march29test/output_empty1.RData")

load("./RealData/april3test/output_full1.RData")
load("./RealData/april3test/output_empty1.RData")
load("./RealData/april4/output_empty1.RData")
load("./RealData/april12/output_empty1.RData")
load("./RealData/comb_april14/output_full1.RData")

load("./RealData/apr26test/output_full1.RData")
model_output$priorPara$a
model_output$priorPara$b

simulation_result = model_output

#MCMC diagnostics: model size
mcmc_iterations = 2000#simulation_result$mcmcIterations
model_size = rep(0, mcmc_iterations)


for (i in 1:(mcmc_iterations)) {
  model_size[i] = sum( as.numeric( simulation_result$result$gamma.p[i,] != 0))# / length(simulation_result$result$gamma.p[i,])
  
}

df = data.frame(model_size)
true_model_size = sum( as.numeric(simulation_result$truePara$beta != 0) )# / length(simulation_result$truePara$beta)


ggplot(data=df, aes(x=1:mcmc_iterations, y=model_size)) + 
  geom_line() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Model size") +
  labs(x="MCMC iterations (after thinning)", y = "Model size") #+ 
#geom_text(aes(label=c("", "True model size")))
#geom_label_repel(aes(label),
#                 nudge_x = 1,
#                 na.rm = TRUE)

#add text for "true model size"



#MCMC diagnostics: log likelihood
loglik = simulation_result$result$log.like
ggplot(data = data.frame(loglik), aes(x=1:mcmc_iterations, y=loglik)) + 
  geom_line() + 
  ggtitle("Log-likelihood") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="MCMC iterations (after thinning)", y = "Log-likelihood")


#



#effective sample size calculation

warmup = ceiling(mcmc_iterations/2)

num_runs = 20
path = "./RealData/comb_april14/"
model_name = c("output_empty", "output_full")
eff_size_loglike = list("output_empty" = c(), "output_full" = c())
eff_size_model_size = list("output_empty" = c(), "output_full" = c())
for (i in 1:num_runs) {
  for (mod in model_name) {
    load( sprintf("%s%s%i.RData", path, mod, i))
    loglike = model_output$result$log.like[-(1:warmup)]
    eff_size = coda::effectiveSize(coda::mcmc(model_output$result$log.like[-(1:warmup)], thin=50, start=50000))
    eff_size_loglike[[mod]] = c(eff_size_loglike[[mod]], eff_size)
    
    eff_size_model_size[[mod]] = c(eff_size_model_size[[mod]], coda::effectiveSize(coda::mcmc(rowSums(model_output$result$gamma.p[-(1:warmup),]), thin=1, start=50000)))
    
    
  }
  
  
}

#box plots of effective sample size, one for model size and one for loglikelihood

my_data <- data.frame(
  Group = rep(names(eff_size_loglike), each = 20),
  Value = unlist(eff_size_loglike)
  )

ggplot(my_data, aes(x=Group, y = Value)) +
  geom_boxplot() + 
  labs(title = "Log-likelihood",
    x = "",#"Simulation scenarios",
    y = "Effective sample size") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12)
  ) + 
  scale_x_discrete(labels=c("Empty graph", "KEGG graph"))
  
#ylim(0,23)


eff_model_size_df <- data.frame(
  Group = rep(names(eff_size_model_size), each = 20),
  Value = unlist(eff_size_model_size)
)

ggplot(eff_model_size_df, aes(x=Group, y = Value)) +
  geom_boxplot() + 
  labs(title = "Model size",
       x = "",#"Simulation scenarios",
       y = "Effective sample size") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12)
  ) + 
  scale_x_discrete(labels=c("Empty graph", "KEGG graph"))








warmup = ceiling(mcmc_iterations/2)
model_size = colMeans(simulation_result$result$post.gamma[-(1:warmup),])
sum(model_size > 0.5)






load("./RealData/march13/plot_info_real.RData")
load("./RealData/april4/plot_info_real.RData")
load("./RealData/april12/plot_info_real.RData")
load("./RealData/comb_april14/plot_info_real.RData")
load("./RealData/dataset_basal.RData")
load("./RealData/march14clin/dataset_w_clin.RData")
ordering = c("output_empty", "output_full")
labels = c("Empty", "KEGG graph")

#model size plot

plot_info$model_size$sim_scenario = factor(plot_info$model_size$sim_scenario, 
                                           levels = ordering)


ggplot(plot_info$model_size, aes(x = sim_scenario, y = model_size)) +
  geom_boxplot() +
  labs(#title = "Model size",
    x = "",#"Simulation scenarios",
    y = "Model size") +
  theme(plot.title = element_text(hjust = 0.5),
        #axis.text.x = element_text(size = 20)
        ) + 
  scale_x_discrete(labels=labels) + 
  geom_line(aes(group=dataset), linetype=2, linewidth=0.3, alpha=0.5) +
  geom_point(aes(group=dataset), alpha=0.5, size=1.2) #+ 
#ylim(0,23)




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
  geom_point(aes(group=dataset), alpha=0.5, size=1.2) #+ 
  #ylim(0,23)





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
num_runs = 20
km_models = km_models_real_data(num_runs=num_runs, path_to_test_indices="./RealData/datasplits", meta)

km_ibs = km_ibs_real_data(km_models, time_point = 10*365, meta, path_to_test_indices = "./RealData/datasplits")
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



km_df = data.frame(dataset=1:num_runs, sim_scenario=rep("kaplan_meier",num_runs), ibs=km_ibs)
#km_df$dataset = rep(NA,20)
#km_df$ibs = rep(NA,20)

cox_ibs = cox_clin_ibs(path_to_test_indices="./RealData/datasplits", X, meta, num_runs=num_runs)

cox_df = data.frame(dataset=1:num_runs, sim_scenario=rep("cox", num_runs), ibs=cox_ibs)

ibs_df = rbind(plot_info$ibs, km_df, cox_df)
ibs_df$sim_scenario = factor(ibs_df$sim_scenario, levels = c("kaplan_meier", "cox", ordering))
data_group = ibs_df[!is.na(ibs_df$dataset),]

plot_ibs = ggplot(ibs_df, aes(x = sim_scenario, y = ibs, label=sim_scenario)) +
  geom_boxplot() +
  geom_line(data=data_group, aes(group=dataset), linetype=2, linewidth=0.3, alpha=0.5) +
  geom_point(data=data_group, aes(group=dataset), alpha=0.5, size=1.2) + 
  labs(#title = "Integrated Brier score",
    x = "", #"Simulation scenarios",
    y = "IBS") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)
        )+
  #geom_segment(aes(x=0.75, xend=1.25, y=km_ibs, yend=km_ibs), color="red") + 
  scale_x_discrete(labels=c("Kaplan-Meier", "Cox clinical", labels))# +
#ylim(0,max(km_ibs))
#geom_text_repel()
#geom_jitter(color="black", alpha = 0.5) + 



plot_ibs



plot_info$beta_mpm$output_empty
vars_included_empty = which(plot_info$beta_mpm$output_empty[1,] != 0)
vars_included_full = which(plot_info$beta_mpm$output_full[1,] != 0)
#overlap_empty = c()
#overlap_full

for (i in 2:20) {
  vars_included_empty = which(plot_info$beta_mpm$output_empty[i,] %in% vars_included_empty)
  vars_included_full = which(plot_info$beta_mpm$output_full[i,] %in% vars_included_full)
  
  
  
}

load("./RealData/april4/plot_info_real_a35.RData")
load("./RealData/graph_w_clin_full.RData")
load("./RealData/march20long/dataset_w_clin.RData")
load("./RealData/comb_april14/plot_info_real.RData")
#load("./RealData/march22test/dataset_w_clin.RData")

load("./RealData/march22test/dataset_w_clin.RData")

#stability in variable selection of MPM

source("plot_helper.R")
num_selected = feature_stability(plot_info$post.gamma, 20, 491)
feature_stab_df = data.frame(x=1:length(num_selected$output_full), empty=num_selected$output_empty/20, full=num_selected$output_full/20)

ggplot(data=feature_stab_df, aes(x=x, y=full)) + 
  geom_bar(stat="identity")


#nodes = unique( which(G != 0,arr.ind=T)[,1] )
#not_nodes = which(!(1:nrow(G) %in% nodes)) 
#which(!(which(num_selected$output_full!= 0) %in% nodes))

#which(nodes %in% which(num_selected$output_empty!=0))


#not_nodes_selected = which(num_selected$output_full[not_nodes]>4)

#full_selected = which(num_selected$output_full>0)
#not_nodes_selected = full_selected[full_selected %in% not_nodes]


#num_selected$output_full[not_nodes][not_nodes_selected]




full_selected = which(num_selected$output_full >= 4) #selected at least 20% of the time (4 out of 20 runs)
selected_in_graph = full_selected[full_selected %in% vars_in_graph]
not_nodes = full_selected[-(which(full_selected %in% vars_in_graph))]


full_selected = which(num_selected$output_empty >= 4)
selected_in_graph = full_selected[full_selected %in% vars_in_graph]
not_nodes = full_selected[-(which(full_selected %in% vars_in_graph))]
if (identical(selected_in_graph, integer(0))) {
  selected_in_graph = c()
  not_nodes = full_selected
}

#find the beta est and hazard ratios of selected features not in graph. 
#for beta est and HR, use only those runs where that var was actually selected

arr_ind_selected = which(plot_info$beta_mpm$output_full[,not_nodes] != 0, arr.ind=T)
beta_est = rep(0, length(not_nodes))
std_err = rep(0, length(not_nodes))
post.gamma = rep(0, length(not_nodes))
post.gamma_std_err = rep(0, length(not_nodes))


for (i in 1:length(not_nodes)) {
  var = not_nodes[i]
  which_runs = which(plot_info$post.gamma$output_full[,var] >= 0.5)
  print(which_runs)
  #which_runs = 1:20
  rel_runs = plot_info$beta_mpm$output_full[which_runs, var]
  rel_runs_gamma = plot_info$post.gamma$output_full[which_runs, var]
  post.gamma[i] = mean(rel_runs_gamma)
  post.gamma_std_err[i] = sd(rel_runs_gamma)/sqrt(length(rel_runs_gamma))
  std_err[i] = sd(rel_runs)/sqrt(length(rel_runs))
  beta_est[i] = mean(rel_runs)
  
  #beta_est[i] = mean(plot_info$beta_conditional_mean$output_full[which_runs, var])
}

beta_est

exp(beta_est)

CI_low = beta_est - 1.96*std_err
CI_high = beta_est + 1.96*std_err

selected_ensembl =sapply(strsplit(colnames(X)[not_nodes], "[.]"), `[`, 1)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mapping_selected = getBM(filters="ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                  "external_gene_name", "description"),values=selected_ensembl,mart= mart)

#which(plot_info$post.gamma$output_empty[,not_nodes] > 0, arr.ind=T)


hr_table = data.frame("gene_name"=rep("", length(beta_est)),
                      "hr"= rep("", length(beta_est)))
for (i in 1:length(beta_est)) {
  #hr_table$hr[i] = paste(round(exp(beta_est[i]), digits=2), " (", round(exp(CI_low[i]), digits=2), ", ", round(exp(CI_high[i]), digits=2), ")", sep="")
  hr_table$hr[i] = paste(round(exp(beta_est[i]), digits=2), " (", round(exp(beta_est[i])*std_err[i], digits=3), ")", sep="")
  
  index = which(mapping_selected$ensembl_gene_id %in% selected_ensembl[i])
  hr_table$gene_name[i] = mapping_selected$external_gene_name[which(selected_ensembl[i] == mapping_selected$ensembl_gene_id)]
} 




library(knitr)
knitr::kable(hr_table, format="latex")




in_graph_ensembl =sapply(strsplit(colnames(X)[selected_in_graph], "[.]"), `[`, 1)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mapping_in_graph= getBM(filters="ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                  "external_gene_name", "description"),values=in_graph_ensembl,mart= mart)


not_selected_graph = vars_in_graph[which(!(vars_in_graph %in% selected_in_graph))]
not_selected_graph_ensembl = sapply(strsplit(colnames(X)[not_selected_graph], "[.]"), `[`, 1)
mapping_not_selected_graph = getBM(filters="ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                            "external_gene_name", "description"),values=not_selected_graph_ensembl,mart= mart)

most_selected = not_nodes[ num_selected$output_full[not_nodes] >= 7 ]




### boxplot of post.gamma of 4 out of graph selected vars by KEGG model





colMeans(plot_info$post.gamma$output_full[,not_nodes])
colMeans(plot_info$post.gamma$output_empty[,not_nodes])
num_scen = 2
scens = c("output_empty", "output_full")
gene_names = c("FBN3", "CSAG3" , "SPRR2E", "MAGEC2",)
num_genes = length(gene_names)

df_post.gamma = data.frame(sim_scenario = rep("none", num_runs*num_scen*num_genes), 
                           dataset = rep(0, num_runs*num_scen*num_genes), 
                           gexp_value = rep(0,num_runs*num_scen*num_genes),
                           gene_name = rep("none", num_runs*num_scen*num_genes)
                           )


for (gene_num in 1:num_genes) {
  

for (scen in 1:num_scen) {
  for (i in 1:num_runs) {
    idx = (scen -1)*num_genes*num_runs +(gene_num-1)*num_runs +  i
    df_post.gamma$sim_scenario[idx] = scens[scen]
    df_post.gamma$dataset[idx] = i
    df_post.gamma$gexp_value[idx] = plot_info$post.gamma[[scens[scen]]][i, not_nodes[gene_num]]
    df_post.gamma$gene_name[idx] = gene_names[gene_num]
      
      
  }
}
}
for (i in 1:length(df_post.gamma$gene_name)) {
  if (df_post.gamma$sim_scenario[i] == "output_full") {
    df_post.gamma$sim_scenario[i] = "KEGG graph"
  } else {
    df_post.gamma$sim_scenario[i] = "Empty graph"
  }
}


ggplot(df_post.gamma, aes(x=gene_name, y = gexp_value, fill=sim_scenario)) +
  #scale_fill_discrete(name = "", limits=c("Empty graph", "KEGG graph")) +
  geom_boxplot() + 
  labs(title = "",
       x = "",#"Simulation scenarios",
       y = "Posterior inclusion probability")+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = rel(1.8)),
        axis.title.y = element_text(size=rel(1.4)),
        legend.text = element_text(size=rel(1.3)),
        legend.title=element_blank()
  )# +

  #geom_boxplot(position=position_dodge(1))





#colMeans(plot_info$beta_mpm$output_full[,not_nodes])

#sum(full_selected > 415)


pl1 = plot_info
ms_tbl = pl1$model_size

ms_tbl$model_size = ms_tbl$model_size - plot_info$model_size$model_size


ms_tbl$model_size - plot_info$model_size$model_size



#look at vars selected when using additional criterion of credible interval not
#covering 0
vars_selected = var_selection_cred_int(plot_info)

nodes = unique( which(G != 0,arr.ind=T)[,1] )
which(vars_selected$output_full[[20]]!= 0) %in% nodes

which_selected = vars_selected$output_full
mean_num_selected = 0
for (i in 1:20) {
  mean_num_selected = mean_num_selected + sum(vars_selected$output_full[[i]])
  
}

mean_num_selected = mean_num_selected/20

var_sel = feature_stability(vars_selected, 20, 489)

feature_stab_df = data.frame(x=1:length(var_sel$output_full), empty=var_sel$output_empty/20, full=var_sel$output_full/20)

ggplot(data=feature_stab_df, aes(x=x, y=full)) + 
  geom_bar(stat="identity")



load("./RealData/march17test/full_small1.RData")

small_thin = model_output

load("./RealData/march17test/full_small2.RData")






#check different hyperparam b for 
load("./RealData/march17test/plot_info_real.RData")
stab = feature_stability(plot_info$beta_mpm, 2)



