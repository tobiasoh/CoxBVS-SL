source("/data/tobiasoh/plot_helper.R")

path = "/data/tobiasoh/SimStudy/RealData/result_basal_"
path = "/data/tobiasoh/RealData/march13/"

#path = "/data/tobiasoh/RealData/march17test/"
#path = "/data/tobiasoh/RealData/march20long/"
path = "/data/tobiasoh/RealData/april3test/"
path = "/data/tobiasoh/RealData/april4/"
path = "/data/tobiasoh/RealData/april12/"
path = "/data/tobiasoh/RealData/comb_april14/"


scenario_names = c("output_empty", "output_full")
#scenario_names = c("full_small", "full_verysmall", "full_normal")
#scenario_names = c("output_full")

#sim_scen = get_simulations(path, scenario_names, num_runs=20, real_data=T)
sim_scen = get_simulations(path, scenario_names, num_runs=20, real_data=T)
#load(sprintf())



warmup = ceiling( (sim_scen[[1]][[1]]$mcmcIterations / sim_scen[[1]][[1]]$thinning) / 2 )

#beta_cred_int = beta_cred_ints(sim_scen, warmup)
beta_cred_int = NULL
post.gamma = posterior_gamma(sim_scen, warmup)


beta_cond_mean = beta_conditional_mean(sim_scen)



model_size = model_size_plot(sim_scen, warmup)

#load("/data/tobiasoh/RealData/march14clin/dataset_w_clin.RData") #load objects X and meta
load("/data/tobiasoh/RealData/march22test/dataset_w_clin.RData") #load objects X and meta
#test_indices = which(!( rownames(X) %in% rownames(sim_scen[[1]][[1]]$X.train) ) )
load("/data/tobiasoh/RealData/datasplits/test_indices1.RData")
X.test = X[test_indices,]
test.time = meta[test_indices,]$time
test.status = as.integer(meta[test_indices,]$vital_status == "Dead")

test_dataset = list("X" = X.test,
                    "time" = test.time,
                    "status" = test.status)


print(cat("Dim beta: ", dim(sim_scen[[1]][[1]]$result$beta.p)))

sim_scen = add_ibs(sim_scen, time_point=365*10, test_dataset, time_step=30, real_data=T, X=X, meta=meta)



ibs_info = ibs_plot2(sim_scen)
#exit()


beta_mpm = beta_mpm_plot(sim_scen)




plot_info = list("ibs" = ibs_info, 
                 "model_size" = model_size,
                 "beta_mpm" = beta_mpm,
                 "brier_score" = list(),
                 "beta_conditional_mean" = beta_cond_mean,
                 "post.gamma" = post.gamma,
                 "beta_cred_int" = beta_cred_int)



#save(plot_info, file=sprintf("%splot_info_feb5_sensitivity_non_uni.RData", path))

#save(plot_info, file=sprintf("%splot_info_feb5_initial.RData", path))

#save(plot_info, file=sprintf("%splot_info_feb3_sensitivity_uni.RData", path))

save(plot_info, file=sprintf("%splot_info_real.RData", path))