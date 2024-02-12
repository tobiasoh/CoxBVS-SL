source("/data/tobiasoh/plot_helper.R")

path = "/data/tobiasoh/SimStudy/real_sim_thin6/"

#scenario_names = c("empty", "true", "non_uniform_max", "non_uniform_max10_", 
#                   "partial_non_uniform", "partial_non_uniform10_")

#scenario_names = c("empty", "true", "partial_uniform", 
#                   "partial_non_uniform", "noise")

#scenario_names = c("empty", "line", "line3_", "line6_", "partial_uniform", 
#                  "true")

scenario_names = c("empty", "noise50_", "noise", "noise200_", "true")

sim_scen = get_simulations(path, scenario_names, num_runs=20)

warmup = ceiling( (sim_scen[[1]][[1]]$mcmcIterations / sim_scen$true[[1]]$thinning) / 2 )

post.gamma = posterior_gamma(sim_scen, warmup)


beta_cond_mean = beta_conditional_mean(sim_scen)



model_size = model_size_plot(sim_scen, warmup)


for (scen in names(sim_scen)) {
  for (i in 1:20) {
    sim_scen[[scen]][i]$metrics$recall = sim_scen[[scen]][i]$metrics$sensitivity / (sim_scen[[scen]][i]$metrics$sensitivity + 1 - sim_scen[[scen]][i]$metrics$specificity)
  }
}

metric_table = create_metric_table(sim_scen)

load("/data/tobiasoh/SimStudy/datasets/test_dataset.RData")
sim_scen = add_ibs(sim_scen, time_point=9, test_dataset)


ibs_info = ibs_plot2(sim_scen)





beta_mpm = beta_mpm_plot(sim_scen)


plot_info = list("metric_table" =metric_table, 
                 "ibs" = ibs_info, 
                 "model_size" = model_size,
                 "beta_mpm" = beta_mpm,
                 "brier_score" = list(),
                 "beta_conditional_mean" = beta_cond_mean,
                 "post.gamma" = post.gamma)



#save(plot_info, file=sprintf("%splot_info_feb5_sensitivity_non_uni.RData", path))

#save(plot_info, file=sprintf("%splot_info_feb5_initial.RData", path))

#save(plot_info, file=sprintf("%splot_info_feb3_sensitivity_uni.RData", path))

save(plot_info, file=sprintf("%splot_info_feb5_sensitivity_noise.RData", path))