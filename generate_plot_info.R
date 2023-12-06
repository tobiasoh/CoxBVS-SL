source("plot_helper.R")

path = "/data/tobiasoh/SimStudy/full_sim/"
sim_scen = get_simulations(path)


metric_table = create_metric_table(sim_scen)

ibs_info = ibs_plot(sim_scen)



warmup = ceiling( (sim_scen$true[[1]]$mcmcIterations / sim_scen$true[[1]]$thinning) / 2 )

model_size = model_size_plot(sim_scen, warmup)

plot_info = list("metric_table" = metric_table, 
                 "ibs" = ibs_info, 
                 "model_size" = model_size)



save(plot_info, file=sprintf("%splot_info.RData", path))