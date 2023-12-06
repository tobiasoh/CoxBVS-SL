library(ggplot2)
library(tidyr)
library(purrr)
library(vcd)

avg_model_size = function(result) {
  mPIP = colMeans(result$post.gamma)
  return( colSums(as.numeric(mPIP > 0.5)) )
  
}

model_size_plot = function(sim_scenarios, warmup) {
  num_runs = length(sim_scenarios["true"])
  
  avg_true = rep(0, num_runs)
  avg_empty = rep(0, num_runs)
  avg_partial_uniform = rep(0, num_runs)
  avg_non_uniform = rep(0, num_runs)
  avg_noise = rep(0, num_runs)
  
  #find avg 
  for (i in 1:num_runs) {
    avg_true = rowSums( sim_scenarios$true[[i]]$result$post.gamma[-(1:warmup),] > 0.5 )#avg_model_size(true[[i]]$result)
    avg_empty = rowSums( sim_scenarios$empty[[i]]$result$post.gamma[-(1:warmup),] > 0.5 ) #avg_model_size(empty[[i]]$result)
    avg_partial_non_uniform = rowSums( sim_scenarios$partial_non_uniform[[i]]$result$post.gamma[-(1:warmup),] > 0.5 )#avg_model_size(partial_non_uniform[[i]]$result)
    avg_partial_uniform = rowSums( sim_scenarios$partial_uniform[[i]]$result$post.gamma[-(1:warmup),] > 0.5 )#avg_model_size(partial_uniform[[i]]$result)
    avg_noise = rowSums( sim_scenarios$noise[[i]]$result$post.gamma[-(1:warmup),] > 0.5 ) #avg_model_size(noise[[i]]$result)
  }
  df = data.frame(true=avg_true, empty=avg_empty, partial_uniform = avg_partial_uniform, 
                  partial_non_uniform=avg_partial_non_uniform, noise=avg_noise)

  
  
  
  
  df_long = pivot_longer(df, cols=colnames(df), cols_vary="slowest")
  colnames(df_long) = c("sim_scenario", "model_size")
  
  return( df_long)
  
  ggplot(df_long, aes(x = sim_scenario, y = model_size)) +
    geom_boxplot() +
    labs(title = "Model size",
         x = "Simulation scenarios",
         y = "Model size") +
    theme(plot.title = element_text(hjust = 0.5))

  
  
  
}

#calculates sensitivity / specificity for each mcmc iter, then avg. over iterations
sensitivity_and_specificity = function(result, p, warmup) {
  mPIP = colMeans(result$post.gamma[-(1:warmup),])
  selected_vars = which( mPIP > 0.5) 
  true_vars = which( truePara$gamma == 1 )
  
  
  sensitivity = sum(true_vars %in% selected_vars) / length(true_vars)
  
  not_in_model = seq(1,p,1)[!(seq(1,p,1) %in% true_vars)]
  
  specificity = sum( !(not_in_model %in% selected_vars) ) / length(not_in_model)
  
  #specificity = sum( !(selected_vars %in% true_vars)) / (p - length(true_vars))#1 - sum(true_vars %in% selected_vars) / (p - length(true_vars))
  
  
  accuracy = ( sum(!(not_in_model %in% selected_vars)) + sum(true_vars %in% selected_vars) ) / p
  
  accuracy = sum( (mPIP > 0.5) == truePara$gamma ) / p
  
  return(list("sensitivity"=sensitivity,
              "specificity" = specificity,
              "accuracy" = accuracy))
  
}



#function that creates a table of sensitivity, specificity and accuracy. 
create_metric_table = function(sim_scenarios) {
  num_runs = length(sim_scenarios$true)
  
  
  df = data.frame(sensitivity=rep(0,num_runs),
                  specificity=rep(0,num_runs),
                  accuracy=rep(0,num_runs))
             
  
  
  df_metric = data.frame(scenario=c("true", "empty"), dat=c(data.frame(df), data.frame(df)))
  
  metric = list("true" = df,
                "empty" = data.frame(df),
                "partial_non_uniform" = data.frame(df),
                "partial_uniform" = data.frame(df),
                "noise" = data.frame(df))
  

  
  
  for (i in 1:num_runs) {
    #might be better to use purrr::map2 here 
    metric$true[i,] = sim_scenarios$true[[i]]$metrics
    metric$empty[i,] = sim_scenarios$empty[[i]]$metrics
    metric$partial_non_uniform[i,] = sim_scenarios$partial_non_uniform[[i]]$metrics
    metric$partial_uniform[i,] = sim_scenarios$partial_uniform[[i]]$metrics
    metric$noise[i,] = sim_scenarios$noise[[i]]$metrics
  }

  accuracy = data.frame(row.names = names(metric), lower = rep(0,5), mean=rep(0, 5),  upper=rep(0,5))
  sensitivity = data.frame(accuracy)
  specificity = data.frame(accuracy)
  
  for (scenario in names(metric)) {
    accuracy[scenario,"mean"] = mean(metric[[scenario]]$accuracy)
    CI =  t.test(metric[[scenario]]$accuracy)$conf.int
    accuracy[scenario, "lower"] = CI[1]
    accuracy[scenario, "upper"] = CI[2]
    
    sensitivity[scenario,"mean"] = mean(metric[[scenario]]$sensitivity)
    CI =  t.test(metric[[scenario]]$sensitivity)$conf.int
    sensitivity[scenario, "lower"] = CI[1]
    sensitivity[scenario, "upper"] = CI[2]
    
    
    specificity[scenario,"mean"] = mean(metric[[scenario]]$specificity)
    #CI =  quantile((metric[[scenario]]$specificity), 0.025)
    specificity[scenario, "lower"] = quantile((metric[[scenario]]$specificity), 0.025)
    specificity[scenario, "upper"] = quantile((metric[[scenario]]$specificity), 0.975)
    
    
  }
  
  
  return( list("sensitivity" = sensitivity, 
               "specificity" = specificity,
               "accuracy" = accuracy) )
  
}


get_simulations = function(path) {
  true = list()
  empty = list()
  partial_uniform = list()
  partial_non_uniform = list()
  noise = list()
  
  #path = "SimStudy/full_test/"
  #path = "/data/tobiasoh/SimStudy/full_sim"
  for (i in 1:20) {
    load(sprintf("%strue%d.RData", path, i))
    #simulation_result$metrics = sensitivity_and_specificity(simulation_result$result, p, warmup)
    #save(simulation_result, file=sprintf("%strue%d.RData", path, i))
    true = append(true, list(simulation_result))

    
    load(sprintf("%sempty%d.RData", path, i))
    #simulation_result$metrics = sensitivity_and_specificity(simulation_result$result, p, warmup)
    #save(simulation_result, file=sprintf("%sempty%d.RData", path, i))
    empty = append(empty, list(simulation_result))

    
    load(sprintf("%spartial_uniform%d.RData", path, i))
    #simulation_result$metrics = sensitivity_and_specificity(simulation_result$result, p, warmup)
    #save(simulation_result, file=sprintf("%spartial_uniform%d.RData", path, i))
    partial_uniform = append(partial_uniform, list(simulation_result))
    
    
    load(sprintf("%spartial_non_uniform%d.RData", path, i))
    #simulation_result$metrics = sensitivity_and_specificity(simulation_result$result, p, warmup)
    #save(simulation_result, file=sprintf("%spartial_non_uniform%d.RData", path, i))
    partial_non_uniform = append(partial_non_uniform, list(simulation_result))

    
    
    load(sprintf("%snoise%d.RData", path, i))
    #simulation_result$metrics = sensitivity_and_specificity(simulation_result$result, p, warmup)
    #save(simulation_result,file=sprintf("%snoise%d.RData", path, i))
    noise = append(noise, list(simulation_result))

    
  }
  
  
  sim_scenarios = list("true"=true, "empty"=empty, "partial_uniform"= partial_uniform,
                       "partial_non_uniform"=partial_non_uniform, "noise"=noise)
  
  return( sim_scenarios )
  
}





ibs_plot = function(sim_scenarios) {
  num_runs = length(sim_scenarios$true)
  true_ibs = rep(0, num_runs)
  empty_ibs = rep(0, num_runs)
  partial_uniform_ibs = rep(0, num_runs)
  partial_non_uniform_ibs = rep(0, num_runs)
  noise_ibs = rep(0, num_runs)
  
  for (i in 1:num_runs) {
    #brier_score = BrierScoreVectorised(sim_scenarios$true[[i]]$result$beta.p[-(1:warmup)], 
    #                                   sim_scenarios$true[[i]]$
    
    
    true_ibs[i] = mean( sim_scenarios$true[[i]]$ibs )
    empty_ibs[i] = mean( sim_scenarios$empty[[i]]$ibs )
    partial_uniform_ibs[i] = mean( sim_scenarios$partial_uniform[[i]]$ibs )
    partial_non_uniform_ibs[i] = mean( sim_scenarios$partial_non_uniform[[i]]$ibs )
    noise_ibs[i] = mean( sim_scenarios$noise[[i]]$ibs )
  }
  
  df = data.frame(true=true_ibs, empty=empty_ibs, partial_uniform=partial_uniform_ibs,
                  partial_non_uniform=partial_non_uniform_ibs,
                  noise=noise_ibs)
  
  df_long = pivot_longer(df, cols=colnames(df), cols_vary="slowest")
  colnames(df_long) = c("sim_scenario", "ibs")
  
  return( df_long)
  
  ggplot(df_long, aes(x = sim_scenario, y = ibs)) +
    geom_boxplot() +
    labs(title = "Integrated Brier score (0-5)",
         x = "Simulation scenarios",
         y = "IBS") +
    theme(plot.title = element_text(hjust = 0.5))
  
}


#for (i in 1:num_runs) {
#  print( sensitivity_and_specificity(sim_scenarios$empty[[i]]$result, p, warmup)$accuracy )
#}
