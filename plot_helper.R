library(ggplot2)
library(tidyr)
library(purrr)
library(vcd)
library(survival)
library(dplyr)

path = "/data/tobiasoh"
#path = "."
source(sprintf("%s/Data_Simulation.R", path))
load(sprintf("%s/Weibull_param.RData", path) )



avg_model_size = function(result) {
  mPIP = colMeans(result$post.gamma)
  return( colSums(as.numeric(mPIP > 0.5)) )
  
}

model_size_plot = function(sim_scenarios, warmup) {
  num_runs = length(sim_scenarios[[1]])
  
  
  
  #avg_true = rep(0, num_runs)
  #avg_empty = rep(0, num_runs)
  #avg_partial_uniform = rep(0, num_runs)
  #avg_non_uniform = rep(0, num_runs)
  #avg_noise = rep(0, num_runs)
  
  df <- data.frame( lapply(setNames(vector("list", length = length(sim_scenarios)), 
                            names(sim_scenarios)), function(x) rep(0, num_runs))
  )
  dataset = 1:num_runs
  df = cbind(df,dataset)
  
  df_ = data.frame(dataset = rep(0,num_runs*length(sim_scenarios)), sim_scenario = rep("true",num_runs*length(sim_scenarios)), model_size = rep(0,num_runs*length(sim_scenarios)))
  scenarios = names(sim_scenarios)
  
  #find avg 
  for (sim_scen in 1:length(sim_scenarios)) {
    
    for (i in 1:num_runs) {
      
      df_$model_size[num_runs*(sim_scen-1) + i] = rowSums( sim_scenarios[[sim_scen]][[i]]$result$post.gamma[-(1:warmup),] > 0.5 )
      df_$dataset[num_runs*(sim_scen-1) + i] = i
      df_$sim_scenario[num_runs*(sim_scen-1) + i] = scenarios[[sim_scen]]
      
      
      #df[[sim_scen]][i] = rowSums( sim_scenarios[[sim_scen]][[i]]$result$post.gamma[-(1:warmup),] > 0.5 )
      
      
    #   avg_true = rowSums( sim_scenarios$true[[i]]$result$post.gamma[-(1:warmup),] > 0.5 )#avg_model_size(true[[i]]$result)
    #   avg_empty = rowSums( sim_scenarios$empty[[i]]$result$post.gamma[-(1:warmup),] > 0.5 ) #avg_model_size(empty[[i]]$result)
    #   avg_partial_non_uniform = rowSums( sim_scenarios$partial_non_uniform[[i]]$result$post.gamma[-(1:warmup),] > 0.5 )#avg_model_size(partial_non_uniform[[i]]$result)
    #   avg_partial_uniform = rowSums( sim_scenarios$partial_uniform[[i]]$result$post.gamma[-(1:warmup),] > 0.5 )#avg_model_size(partial_uniform[[i]]$result)
    #   avg_noise = rowSums( sim_scenarios$noise[[i]]$result$post.gamma[-(1:warmup),] > 0.5 ) #avg_model_size(noise[[i]]$result)
    #   
    #   avg_line = rowSums( sim_scenarios$line[[i]]$result$post.gamma[-(1:warmup),] > 0.5 ) #avg_model_size(noise[[i]]$result)
    #   avg_non_uniform_max = rowSums( sim_scenarios$non_uniform_max[[i]]$result$post.gamma[-(1:warmup),] > 0.5 ) #avg_model_size(noise[[i]]$result)
     
    }
  }
    
  #df = data.frame(true=avg_true, empty=avg_empty, partial_uniform = avg_partial_uniform, 
  #                partial_non_uniform=avg_partial_non_uniform, noise=avg_noise,
  #                line=avg_line, non_uniform_max=avg_non_uniform_max)

  
  
  
  
  #df_long = pivot_longer(df, cols=colnames(df), cols_vary="slowest")
  #colnames(df_long) = c("sim_scenario", "model_size")
  
  return( df_)

  
  
  
}

#calculates sensitivity / specificity for each mcmc iter, then avg. over iterations
sensitivity_and_specificity = function(simulation_result, p, warmup) {
  
  mPIP = colMeans(simulation_result$result$post.gamma[-(1:warmup),])
  selected_vars = which( mPIP > 0.5) 
  true_vars = which( simulation_result$truePara$gamma == 1 )

  
  
  sensitivity = sum(true_vars %in% selected_vars) / length(true_vars)
  
  not_in_model = seq(1,p,1)[!(seq(1,p,1) %in% true_vars)]
  
  specificity = sum( !(not_in_model %in% selected_vars) ) / length(not_in_model)
  
  #specificity = sum( !(selected_vars %in% true_vars)) / (p - length(true_vars))#1 - sum(true_vars %in% selected_vars) / (p - length(true_vars))
  
  
  accuracy = ( sum(!(not_in_model %in% selected_vars)) + sum(true_vars %in% selected_vars) ) / p
  
  accuracy = sum( (mPIP > 0.5) == simulation_result$truePara$gamma ) / p
  
  recall = sensitivity / (sensitivity + 1 - specificity)
  
  return(list("sensitivity"=sensitivity,
              "specificity" = specificity,
              "accuracy" = accuracy,
              "recall" = recall))
  
}



#function that creates a table of sensitivity, specificity and accuracy. 
create_metric_table = function(sim_scenarios) {
  num_runs = length(sim_scenarios[[1]])
  
  
  df = data.frame(sensitivity=rep(0,num_runs),
                  specificity=rep(0,num_runs),
                  accuracy=rep(0,num_runs),
                  recall=rep(0,num_runs))

  metric <- lapply(setNames(vector("list", length = length(sim_scenarios)), 
                            names(sim_scenarios)), function(x) data.frame(df))
  
  #metric <- vector("list", length = length(sim_scenarios))
  names(metric) = names(sim_scenarios)
  #print(metric)

  #metric = list("true" = df,
  #              "empty" = data.frame(df),
  #              "partial_non_uniform" = data.frame(df),
  #              "partial_uniform" = data.frame(df),
  #              "noise" = data.frame(df),
  #              "line"= data.frame(df),
  #              "non_uniform_max" = data.frame(df))
  

  
  for (sim_scen in names(metric)) {  
    for (i in 1:num_runs) {
      
      metric[[sim_scen]][i,] = sim_scenarios[[sim_scen]][[i]]$metrics
      
      #might be better to use purrr::map2 here 
      #metric$true[i,] = sim_scenarios$true[[i]]$metrics
      #metric$empty[i,] = sim_scenarios$empty[[i]]$metrics
      #metric$partial_non_uniform[i,] = sim_scenarios$partial_non_uniform[[i]]$metrics
      #metric$partial_uniform[i,] = sim_scenarios$partial_uniform[[i]]$metrics
      #metric$noise[i,] = sim_scenarios$noise[[i]]$metrics
      
      #metric$line[i,] = sim_scenarios$noise[[i]]$metrics
      #metric$non_uniform_max[i,] = sim_scenarios$noise[[i]]$metrics
    }
  }

  accuracy = data.frame(row.names = names(metric), lower = rep(0,length(metric)), mean=rep(0, length(metric)),  upper=rep(0,length(metric)))
  sensitivity = data.frame(accuracy)
  specificity = data.frame(accuracy)
  recall = data.frame(accuracy)
  
  for (scenario in names(metric)) {
    accuracy[scenario,"mean"] = mean(metric[[scenario]]$accuracy)
    #CI =  t.test(metric[[scenario]]$accuracy)$conf.int
    accuracy[scenario, "lower"] = accuracy[scenario, "mean"] - sd(metric[[scenario]]$accuracy) #quantile((metric[[scenario]]$accuracy), 0.025)
    accuracy[scenario, "upper"] = accuracy[scenario, "mean"]  + sd(metric[[scenario]]$accuracy) #quantile((metric[[scenario]]$accuracy), 0.975)
    
    sensitivity[scenario,"mean"] = mean(metric[[scenario]]$sensitivity)
    #CI =  t.test(metric[[scenario]]$sensitivity)$conf.int
    sensitivity[scenario, "lower"] = sensitivity[scenario,"mean"] - sd(metric[[scenario]]$sensitivity)   #quantile((metric[[scenario]]$sensitivity), 0.025)
    sensitivity[scenario, "upper"] = sensitivity[scenario,"mean"] + sd(metric[[scenario]]$sensitivity)  #quantile((metric[[scenario]]$sensitivity), 0.975)
    
    
    specificity[scenario,"mean"] = mean(metric[[scenario]]$specificity)
    #CI =  quantile((metric[[scenario]]$specificity), 0.025)
    specificity[scenario, "lower"] = specificity[scenario,"mean"] - sd(metric[[scenario]]$specificity) #quantile((metric[[scenario]]$specificity), 0.025)
    specificity[scenario, "upper"] = specificity[scenario,"mean"] + sd(metric[[scenario]]$specificity) #quantile((metric[[scenario]]$specificity), 0.975)
    
    recall[scenario,"mean"] = mean(metric[[scenario]]$recall)
    #CI =  quantile((metric[[scenario]]$specificity), 0.025)
    recall[scenario, "lower"] = recall[scenario,"mean"] - sd(metric[[scenario]]$recall) #quantile((metric[[scenario]]$specificity), 0.025)
    recall[scenario, "upper"] = recall[scenario,"mean"] + sd(metric[[scenario]]$recall) #quantile((metric[[scenario]]$specificity), 0.975)
    
    
    
  }
  
  return( list("sensitivity" = sensitivity, 
               "specificity" = specificity,
               "accuracy" = accuracy,
               "recall" = recall) )
  
}


get_simulations = function(path, scenario_names, num_runs) {
  #true = list()
  #empty = list()
  #partial_uniform = list()
  #partial_non_uniform = list()
  #noise = list()
  #line = list()
  #non_uniform_max = list()
  
  sim_scenarios <- vector("list", length = length(scenario_names))
  names(sim_scenarios) = scenario_names

  sim_scenarios = setNames(replicate(length(scenario_names), list()), scenario_names)

  for (sim_scen in names(sim_scenarios)) {
    for (i in 1:num_runs) {
      load(sprintf("%s%s%d.RData", path, sim_scen, i))
      sim_scenarios[[sim_scen]] = c(sim_scenarios[[sim_scen]], list(simulation_result))
    }
  }
  
  
  
  return( sim_scenarios )
  
}




ibs_plot2 = function(sim_scenarios) {
  num_runs = length(sim_scenarios[[1]])
  true_ibs = rep(0, num_runs)
  empty_ibs = rep(0, num_runs)
  partial_uniform_ibs = rep(0, num_runs)
  partial_non_uniform_ibs = rep(0, num_runs)
  noise_ibs = rep(0, num_runs)
  line_ibs = rep(0, num_runs)
  non_uniform_max_ibs = rep(0,num_runs)
  
  df_ = data.frame(dataset = rep(0,num_runs*length(sim_scenarios)), sim_scenario = rep("true",num_runs*length(sim_scenarios)), ibs = rep(0,num_runs*length(sim_scenarios)))
  scenarios = names(sim_scenarios)
  
  for (scenario in 1:length(scenarios)) {
    for (i in 1:num_runs) {
      
      df_$ibs[num_runs*(scenario-1) + i] = mean( sim_scenarios[[scenarios[[scenario]]]][[i]]$ibs )
      df_$dataset[num_runs*(scenario-1) + i] = i
      df_$sim_scenario[num_runs*(scenario-1) + i] = scenarios[[scenario]]
    }
  }

  
  return( df_)
  
  
  
}



posterior_gamma = function(sim_scenarios, warmup) {
  num_runs = length(sim_scenarios[[1]])
  num_vars = dim(sim_scenarios[[1]]$result$beta.p)[2]
  
  scenarios = names(sim_scenarios)
  
  p = length(sim_scenarios[[1]][[1]]$truePara$beta)
  matr = matrix(0, nrow=num_runs, ncol = p)
  post.gamma = lapply(setNames(vector("list", length = length(sim_scenarios)), 
                       scenarios), function(x) matr)
  
  
  for (scenario in scenarios) {
    for (i in 1:num_runs) {
      
      post.gamma[[scenario]][i,] = colMeans(sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup),]) 
      
    }
  }
  
  
  return(post.gamma)
  
  
  
}
  




beta_mpm_plot = function(sim_scenarios) {
  num_runs = length(sim_scenarios[[1]])

  scenarios = names(sim_scenarios)

  warmup = ceiling( sim_scenarios[[1]][[1]]$mcmcIterations / (2*sim_scenarios[[1]][[1]]$thinning ) )
  p = dim(sim_scenarios[[1]][[1]]$result$beta.p)[2]

  #df = data.frame(sim_scenario=rep("none", num_runs*length(scenarios)),
  #                beta_mpm = matrix(0, nrow=p, ncol=num_runs*length(scenarios)))
  
  matr = matrix(0, nrow=num_runs, ncol = p)
  df = lapply(setNames(vector("list", length = length(sim_scenarios)), 
                  names(sim_scenarios)), function(x) matr)
  #df = list("true" = matr, "empty"=matr, "noise"=matr, "partial_uniform"=matr, "partial_non_uniform"=matr, "line"=matr, "non_uniform_max"=matr)
  
  
  for (scenario in scenarios) {
    for (i in 1:num_runs) {
      beta_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$beta[-(1:warmup),]) # instead of posterior mean of betas, better to use MPM betas here
      gamma_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup),])
      beta_mpm <- beta_mean / gamma_mean
      beta_mpm[gamma_mean <= 0.5] <- 0
      
      df[[scenario]][i,] = beta_mpm
    }
  }
  
  return(df)
  
  
}


beta_conditional_mean = function(sim_scenarios) {
  num_runs = length(sim_scenarios[[1]])
  
  scenarios = names(sim_scenarios)
  
  warmup = ceiling( sim_scenarios[[1]][[1]]$mcmcIterations / (2*sim_scenarios[[1]][[1]]$thinning ) )
  p = dim(sim_scenarios[[1]][[1]]$result$beta.p)[2]
  
  #df = data.frame(sim_scenario=rep("none", num_runs*length(scenarios)),
  #                beta_mpm = matrix(0, nrow=p, ncol=num_runs*length(scenarios)))
  
  matr = matrix(0, nrow=num_runs, ncol = p)
  df = lapply(setNames(vector("list", length = length(sim_scenarios)), 
                       names(sim_scenarios)), function(x) matr)
  #df = list("true" = matr, "empty"=matr, "noise"=matr, "partial_uniform"=matr, "partial_non_uniform"=matr, "line"=matr, "non_uniform_max"=matr)
  
  
  for (scenario in scenarios) {
    for (i in 1:num_runs) {
      #beta_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$beta[-(1:warmup),]) # instead of posterior mean of betas, better to use MPM betas here
      #gamma_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup),])
      #beta_mpm <- beta_mean / gamma_mean
      #beta_mpm[gamma_mean <= 0.5] <- 0
      # gamma_1 = which( sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup),] == 1, arr.ind=T )
      # beta_no_warmup = sim_scenarios[[scenario]][[i]]$result$beta.p[-(1:warmup),]
      # beta_conditional = beta_no_warmup[gamma_1]
      # 
      # #beta_no_warmup * sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup)]
      # 
      # 
      # print(dim(beta_no_warmup))
      # print(dim(gamma_1))
      # 
      # print(dim(beta_conditional))
      # beta_cond_mean = colMeans(beta_conditional)
      # 
      # print(dim(beta_cond_mean))
      # print(dim(df[[scenario]][i,]))
      # 
      # df[[scenario]][i,] = beta_cond_mean
      # 
      
      
      for (bet in 1:p) {
        gamma_1 = which( sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup),][,bet] == 1)
        if ( identical(gamma_1, integer(0)) ) {
          df[[scenario]][i,bet] = 0
        } else {
          beta_conditional = sim_scenarios[[scenario]][[i]]$result$beta.p[-(1:warmup),][,bet][gamma_1]
          df[[scenario]][i,bet] = mean(beta_conditional)
        }
        
        
      }
    }
  }
  
  return(df)
  
  
}


create_test_dataset = function(n, p, truePara, Surv.e, Surv.c) {
  sim_data = sim_data_fun(n=n, p=p, surv.e=Surv.e, surv.c=Surv.c, beta1.p=truePara$beta, beta2.p=truePara$beta, truePara$sigma)
  x.test = sim_data$Testdata[[1]]$X.test
  time.test = sim_data$Testdata[[1]]$time.test
  status.test = sim_data$Testdata[[1]]$status.test
  
  return( list( "X" = x.test,
                "time" = time.test,
                "status" = status.test))
}


add_ibs = function(sim_scenarios, time_point, test_dataset) {
  num_runs = length(sim_scenarios[[1]])
  x.test = test_dataset$X
  time.test = test_dataset$time
  status.test = test_dataset$status
  
 
  
  warmup = ceiling( sim_scenarios$true[[1]]$mcmcIterations / (2*sim_scenarios$true[[1]]$thinning ) )
  
  km_models = kaplan_meier(num_runs, "/data/tobiasoh/SimStudy/datasets")
  
  scenarios = names(sim_scenarios)

  for (scenario in scenarios) {
    for (i in 1:num_runs) {
      
      if (scenario != "kaplan_meier") {
        
      
      beta_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$beta[-(1:warmup),]) # instead of posterior mean of betas, better to use MPM betas here
      gamma_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup),])
      beta_mpm <- beta_mean / gamma_mean
      beta_mpm[gamma_mean <= 0.5] <- 0
      lp <- as.vector(sim_scenarios[[scenario]][[i]]$X.train %*% beta_mpm)
      times <- sim_scenarios[[scenario]][[i]]$survival_data$time.train
      status <- sim_scenarios[[scenario]][[i]]$survival_data$status.train
      data_train <- data.frame(time = times, status = status, lp = lp)
      bayes_train <- coxph(Surv(time, status) ~ lp, 
                           data = data_train, y = TRUE, x = TRUE)
      
      lp_test <- as.vector(x.test %*% beta_mpm)
      data_test <- data.frame(time = time.test, status = status.test, lp = lp_test)
      bayes_test <- coxph(Surv(time, status) ~ lp, 
                           data = data_test, y = TRUE, x = TRUE)
      
      }
      else {
        #data_train <- data.frame(time = times, status = status, lp = 1)
        data_test = data.frame(time = time.test, status = status.test, lp = 1)
        bayes_test = survfit(km_models[i], data_test)
      }

      
      library(riskRegression)
      Brier_train <- riskRegression::Score(list("Brier_test" = bayes_test), 
                                           formula = Surv(time, status) ~ 1, 
                                           data = data_test, conf.int = FALSE, 
                                           metrics = "brier", summary = "ibs", 
                                           times = seq(0, time_point, 0.1))$Brier$score#sort(unique(data_train$time)))$Brier$score
      #times = time_star)$Brier$score
      Brier_score <- Brier_train[Brier_train$model != "Null model", ]
      plot(Brier_score$Brier ~ Brier_score$times, type = "S", lty = 1, 
           xlab = "time", ylab = "Brier score", main = "riskRegression")
      
      sim_scenarios[[scenario]][[i]]$ibs <- Brier_score$IBS[time_point]#which.max(Brier_score$times)]
      
      
      
    }
    
  }
  return(sim_scenarios)
  
}



time_dep_brier = function() {
  
  # make a time-dependent brier score plot of all scenarios for dataset 1
  source("./Data_Simulation.R")
  
  load("./SimStudy/real_sim/true1.RData")
  true1 = simulation_result
  load("./SimStudy/real_sim/empty1.RData")
  empty1 = simulation_result
  load("./SimStudy/real_sim/noise1.RData")
  noise1 = simulation_result
  
  load("./SimStudy/real_sim/partial_uniform1.RData")
  partial_uniform1 = simulation_result
  load("./SimStudy/real_sim/partial_non_uniform1.RData")
  partial_non_uniform1 = simulation_result
  
  n = true1$n
  p = length(true1$truePara$beta)
  truebet = true1$truePara$beta
  sim_data = sim_data_fun(n=n, p=p, surv.e=Surv.e, surv.c=Surv.c, beta1.p=truebet, beta2.p=truebet, true1$truePara$sigma)
  x.test = sim_data$Testdata[[1]]$X.test
  time.test = sim_data$Testdata[[1]]$time.test
  status.test = sim_data$Testdata[[1]]$status.test

  
  warmup = ceiling( true1$mcmcIterations / (2*true1$thinning ) )
  
  
  scenarios = list("true" = true1, "empty" = empty1, "noise" =noise1, "partial_uniform"=partial_uniform1, "partial_non_uniform"=partial_non_uniform1)
  
  
  for (scenario in names(scenarios)) {
    
    beta_mean <- colMeans(scenarios[[scenario]]$result$beta.p[-(1:warmup),]) # instead of posterior mean of betas, better to use MPM betas here
    gamma_mean <- colMeans(scenarios[[scenario]]$result$gamma.p[-(1:warmup),])
    beta_mpm <- beta_mean / gamma_mean
    beta_mpm[gamma_mean <= 0.5] <- 0
    lp <- as.vector(scenarios[[scenario]]$X.train %*% beta_mpm)
    times <- scenarios[[scenario]]$survival_data$time.train
    status <- scenarios[[scenario]]$survival_data$status.train
    data_train <- data.frame(time = times, status = status, lp = lp)
    bayes_train <- coxph(Surv(time, status) ~ lp, 
                         data = data_train, y = TRUE, x = TRUE)
    
    lp_test <- as.vector(x.test %*% beta_mpm)
    data_test <- data.frame(time = time.test, status = status.test, lp = lp_test)
    bayes_test <- coxph(Surv(time.test, status.test) ~ lp, 
                        data = data_test, y = TRUE, x = TRUE)

    library(riskRegression)
    Brier_test <- riskRegression::Score(list("Brier_test" = bayes_test), 
                                         formula = Surv(time, status) ~ 1, 
                                         data = data_test, conf.int = FALSE, 
                                         metrics = "brier", summary = "ibs", 
                                         times = sort(unique(time.test)))$Brier$score
    #times = time_star)$Brier$score
    Brier_score <- Brier_test[Brier_test$model != "Null model", ]
    plot(Brier_score$Brier ~ Brier_score$times, type = "S", lty = 1, 
         xlab = "time", ylab = "Brier score", main = "riskRegression")
    
    scenarios[[scenario]]$brier_score = Brier_score$Brier
    
    
    
    
  }
  
  
  # Create a sample data frame with multiple lines
  ibs_df <- data.frame(
    true = scenarios$true$brier_score,
    empty = scenarios$empty$brier_score,
    noise = scenarios$noise$brier_score,
    partial_uniform = scenarios$partial_uniform$brier_score,
    partial_non_uniform = scenarios$partial_non_uniform$brier_score
  )
  #return( ibs_df )
  
  # Create a line plot with colors and a legend
  time_dep_plot = ggplot(ibs_df, aes(x=unique(sort(time.test)))) + #aes(true=true, empty=empty, noise=noise, partial_uniform=partial_uniform, partial_non_uniform=partial_non_uniform)) +
    geom_line(aes(y=empty, color="Empty")) + 
    geom_line(aes(y=true, color="True")) +
    geom_line(aes(y=partial_uniform, color="Partial uniform")) + 
    geom_line(aes(y=partial_non_uniform, color="Partial non-uniform")) + 
    geom_line(aes(y=noise, color="Noise")) + 
    scale_colour_manual("", 
                        breaks = c("Empty", "True", "Partial uniform", "Partial non-uniform", "Noise"),
                        values = c("darkgreen", "darkred", "orange", "purple", "blue")) + 
    labs(title = "Time-dependent Brier score",
         x = "Time",
         y = "Brier score") +
    #scale_color_manual(values = c("true" = "red", "empty" = "blue", "noise" = "green", "partial_uniform" = "black", "partial_non_uniform" = "yellow")) +
    
    #xlab("") + 
    theme_minimal() + 
    theme(legend.position = c(0.3,0.3), plot.title = element_text(hjust = 0.5))
  
  time_dep_plot
  
  return(time_dep_plot)
  
}




kaplan_meier = function(num_runs, path_datasets) {
  
  km_models = list()
  
  for (i in 1:num_runs) {
    load( sprintf("%s/dataset%d.RData", path_datasets, i) )
    km_fit <- survfit(Surv(time.train, status.train) ~ 1, data=dataset)
    #append(km_models, list(km_fit))
    km_models[[i]] = km_fit
    
  }
  
  return( km_models )
  
  
}
  

kaplan_meier_ibs = function(km_models, time_point, path_test_dataset) {
  load(path_test_dataset)
  
  x.test = test_dataset$X
  time.test = test_dataset$time
  status.test = test_dataset$status
  
  km_ibs = rep(0, length(km_models))
  data_test = data.frame(time = time.test, status = status.test, lp = 1)
  
  
  for (i in 1:length(km_models)) {
    #bayes_test = predict(km_models[[i]], newdata=Surv(time.test, status.test), type="survival")
    
    
    library(riskRegression)
    Brier_train <- riskRegression::Score(list("Brier_test" = km_models[[i]]), 
                                         formula = Surv(time, status) ~ 1, 
                                         data = data_test, conf.int = FALSE, 
                                         metrics = "brier", summary = "ibs", 
                                         times = seq(0, time_point, 0.1))$Brier$score#sort(unique(data_train$time)))$Brier$score
    #times = time_star)$Brier$score
    Brier_score <- Brier_train[Brier_train$model != "Null model", ]
    #plot(Brier_score$Brier ~ Brier_score$times, type = "S", lty = 1, 
    #     xlab = "time", ylab = "Brier score", main = "riskRegression")
    
   km_ibs[i] <- Brier_score$IBS[time_point]#which.max(Brier_score$times)]
    
  }
  return( km_ibs)
}

