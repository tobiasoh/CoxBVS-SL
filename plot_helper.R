library(ggplot2)
library(tidyr)
library(purrr)
library(vcd)
library(survival)
library(dplyr)

path = "/data/tobiasoh"
path = "."
source(sprintf("%s/Data_Simulation.R", path))
load(sprintf("%s/Weibull_param.RData", path) )



avg_model_size = function(result) {
  mPIP = colMeans(result$post.gamma)
  return( colSums(as.numeric(mPIP > 0.5)) )
  
}

model_size_plot = function(sim_scenarios, warmup) {
  num_runs = length(sim_scenarios[[1]])
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
      
      #df_$model_size[num_runs*(sim_scen-1) + i] = rowSums( sim_scenarios[[sim_scen]][[i]]$result$post.gamma[-(1:warmup),] > 0.5 )
      df_$model_size[num_runs*(sim_scen-1) + i] = sum( colMeans(sim_scenarios[[sim_scen]][[i]]$result$post.gamma[-(1:warmup),]) > 0.5 )#mean( rowSums( sim_scenarios[[sim_scen]][[i]]$result$post.gamma[-(1:warmup),] > 0.5 ) )
      df_$dataset[num_runs*(sim_scen-1) + i] = i
      df_$sim_scenario[num_runs*(sim_scen-1) + i] = scenarios[[sim_scen]]
  
  
    }
  }
  
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


  
  for (sim_scen in names(metric)) {  
    for (i in 1:num_runs) {
      
      metric[[sim_scen]][i,] = sim_scenarios[[sim_scen]][[i]]$metrics
      
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


get_simulations = function(path, scenario_names, num_runs, real_data=F) {

  sim_scenarios <- vector("list", length = length(scenario_names))
  names(sim_scenarios) = scenario_names

  sim_scenarios = setNames(replicate(length(scenario_names), list()), scenario_names)

  for (sim_scen in names(sim_scenarios)) {
    for (i in 1:num_runs) {
      load(sprintf("%s%s%d.RData", path, sim_scen, i))
      if (real_data) {
        simulation_result = model_output
      }
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
      #print(sim_scenarios[[scenarios[[scenario]]]][[i]]$ibs)
      
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
  p = ncol(sim_scenarios[[1]][[1]]$result$gamma.p)
  #print(num_runs)
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

  matr = matrix(0, nrow=num_runs, ncol = p)
  df = lapply(setNames(vector("list", length = length(sim_scenarios)), 
                       names(sim_scenarios)), function(x) matr)

  
  for (scenario in scenarios) {
    for (i in 1:num_runs) {

      
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


add_ibs = function(sim_scenarios, time_point, test_dataset, time_step=0.1, real_data=F, X=NA, meta=NA) {
  num_runs = length(sim_scenarios[[1]])
  x.test = test_dataset$X
  time.test = test_dataset$time
  status.test = test_dataset$status
  
 
  
  warmup = ceiling( sim_scenarios[[1]][[1]]$mcmcIterations / (2*sim_scenarios[[1]][[1]]$thinning ) )
  
  #km_models = kaplan_meier(num_runs, "/data/tobiasoh/SimStudy/datasets")
  
  scenarios = names(sim_scenarios)

  for (scenario in scenarios) {
    for (i in 1:num_runs) {
      if (real_data) {
        test_indices = sim_scenarios[[scenario]][[i]]$test_indices
        sim_scenarios[[scenario]][[i]]$X.train = X[-test_indices,]
        sim_scenarios[[scenario]][[i]]$survival_data$time.train = meta[-test_indices,]$time
        sim_scenarios[[scenario]][[i]]$survival_data$status.train = as.integer(meta[-test_indices,]$vital_status == "Dead")
        
        x.test = X[test_indices,]
        time.test = meta[test_indices,]$time
        status.test = as.integer(meta[test_indices,]$vital_status == "Dead")
        
        #print(x.test)
        #print(time.test)
        #print(status.test)
        
        
        #max_event_time = max(time.test[which(status.test==1)])
        
        #print("x.test, time.test, status.test")
        #print(sum(is.na(x.test)))
        #print(sum(is.na(time.test)))
        #print(sum(is.na(status.test)))
      }
      #print(sum(is.na(sim_scenarios[[scenario]][[i]]$X.train)))
      #print(sum(is.na(sim_scenarios[[scenario]][[i]]$survival_data$time.train)))
      #print(sum(is.na(sim_scenarios[[scenario]][[i]]$survival_data$status.train)))
      
      
      
      if (scenario != "kaplan_meier") {
        
      
      beta_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$beta.p[-(1:warmup),]) # instead of posterior mean of betas, better to use MPM betas here
      gamma_mean <- colMeans(sim_scenarios[[scenario]][[i]]$result$gamma.p[-(1:warmup),])
      beta_mpm <- beta_mean / gamma_mean
      beta_mpm[gamma_mean <= 0.5] <- 0
      
      #lp <- as.vector(sim_scenarios[[scenario]][[i]]$X.train %*% beta_mpm)
      #times <- sim_scenarios[[scenario]][[i]]$survival_data$time.train
      #status <- sim_scenarios[[scenario]][[i]]$survival_data$status.train
      #data_train <- data.frame(time = times, status = status, lp = lp)
      #bayes_train <- coxph(Surv(time, status) ~ lp, 
      #                     data = data_train, y = TRUE, x = TRUE)
      
      print(dim(x.test))
      print(length(beta_mpm))
      
      lp_test <- as.vector(x.test %*% beta_mpm)
      data_test <- data.frame(time = time.test, status = status.test, lp = lp_test)
      
      bayes_test <- coxph(Surv(time, status) ~ lp, 
                          data = data_test, y = TRUE, x = T)
      if (all(lp_test == 0)) {
        bayes_test = coxph(Surv(time, status) ~ 1, 
                           data = data_test, y = TRUE, x = T)
      }
      
      
      }
      else {
        #data_train <- data.frame(time = times, status = status, lp = 1)
        data_test = data.frame(time = time.test, status = status.test, lp = 1)
        bayes_test = survfit(km_models[i], data_test)
      }
      
      
      print(sum(is.na(bayes_test)))
      #print(data_test$lp)
      print(sum(is.na(beta_mpm)))
      print(sum(is.na(coef(bayes_test))))
      print(sum(which(is.na(coef(bayes_test)))))
      print(is.na(coef(bayes_test)))
      
      print(cat("Dataset nr ", i))
      print(cat( "Coef Bayes_test: ", coef(bayes_test)))
      
      max_event_time = max(time.test[which(status.test==1)])
      
      library(riskRegression)
      Brier_train <- riskRegression::Score(list("Brier_test" = bayes_test), 
                                           formula = Surv(time, status) ~ 1, 
                                           data = data_test, conf.int = FALSE, 
                                           metrics = "brier", summary = "ibs", 
                                           times = seq(0, max_event_time, time_step))$Brier$score#sort(unique(data_train$time)))$Brier$score
      #times = time_star)$Brier$score
      Brier_score <- Brier_train[Brier_train$model != "Null model", ]
      #plot(Brier_score$Brier ~ Brier_score$times, type = "S", lty = 1, 
      #     xlab = "time", ylab = "Brier score", main = "riskRegression")

      
      sim_scenarios[[scenario]][[i]]$ibs <- Brier_score$IBS[which.max(Brier_score$times)]
      #print(sprintf("IBS for %s, ds %d: %f", scenario, i, sim_scenarios[[scenario]][[i]]$ibs))
      
      
      
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
    
   #km_ibs[i] <- Brier_score$IBS[which(max(Brier_score$times))]#[time_point]#which.max(Brier_score$times)]
    km_ibs[i] <- Brier_score$IBS[length(Brier_score$IBS)]
  }
  return( km_ibs)
}



km_models_real_data = function(num_runs, path_to_test_indices, meta) {
  
  km_models = list()
  
  for (i in 1:num_runs) {
    load( sprintf("%s/test_indices%d.RData", path_to_test_indices, i) )
    time.train = meta[-test_indices,]$time
    print(time.train)
    status.train = as.integer( meta[-test_indices,]$vital_status == "Dead" )
    km_fit <- survfit(Surv(time.train, status.train) ~ 1)
    #append(km_models, list(km_fit))
    km_models[[i]] = km_fit
    
  }
  
  return( km_models )
}


km_ibs_real_data = function(km_models, time_point, meta, path_to_test_indices) {

  
  
  
  
  km_ibs = rep(0, length(km_models))

  
  for (i in 1:length(km_models)) {
    load(sprintf("%s/test_indices%d.RData", path_to_test_indices, i))
    time.test = meta[test_indices,]$time
    status.test = as.integer( meta[test_indices,]$vital_status == "Dead" )
    data_test = data.frame(time=time.test, status=status.test, lp=1)
    max_event_time = max(time.test[which(status.test==1)])
    #bayes_test = predict(km_models[[i]], newdata=Surv(time.test, status.test), type="survival")
    
    
    library(riskRegression)
    Brier_train <- riskRegression::Score(list("Brier_test" = km_models[[i]]), 
                                         formula = Surv(time, status) ~ 1, 
                                         data = data_test, conf.int = FALSE, 
                                         metrics = "brier", summary = "ibs", 
                                         times = seq(0, time_point, 30))$Brier$score#sort(unique(data_train$time)))$Brier$score
    #times = time_star)$Brier$score
    Brier_score <- Brier_train[Brier_train$model != "Null model", ]
    #plot(Brier_score$Brier ~ Brier_score$times, type = "S", lty = 1, 
    #     xlab = "time", ylab = "Brier score", main = "riskRegression")
    
    km_ibs[i] <- Brier_score$IBS[which.max(Brier_score$times)] #[time_point]#
    
  }
  return( km_ibs)
  
}


cox_clin_ibs = function(path_to_test_indices, X, meta, num_runs=20) {
  
  cox_ibs = rep(0, num_runs)
  
  for (i in 1:num_runs) {
    load(sprintf("%s/test_indices%d.RData", path_to_test_indices, i))
    
    time.train = meta[-test_indices,]$time
    status.train = as.integer( meta[-test_indices,]$vital_status == "Dead")
    age = X[-test_indices,490]
    treatment = X[-test_indices, 491]
    data_train = data.frame(time=time.train, status=status.train, age=age, treatment = treatment)
    model_train = coxph(Surv(time, status) ~ age + treatment, data=data_train, x=T)
    
    x.test = X[test_indices, c(490,491)]
    time.test = meta[test_indices,]$time
    status.test = as.integer( meta[test_indices,]$vital_status == "Dead" )
    lp_test = model_train$coefficients %*% t(x.test)
    data_test = data.frame(time=time.test, status=status.test, age=X[test_indices,490], treatment=X[test_indices,491])#lp=lp_test)
    max_event_time = max(time.test[which(status.test==1)])
    #bayes_test = predict(km_models[[i]], newdata=Surv(time.test, status.test), type="survival")
    
    
    
    library(riskRegression)
    Brier_train <- riskRegression::Score(list("Brier_test" = model_train), 
                                         formula = Surv(time, status) ~ 1, 
                                         data = data_test, conf.int = FALSE, 
                                         metrics = "brier", summary = "ibs", 
                                         times = seq(0, max_event_time, 30))$Brier$score#sort(unique(data_train$time)))$Brier$score
    #times = time_star)$Brier$score
    Brier_score <- Brier_train[Brier_train$model != "Null model", ]
    #plot(Brier_score$Brier ~ Brier_score$times, type = "S", lty = 1, 
    #     xlab = "time", ylab = "Brier score", main = "riskRegression")
    
    cox_ibs[i] <- Brier_score$IBS[which.max(Brier_score$times)] #[time_point]#
    
  }
  
  return( cox_ibs )
}

feature_stability = function(data, num_runs, num_vars, cutoff=0.5) {
  #num_runs = #nrow(data[[1]])
  num_selected = rep(0, num_vars)#rep(0, ncol(data[[1]]))
  
  num_selected <- lapply(setNames(vector("list", length(data)), names(data)), function(x) rep(0, length(num_selected)))
  
  
  for (scenario in names(data)) {
    for (i in 1:num_runs) {
      selected = as.integer( data[[scenario]][[i]] > cutoff )
      selected = as.integer( data[[scenario]][i,] > cutoff)
      num_selected[[scenario]] = num_selected[[scenario]] + selected
    }
  }
  
  return( num_selected )
  

  
}


beta_cred_ints = function(sim_scenarios, warmup, quantiles = c(0.025, 0.975) ) {
  num_vars = ncol( sim_scenarios[[1]][[1]]$result$beta.p)
  num_runs = length(sim_scenarios[[1]])
  cred_int = data.frame(lower = rep(0, num_vars), upper=rep(0,num_vars))
  
  cred_int = matrix(0, nrow=num_vars, ncol=2)
  
  cred_int_all_runs = rep(cred_int, num_runs)
  cred_int_all_runs = lapply(vector("list", length=num_runs), function(x) cred_int)
  
  #print( dim(cred_int_all_runs[[1]]) )
                        
  df = lapply(setNames(vector("list", length = length(sim_scenarios)), 
                       names(sim_scenarios)), function(x) cred_int_all_runs)
  
  #print( dim(cred_int_all_runs))
  
  for (scen in names(sim_scenarios)) {
    for (i in 1:num_runs) {
      beta_mcmc = sim_scenarios[[scen]][[i]]$result$beta.p[-c(1:warmup),]
      #print(dim((df[[scen]])))
      #print(dim(df))
      df[[scen]][[i]] = t( apply( beta_mcmc, 2, function(x) quantile(x, probs=quantiles) ) )
      
    }
  }
  #print("Final dim of df")
  #print(dim(df[[1]]))
  
  return( df )
  
}

if (F) {
  cred_int = data.frame(lower=rep(0,489), upper=rep(0,489))
  cred_int = matrix(0, nrow=489, ncol=2)
  cred_int_all_runs = rep(cred_int, 20)
  cred_int_all_runs = lapply(vector("list", length=20), function(x) cred_int)
  df = lapply(setNames(vector("list", length = 2), 
                       c("empty", "full")), function(x) cred_int_all_runs)
}


time_dep_brier_func = function(sim_scenarios, warmup) {
  
}


var_selection_cred_int = function(plot_info) {
  num_runs = nrow(plot_info$beta_mpm[[1]])
  vars_selected = plot_info$beta_cred_int
  for (scenario in names(vars_selected)) {
    for (i in 1:num_runs) {
      vars_mpm = as.integer(plot_info$beta_mpm[[scenario]][i,] != 0)
      vars_cred_int = as.integer(!(plot_info$beta_cred_int[[scenario]][[i]][,1] < 0 & plot_info$beta_cred_int[[scenario]][[i]][,2] > 0))
      vars_selected[[scenario]][[i]] = vars_mpm * vars_cred_int
      
      
    }
  }
  
  return( vars_selected )
  
}

