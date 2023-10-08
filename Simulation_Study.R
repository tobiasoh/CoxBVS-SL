source("func_MCMC.R")
source("func_MCMC_cox.R")
#source("func_MCMC_graph.R")

source("Data_Simulation.R")

make_dataset = function(n,
                        truePara)
{
  p = length(truePara$beta)
  #simulated gene expression data, for two subgroups, split into test and training data
  sim_data = sim_data_fun(n=n, p=p, surv.e=Surv.e, surv.c=Surv.c, beta1.p = truePara$beta, beta2.p = truePara$beta, cov_matrix=truePara$sigma)

  X.train = sim_data$Traindata[[1]]$X.train
  time.train = sim_data$Traindata[[1]]$time.train
  status.train = sim_data$Traindata[[1]]$status.train
  
  
  
  #for X.train and X.test, should I use scaled or unscaled
  X.train = sim_data$Traindata[[1]]$X.train
  time.train = sim_data$Traindata[[1]]$time.train
  status.train = sim_data$Traindata[[1]]$status.train
  
  X.test = sim_data$Testdata[[1]]$X.test
  time.test = sim_data$Testdata[[1]]$time.test
  status.test = sim_data$Testdata[[1]]$status.test
  
  dataset = list("X.train" = X.train, 
                 "time.train" = time.train,
                 "status.train" = status.train,
                 "X.test" = X.test,
                 "time.test" = time.test,
                 "status.test" = status.test,
                 "n" = n,
                 "truePara" = truePara)
  
  return(dataset)
  
}

simulate = function(data,
                    priorPara,
                    initial,
                    num.reps,
                    seed)
{
  
  survObj = list("t" = data$time.train, "c" = data$status.train, "X" = data$X.train,
                 "SSig" = t(data$X.train)%*%data$X.train ) # sample covariance matrix
  survObj$n = length(survObj$t) 
  survObj$p = ncol(survObj$X)
  
  model.type = "Pooled"
  #time.train = sim_surv_data[[2]] #difference between this and previously defined time.train??
  #status.train = sim_surv_data[[1]] 
  
  #log.like  = coxph( Surv(time.train, status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariates
  #initial$log.like.ini = log.like
  
  
  
  #defining parameters for the main MCMC simulation function
  nu0 = 0.05
  nu1 = 5
  lambda = 3
  S = 1
  
  
  # estimate shape and scale parameter of Weibull distribution
  # (hyperparameters for prior of H* (mean function in Gamma process prior for baseline hazard))
  surv_obj = Surv(time=pmin(survObj$t, survObj$c), event = as.numeric(survObj$t <= survObj$c))
  surv_obj = Surv(data$time.train, data$status.train)
  fit     = survreg(surv_obj ~ 1, dist = "weibull", x = TRUE, y = TRUE)
  kappa0  = 1/exp(fit$icoef["Log(scale)"]) 
  eta0    = exp(fit$coefficients)^(-kappa0)  
  priorPara$kappa0 = kappa0
  priorPara$eta0 = eta0
  
  fit = func_MCMC(survObj = survObj,
                  priorPara = priorPara, 
                  initial = initial, 
                  num.reps = num.reps, 
                  S = 1, 
                  method = "Pooled",#"Pooled" , #"Pooled", "CoxBVSSL", "Sub-struct" 
                  MRF_2b = FALSE, 
                  seed=2036237
  )
  
  
  
  return(fit) 
}

run_simulation = function(n, #sample size
                          truePara, #true parameters used to simulate the data
                          priorPara, #prior parameters
                          initial, #initial values of the parameters
                          num.reps=100,
                          seed)
{
  p = length(truePara$beta)
  #simulated gene expression data, for two subgroups, split into test and training data
  sim_data = sim_data_fun(n=n, p=p, surv.e=Surv.e, surv.c=Surv.c, beta1.p = truePara$beta, beta2.p = truePara$beta, cov_matrix=truePara$sigma)
  
  
  #simulated survival data, based on the simulated gene expression data from above
  #sim_surv_data = sim.surv(X = sim_data$Traindata[[2]]$X.train, beta = truePara$beta, surv.e=Surv.e[[1]], surv.c = Surv.c[[1]], n=n)
  
  X.train = sim_data$Traindata[[1]]$X.train
  time.train = sim_data$Traindata[[1]]$time.train
  status.train = sim_data$Traindata[[1]]$status.train
  
  survObj = list("t" = time.train, "c" = status.train, "X" = X.train,
                 "SSig" = t(X.train)%*%X.train ) # sample covariance matrix
  survObj$n = length(survObj$t) 
  survObj$p = ncol(survObj$X)
  
  model.type = "Pooled"
  #time.train = sim_surv_data[[2]] #difference between this and previously defined time.train??
  #status.train = sim_surv_data[[1]] 
  
  log.like  = coxph( Surv(time.train, status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariates
  initial$log.like.ini = log.like
  
  
  
  #defining parameters for the main MCMC simulation function
  nu0 = 0.05
  nu1 = 5
  lambda = 3
  S = 1
  
  
  # estimate shape and scale parameter of Weibull distribution
  # (hyperparameters for prior of H* (mean function in Gamma process prior for baseline hazard))
  fit     = survreg(surv_obj ~ 1, dist = "weibull", x = TRUE, y = TRUE)
  kappa0  = 1/exp(fit$icoef["Log(scale)"]) 
  eta0    = exp(fit$coefficients)^(-kappa0)  
  priorPara$kappa0 = kappa0
  priorPara$eta0 = eta0
  
  fit = func_MCMC(survObj = survObj,
                                   priorPara = priorPara, 
                                   initial = initial, 
                                   num.reps = num.reps, 
                                   S = 1, 
                                   method = "Pooled",#"Pooled" , #"Pooled", "CoxBVSSL", "Sub-struct" 
                                   MRF_2b = FALSE, 
                                   seed=2036237
                   )
  
  
  
  return(fit)
  
}










# 
# 
# 
# 
# 
# ga.ini = 0
# # set starting values of parameters to be updated
# if( 0 < ga.ini & ga.ini < 1 ){
#   Gamma.ini = rep(0, p)
#   Gamma.ini[ sample(1:p, round(p*ga.ini)) ] = 1
# }else{
#   Gamma.ini = rep(ga.ini, p)
# }
# 
# 
# Beta.ini = numeric(p)
# Beta.ini[ which(Gamma.ini == 0) ] = runif( sum(Gamma.ini == 0), -0.02, 0.02) 
# Beta.ini[ which(Gamma.ini == 1) ] = runif( sum(Gamma.ini == 1), -2, 2)
# 
# Beta.ini = rep(0, p)
# 
# if(model.type == "Pooled"){
#   beta.ini  = Beta.ini 
#   gamma.ini = Gamma.ini
#   log.like  = coxph( Surv(time.train, status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariables
#   
# }
# 
# 
# 
# 
# 
# G_rs = diag(1, p, p) # subgraph between subgroups (since all subgroups have the same starting values for gamma)
# G_ss                  = matrix(0, p, p) # subgraph within each subgroup 
# #updiag                = sample( which(upper.tri(G_ss)), round(((p^2-p)/2)*g.ini) )   
# #G_ss[updiag]          = 1
# G_ss[lower.tri(G_ss)] = t(G_ss)[lower.tri(G_ss)]
# 
# G.ini = do.call("cbind", rep(list(do.call("rbind", rep(list(G_rs), S))), S) )
# for(g in 1:S){
#   G.ini[(g-1)*p+(1:p),(g-1)*p+(1:p)] = G_ss
# }
# 
# # variances for prior of precision matrices
# v_ini                     = priorPara$V1
# v_ini[ which(G_ss == 0) ] = nu0^2
# V.ini                     = rep( list(v_ini), S ) 
# 
# 
# C.ini = Sig.ini = rep( list(diag(1, p, p)), S) 
# 
# 
# initial = list("G.ini" = G.ini, "V.ini" = V.ini, "C.ini" = C.ini,  "Sig.ini" = Sig.ini, 
#                "gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = log.like)
# 
# 
# mcmc_iterations = 100
# #run main MCMC simulation:
# fit = func_MCMC(survObj = survObj,
#                 priorPara = priorParaPooled, 
#                 initial = initial, 
#                 num.reps = mcmc_iterations, 
#                 S = 1, 
#                 method = "Pooled",#"Pooled" , #"Pooled", "CoxBVSSL", "Sub-struct" 
#                 MRF_2b = FALSE, 
#                 seed=2036237
# )
# 
# 
# 
# avg_gamma1 = rep(0, p)
# avg_beta.p = rep(0,p)
# avg_accept = rep(0,p)
# for (i in 1:p) {
#   avg_gamma1[i] = mean( fit$gamma.p[2:(mcmc_iterations+1),i] )
#   avg_beta.p[i] = mean( fit$beta.p[2:(mcmc_iterations+1),i] )
#   #avg_accept[i] = mean( fit$accept.RW[2:(mcmc_iterations+1), i] )
#   
# }
# 
# model_size = rep(0, mcmc_iterations)
# for (i in 1:mcmc_iterations) {
#   model_size[i] = sum(fit$gamma.p[i,])/mcmc_iterations
#   
# }
