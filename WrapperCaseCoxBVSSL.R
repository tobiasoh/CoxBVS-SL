#' @title Wrapper function for inference of a specific type of Cox model based on simulated data using MCMC sampling.
#' Difference compared to 'WrapperSimCoxBVSSL.R': data is stored in 'data' object and
#' 'instance' contains only training and test indices.
#'
#' @description
#' Prepare training data (standarize covariate data for Pooled model); define hyperparameters for
#' all priors and set starting values for MCMC sampler.
#' 
#' @param job [Job]
#' Formal argument needed for algorithm function in 'batchtools' registry. 
#' Contains information about the job.
#' @param instance [any]
#' Formal argument needed for algorithm function in 'batchtools' registry.
#' Passes the return value of the evaluated 'problem' function.
#' @param model.type [character(1)]
#' Type of Cox model, can be one of c("CoxBVSSL", "Pooled, "Subgroup", "Sub-struct").
#' "CoxBVSSL": CoxBVS-SL model 
#' "Pooled": Standard Pooled model
#' "Subgroup": Standard Subgroup model
#' "Sub-struct": Sub-struct model
#'  @param ga.ini [numeric(1)]
#'  Initial values for latent variable selection indicators gamma, should be a value in [0,1]
#'  0 or 1: set all indicators to 0 or 1
#'  value in (0,1): rate of indicators randomly set to 1, the remaining indicators are 0
#'  @param g.ini [numeric(1)]
#'  Initial values for latent edge inclusion indicators in graph, should be a value in [0,1]
#'  0 or 1: set all random edges to 0 or 1
#'  value in (0,1): rate of indicators randomly set to 1, the remaining indicators are 0
#'  @param nu0 [numeric(1)]
#'  Hyperparameter in prior of precision matrix, default value chosen in accordance with Peterson et al. (2016)
#'  @param nu1 [numeric(1)]
#'  Hyperparameter in prior of precision matrix, default value chosen in accordance with Peterson et al. (2016)
#'  @param lambda [numeric(1)]
#'  Hyperparameter in prior of precision matrix, default value chosen in accordance with Peterson et al. (2016)
#'  @param a0 [numeric(1)]
#'  Hyperparameter in MRF prior, default value chosen based on sensitivity analysis
#'  @param b0 [numeric(1)]
#'  Hyperparameter in MRF prior, default value chosen based on sensitivity analysis
#'  @param ga.pi [numeric(1)]
#'  Hyperparameter in Bernoulli prior of variable selection indicators gamma for the standard approaches, 
#'  default value fits to default hyperparameters in MRF prior
#'  @param MRF2b [logical(1)]
#'  TRUE: two different hyperparameters b in MRF prior (values b01 and b02)
#'  FALSE: one hyperparamter b in MRF prior
#'  @param b01 [numeric(1)]
#'  Hyperparameter b in MRF prior for subgraph G_ss
#'  @param b02 [numeric(1)]
#'  Hyperparameter b in MRF prior for subgraph G_rs
#'  @param n.iter [numeric(1)]
#'  Total number of MCMC samples for MCMC chain to run
#'  @param m.seed [numeric(1)]
#'  Main seed for whole experiment; job IDs are added to this seed for each MCMC chain 
#'  
#'  
#' @return [list]
#'  Named list with MCMC output and hyperparameters of all priors.
#'  
#' "gamma.p": MCMC output of latent variable inclusion indicators.
#' List of length S (S = number of subgroups) with results for each subgroup, except for Pooled model.
#' Each list element contains a (n.iter+1)xp matrix, where rows correspond to 
#' starting value + MCMC iterations and columns correspond to covariates.
#' 
#' "beta.p": MCMC output of regression coefficients; same dimension as "gamma.p".
#' 
#' "h.p": MCMC output of the increment in the cumulative baseline hazard in each interval;
#' same dimension as "gamma.p" except for number of columns (here: different intervals).
#' 
#' "post.gamma": MCMC output of acceptance probability for gamma=1; same dimension als "gamma.p".
#' 
#' "G.p": MCMC output of graph (missing for Pooled and Subgroup model). 
#' List of length (n.iter+1) with starting value and MCMC iterations. 
#' Each list element is a (pS x pS) adjacency matrix of the graph.
#' 
#' "V.p": MCMC output for variances in prior of precision matrix (missing for Pooled and Subgroup model).
#' List of length (n.iter+1) with starting value and MCMC iterations. 
#' Each list element is a list of length S for the S subgroups, with a (p x p)-matrix containing the variances.
#'      
#' "C.p": MCMC output of the precision matrix (missing for Pooled and Subgroup model); same dimension as "V.p".
#'        
#' "Sig.p": MCMC output of the covariance matrix (missing for Pooled and Subgroup model); same dimension as "V.p".
#'  
#' "s": List of length S for the S subgroups, except for Pooled with unique observed event times 
#' for each interval in the grouped data likelihood.
#'        
#' "seed": Individual job seed for each MCMC chain; "m.seed" + job ID.
#' 
#' Hyperparameters:
#' "eta0": hyperparameter in Weibull distribution of H* (in cumulative baseline hazard).       
#' "kappa0": hyperparameter in Weibull distribution of H* (in cumulative baseline hazard).    
#' "c0": hyperparameter in Gamma process prior of cumulative baseline hazard.         
#' "pi.ga": hyperparameter in Bernoulli prior of latent variable inclusion indicators.         
#' "tau": hyperparameter (standard deviation) in prior of regression coefficients.
#' "cb": hyperparameter (standard deviation) in prior of regression coefficients.        
#' "V0": hyperparameter in prior of precision matrix (small variance); only for CoxBVSSL and Sub-struct model.      
#' "V1": hyperparameter in prior of precision matrix (large variance); only for CoxBVSSL and Sub-struct model.          
#' "lambda": hyperparameter in prior of precision matrix; only for CoxBVSSL and Sub-struct model.       
#' "a": hyperparameter in MRF prior; only for CoxBVSSL and Sub-struct model.  
#' "b": hyperparameter in MRF prior; only for CoxBVSSL and Sub-struct model.
#' "pi.G": hyperparameter in Bernoulli prior of edge inclusion indicators; only for CoxBVSSL and Sub-struct model.     
#'       
#' "log.jpost": vector of length n.iter with MCMC output of log joint posterior distribution
#' 
#' "log.like": MCMC output of grouped data log-likelihood.
#' (n.iter x S)-matrix with log-likelihood for the S subgroups, except for Pooled model.  
#' 
#' "accept.RW": MCMC output of acceptance of proposal for regression coefficient in random walk Metropolis-Hastings algorithm.



WrapperCaseCoxBVSSL = function(job, data, instance, model.type, ga.ini = 0, g.ini = 0,   
                               nu0 = 0.6, nu1 = 360, lambda = 1, a0 = -1.75, b0 = 0.5, ga.pi = 0.2, 
                               MRF2b = FALSE, b01 = 0.5, b02 = 0.5, n.iter = 20000, m.seed, ...){  

  library(survival)
  
  source("func_MCMC.R")       # main function for MCMC sampling  
  source("func_MCMC_cox.R")   # helper functions for MCMC sampling of Cox model parameters
  source("func_MCMC_graph.R") # helper functions for MCMC sampling of structure learning (graph and precision matrix)
  
  seed = m.seed + job$job.id    # seed for MCMC chain 
  set.seed(seed)
  
  train     = instance$train      # row indices of training data set

  study.ids = split(1:nrow(data$surv), data$group) 
  # training data:
  Data      = lapply(study.ids, function(x){ids = intersect(x,train); dat = list("X.train" = data$genes[ids,], "time.train" = data$surv[ids,"time"], "status.train" = data$surv[ids,"status"]); dat})

  if(model.type == "Pooled"){ # for Pooled model combine/ pool data of all subgroups
    
    time.train   = unlist( lapply(Data, function(x) x$time.train) )
    status.train = unlist( lapply(Data, function(x) x$status.train) )
    
    X.train = do.call( "rbind", lapply(Data, function(x) x$X.train) )  # unscaled covariates
    X.train = scale( X.train, scale=TRUE ) # scale combined covariates of all subgroups (instead of scaling each subgroup separately)
   
    # estimate shape and scale parameter of Weibull distribution
    # (hyperparameters for prior of H* (mean function in Gamma process prior for baseline hazard))
    fit     = survreg(Surv(time.train, status.train, type = c('right')) ~ 1, dist = "weibull", x = TRUE, y = TRUE)
    kappa0  = 1/exp(fit$icoef["Log(scale)"]) 
    eta0    = exp(fit$coefficients)^(-kappa0)  
    
    # survival object and covariate data
    survObj   = list("t" = time.train, "c" = status.train, "X" = X.train,
                     "SSig" = t(X.train)%*%X.train ) # sample covariance matrix
    survObj$n = length(survObj$t) 
    survObj$p = ncol(survObj$X)
    
  }else{ # for Subgroup or CoxBVS-SL model keep the data of all subgroups separate
    
    for(g in 1:length(Data)){ # separate scaling of covariates
      Data[[g]]$X.train = scale(Data[[g]]$X.train, scale=TRUE)
    }
    
    # estimate shape and scale parameter of Weibull distribution
    # (hyperparameters for prior of H* (mean function in Gamma process prior for baseline hazard))
    fit     = lapply(Data, function(x){survreg(Surv(x$time.train, x$status.train, type = c('right')) ~ 1, dist = "weibull", x = TRUE, y = TRUE) } )
    kappa0  = lapply(fit, function(x){ 1/exp(x$icoef["Log(scale)"])} )
    eta0    = lapply(1:length(fit), function(g){ exp(fit[[g]]$coefficients)^(-kappa0[[g]]) } )  
    
    # survival object and covariate data
    survObj = list("t" = lapply(Data, function(x)x$time.train), 
                   "c" = lapply(Data, function(x)x$status.train), 
                   "X" = lapply(Data, function(x) x$X.train ), # each subgroup scaled separately
                   "n" = lapply(Data, function(x) length(x$time.train) ) # sample size per subgroup
    )
    survObj$SSig = lapply(survObj$X, function(x) t(x)%*%x)  # sample covariance matrix per subgroup
    survObj$p = ncol(survObj$X[[1]])

  }
  
  p = survObj$p    # same number of covariates p in all subgroups
  S = length(Data) # number of subgroups
  
  # set hyperparamters of all piors
  
  if( model.type == "Sub-struct" ){ MRF2b = TRUE; b02 = 0 }
  if( MRF2b ){ b0 = c(b01, b02) }
  
  priorPara = list("eta0"   = eta0,                   # prior of baseline hazard
                   "kappa0" = kappa0,                 # prior of baseline hazard
                   "c0"     = 2,                      # prior of baseline hazard
                   "tau"    = 0.0375,                 # standard deviation for prior of regression coefficients
                   "cb"     = 20,                     # standard deviation for prior of regression coefficients
                   "pi.ga"  = ga.pi,                  # prior variable selection probability for standard Cox models
                   "a"      = a0,                     # hyperparameter in MRF prior
                   "b"      = b0,                     # hyperparameter in MRF prior
                   "pi.G"   = 2/(p-1),                # prior edge selection probability for all edges
                   "V0"     = nu0^2 * matrix(1,p,p),  # matrix with small variances for prior of precision matrices
                   "V1"     = nu1^2 * matrix(1,p,p),  # matrix with large variances for prior of precision matrices
                   "lambda" = lambda                  # for prior of diagonal elements in precision matrices
  )
  
  # set starting values of parameters to be updated
  if( 0 < ga.ini & ga.ini < 1 ){
    Gamma.ini = rep(0, p)
    Gamma.ini[ sample(1:p, round(p*ga.ini)) ] = 1
  }else{
    Gamma.ini = rep(ga.ini, p)
  }
  
  Beta.ini = numeric(p)
  Beta.ini[ which(Gamma.ini == 0) ] = runif( sum(Gamma.ini == 0), -0.02, 0.02) 
  Beta.ini[ which(Gamma.ini == 1) ] = runif( sum(Gamma.ini == 1), -2, 2)
  
  if(model.type == "Pooled"){
    beta.ini  = Beta.ini 
    gamma.ini = Gamma.ini
    log.like  = coxph( Surv(time.train, status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariables
    
  }else{
    beta.ini  = rep( list(Beta.ini), S ) 
    gamma.ini = rep( list(Gamma.ini), S ) 
    log.like  = lapply(Data, function(x){coxph( Surv(x$time.train, x$status.train, type = c('right')) ~ 1 )$loglik} )
  }
  
  initial = list("gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = log.like)
  
  if(model.type %in% c("CoxBVSSL", "Sub-struct") ){
    
    # starting value for graph G (pS x pS matrix across all subgroups)
    if( g.ini %in% c(0,1) ){  
      G_ss       = matrix(g.ini, p, p) # subgraph within each subgroup 
      diag(G_ss) = 0
      
      G_rs       = diag(g.ini, p, p)   # subgraph between subgroups
      if( model.type == "Sub-struct" ) G_rs = matrix(0, p, p)
      
      G.ini      = do.call("cbind", rep(list(do.call("rbind", rep(list(G_rs), S))), S) ) # complete graph G
      for(g in 1:S){
        G.ini[(g-1)*p+(1:p),(g-1)*p+(1:p)] = G_ss
      }
      # variances for prior of precision matrices:
      if( g.ini == 0 ) V.ini = rep( list(priorPara$V0), S ) 
      else V.ini = rep( list(priorPara$V1), S )   
      
    }else{ # g.ini is rate of selected edges in subgraphs within each subgroup
      
      G_rs = diag(1, p, p) # subgraph between subgroups (since all subgroups have the same starting values for gamma)
      G_ss                  = matrix(0, p, p) # subgraph within each subgroup 
      updiag                = sample( which(upper.tri(G_ss)), round(((p^2-p)/2)*g.ini) )   
      G_ss[updiag]          = 1
      G_ss[lower.tri(G_ss)] = t(G_ss)[lower.tri(G_ss)]
      
      if( model.type == "Sub-struct" ) G_rs = matrix(0, p, p)
      
      G.ini = do.call("cbind", rep(list(do.call("rbind", rep(list(G_rs), S))), S) )
      for(g in 1:S){
        G.ini[(g-1)*p+(1:p),(g-1)*p+(1:p)] = G_ss
      }
      
      # variances for prior of precision matrices
      v_ini                     = priorPara$V1
      v_ini[ which(G_ss == 0) ] = nu0^2
      V.ini                     = rep( list(v_ini), S ) 
      
    }
    
    # list with precision and covariance matrix per subgroup 
    C.ini = Sig.ini = rep( list(diag(1, p, p)), S) 
    
    initial = list("G.ini" = G.ini, "V.ini" = V.ini, "C.ini" = C.ini,  "Sig.ini" = Sig.ini, 
                   "gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = log.like)
  }
  
  res = func_MCMC(survObj = survObj, priorPara = priorPara, initial = initial, 
                  num.reps = n.iter, S = S, method = model.type, MRF_2b = MRF2b, seed = seed )
  
  return(res)
}
