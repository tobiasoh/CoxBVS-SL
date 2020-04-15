#' Simulate survival and gene expression data
#' Contains two functions for data simulation

# Weibull parameters for simulation of survival data
load("Weibull_param.RData")

#' Helper function for simulation of survival outcome
#' 
#' @param X [matrix]
#' nxp matrix with simulated gene expression data
#' @param beta [vector]
#' Vector with true effects of all genes 
#' @param surv.e [list]
#' List of length 2 for two different subgroups, each with 'scale' and 'shape' parameter for 
#' Weibull distribution of event times 
#' @param surv.c [list]
#' List of length 2 for two different subgroups, each with 'scale' and 'shape' parameter for 
#' Weibull distribution of censoring times
#' @param n [integer(1)]
#' Number of observations
#' 
#' @return [list]
#' List with outcome (status and time vectors) 

sim.surv = function(X, beta, surv.e, surv.c, n){
  library(survival)
  
  # simulate event times from Weibull distribution
  dt = (-log(runif(n)) * (1/surv.e$scale) * exp(-X %*% beta))^(1/surv.e$shape)
  
  # simulate censoring times from Weibull distribution
  cens = rweibull(n, shape = surv.c$shape, scale = ((surv.c$scale)^(-1/surv.c$shape)))
  
  # observed time and status for each observation
  status = ifelse(dt <= cens, 1, 0)
  time = pmin(dt, cens)
  
  return(list(as.numeric(status), as.numeric(time)))
}


#' Generate one training and one test data set for each of the two subgroups with 
#' multivariate normally-distributed gene expression data and Weibull-distributed survival data
#' 
#' @param n [integer(1)]
#' Number of observations of each subgroup
#' @param p [integer(1)]
#' Number of genes in total (prognostic genes + noise) of each subgroup
#' @param surv.e [list]
#' List of length 2 for two different subgroups, each with 'scale' and 'shape' parameter for 
#' Weibull distribution of event times 
#' @param surv.c [list]
#' List of length 2 for two different subgroups, each with 'scale' and 'shape' parameter for 
#' Weibull distribution of censoring times
#' @param beta1.p [vector]
#' Vector with true effects of first prognostic genes (remaining genes have no effect) in
#' subgroup 1
#' @param beta2.p [vector]
#' Vector with true effects of first prognostic genes (remaining genes have no effect) in
#' subgroup 2
#' 
#' @return [list]
#' List of length 2 with training and test data for both subgroups, including unscaled and scaled
#' gene expression data, survival time and status.

sim_data_fun <- function(n, p, surv.e, surv.c, beta1.p, beta2.p){
  
  library(survival)
  library(mvtnorm) 
  library(MASS)
  
  p.e = length(beta1.p) # Number of prognostic variables 
  S=2 # Number of subgroups
  
  # True effects in each subgroup
  beta1 = c( beta1.p, rep(0,p-p.e) )
  beta2 = c( beta2.p, rep(0,p-p.e) )
  
  # Covariance matrix in both subgroups
  sigma = diag(p)
  block = matrix(rep(.5,9), nrow=3); diag(block) = 1
  sigma[1:3, 1:3] = sigma[4:6, 4:6] = sigma[7:9, 7:9] = block

  # Sample gene expression data from multivariate normal distribution
  X.train1 = mvrnorm(n, rep(0,p), sigma)  
  X.train2 = mvrnorm(n, rep(0,p), sigma)
  X.test1 = mvrnorm(n, rep(0,p), sigma)  
  X.test2 = mvrnorm(n, rep(0,p), sigma)
  
  # Simulate survival data for both subgroups
  surv1 = sim.surv(X.train1, beta1, surv.e[[1]], surv.c[[1]], n) 
  surv2 = sim.surv(X.train2, beta2, surv.e[[2]], surv.c[[2]], n) 
  surv1t = sim.surv(X.test1, beta1, surv.e[[1]], surv.c[[1]], n) 
  surv2t = sim.surv(X.test2, beta2, surv.e[[2]], surv.c[[2]], n) 
  
  # Combine training and test data of both subgroups
  data1 = list("X.train" = X.train1, "time.train" = surv1[[2]], "status.train" = surv1[[1]])
  data2 = list("X.train" = X.train2, "time.train" = surv2[[2]], "status.train" = surv2[[1]])
  Data = list(data1, data2)
  
  data1 = list("X.test" = X.test1, "time.test" = surv1t[[2]], "status.test" = surv1t[[1]])
  data2 = list("X.test" = X.test2, "time.test" = surv2t[[2]], "status.test" = surv2t[[1]])
  Data.test = list(data1, data2)
  
  # Scale covariates using parameters of training data
  sd.train = lapply(Data, function(X) apply(X$X.train,2,sd))  
  for(g in 1:length(Data)){
    Data.test[[g]]$X.test <- scale(Data.test[[g]]$X.test, scale=sd.train[[g]])
    Data[[g]]$X.train <- scale(Data[[g]]$X.train, scale=sd.train[[g]])
  }
  
  # Unscaled covariates (for Pooled model):
  Data[[1]]$X.train.unsc = X.train1
  Data[[2]]$X.train.unsc = X.train2
  Data.test[[1]]$X.test.unsc = X.test1
  Data.test[[2]]$X.test.unsc = X.test2
  
  return(list("Traindata" = Data, "Testdata" = Data.test))
}
