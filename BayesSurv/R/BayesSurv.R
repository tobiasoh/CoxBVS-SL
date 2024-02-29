#' BayesSurv
#' @title Function to Fit the Bayesian Cox Lasso Model (information below to be updated)
#'
#' @description
#' This a 
#'
#' @name BayesSurv
#' @docType package
#' @useDynLib BayesSurv
#' @aliases BayesSurv-package
#' @importFrom Rcpp evalCpp
#' @importFrom xml2 as_xml_document write_xml
#' @importFrom stats rexp rgamma runif
#' @importFrom utils write.table
#'
#' @param survObj a list containing observed data from \code{n} subjects;
#' \code{t}, \code{di}, \code{x}. See details for more information
#' @param p number of covariates for variable selection
#' @param q number of mandatory covariates
#' @param hyperpar a list containing prior parameter values; among
#' \code{c('groupInd', 'beta.ini', 'eta0', 'kappa0', 'c0', 'r', 'delta',
#' 'lambdaSq', 'sigmaSq', 'tauSq', 's', 'h', 'beta.prop.var',
#' 'beta.clin.var')}. See details for more information
#' @param nIter the number of iterations of the chain
#' @param burnin number of iterations to discard at the start of the chain.
#' Default is 0
#' @param thin thinning MCMC intermediate results to be stored
#' @param rw when setting to "TRUE", the conventional random walk Metropolis
#' Hastings algorithm is used. Otherwise, the mean and the variance of the
#' proposal density is updated using the jumping rule described in
#' Lee et al. (2011)
#' @param outFilePath path to where the output files are to be written
#' @param tmpFolder the path to a temporary folder where intermediate data
#' files are stored (will be erased at the end of the chain). It is specified
#' relative to \code{outFilePath}
#'
#' @details
#' \tabular{ll}{
#' \code{t} \tab a vector of \code{n} times to the event \cr
#' \code{di} \tab a vector of \code{n} censoring indicators for the event time (1=event occurred, 0=censored) \cr
#' \code{x} \tab covariate matrix, \code{n} observations by \code{p} variables\cr
#' \code{groupInd} \tab a vector of \code{p} group indicator for each variable\cr
#' \code{beta.ini} \tab the starting values for coefficients \eqn{\beta}\cr
#' \code{eta0} \tab scale parameter of gamma process prior for the cumulative baseline hazard, \eqn{\eta_0 > 0}\cr
#' \code{kappa0} \tab shape parameter of gamma process prior for the cumulative baseline hazard, \eqn{\kappa_0 > 0}\cr
#' \code{c0} \tab the confidence parameter of gamma process prior for the cumulative baseline hazard, \eqn{c_0 > 0}\cr
#' \code{r} \tab the shape parameter of the gamma prior for \eqn{\lambda^2}\cr
#' \code{delta} \tab the rate parameter of the gamma prior for \eqn{\lambda^2}\cr
#' \code{lambdaSq} \tab the starting value for \eqn{\lambda^2}\cr
#' \code{sigmaSq} \tab the starting value for \eqn{\sigma^2}\cr
#' \code{tauSq} \tab the starting values for \eqn{\tau^2}\cr
#' \code{s} \tab the set of time partitions for specification of the cumulative baseline hazard function\cr
#' \code{h} \tab the starting values for \eqn{h}\cr
#' \code{beta.prop.var} \tab the variance of the proposal density for \eqn{\beta} in a random walk M-H sampler\cr
#' \code{beta.clin.var} \tab the starting value for the variance of \eqn{\beta}\cr
#' }
#'
#' @return An object of class \code{BayesSurv} is saved as
#' \code{obj_BayesSurv.rda} in the output file, including the following components:
#' \itemize{
#' \item input - a list of all input parameters by the user
#' \item output - a list of the all output estimates:
#' \itemize{
#' \item "\code{beta.p}" - a matrix with MCMC intermediate estimates of the regression coefficients.
#' \item "\code{h.p}" - a matrix with MCMC intermediate estimates of the increments in the cumulative baseline hazard in each interval.
#' \item "\code{tauSq.p}" - a vector MCMC intermediate estimates of the hyperparameter "tauSq".
#' \item "\code{sigmaSq.p}" - a vector MCMC intermediate estimates of the hyperparameter "sigmaSq".
#' \item "\code{lambdaSq.p}" - a vector MCMC intermediate estimates of the hyperparameter "lambdaSq".
#' \item "\code{accept.rate}" - a vector acceptance rates of individual regression coefficients.
#' }
#' \item call - the matched call.
#' }
#'
#'
#' @references Lee KH, Chakraborty S, and Sun J (2011). Bayesian Variable
#' Selection in Semiparametric Proportional Hazards Model for High Dimensional
#' Survival Data. \emph{The International Journal of Biostatistics}, 7(1):1-32.
#' @references Zucknick M, Saadati M, and Benner A (2015). Nonidentical twins:
#' Comparison of frequentist and Bayesian lasso for Cox models.
#' \emph{Biometrical Journal}, 57(6):959–81.
#'
#' @examples
#'
#' # Load the example dataset
#'
#' @export
BayesSurv <- function(data,
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
  
  # check the formula
  cl <- match.call()
  
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
  
  ret <- list(input = list(), output = list(), call = cl)
  class(ret) <- "BayesSurv"
  
  # Copy the inputs
  ret$input["p"] <- survObj$p
  ret$input["nIter"] <- num.reps
  ret$input["burnin"] <- 1
  ret$input["thin"] <- 1
  ret$input["rw"] <- NULL
  ret$input$priorPara <- priorPara
  
  ret$output = func_MCMC(survObj = survObj,
                  priorPara = priorPara, 
                  initial = initial, 
                  num.reps = num.reps, 
                  S = 1, 
                  method = "Pooled",#"Pooled" , #"Pooled", "CoxBVSSL", "Sub-struct" 
                  MRF_2b = FALSE, 
                  seed=seed
  )
  
  return(ret) 
}

