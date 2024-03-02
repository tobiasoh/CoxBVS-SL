#' BayesSurv
#' @title Fit the Bayesian Cox Model
#'
#' @description
#' This is the main function to fit a Bayesian Cox model with graph-structured 
#' selection priors for sparse identification of high-dimensional covariates.
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
#' \code{t}, \code{di}, \code{X}. See details for more information
#' @param priorPara a list containing prior parameter values
#' @param initial a list containing prior parameters' initial values
#' @param nIter the number of iterations of the chain
#' @param seed random seed
#'
#'
#' @return An object of class \code{BayesSurv} is saved as
#' \code{obj_BayesSurv.rda} in the output file, including the following components:
#' \itemize{
#' \item input - a list of all input parameters by the user
#' \item output - a list of the all output estimates:
#' \itemize{
#' \item "\code{gamma.p}" - a matrix with MCMC intermediate estimates of the indicator variables of regression coefficients.
#' \item "\code{beta.p}" - a matrix with MCMC intermediate estimates of the regression coefficients.
#' \item "\code{h.p}" - a matrix with MCMC intermediate estimates of the increments in the cumulative baseline hazard in each interval.
#' }
#' \item call - the matched call.
#' }
#'
#'
#' @examples
#'
#' # Load the example dataset
#'
#' @export
BayesSurv <- function(survObj,
                      priorPara,
                      initial,
                      nIter = 1,
                      seed = 123) {
  # survObj <- list(
  #   "t" = data$time.train, "c" = data$status.train, "X" = data$X.train,
  #   "SSig" = t(data$X.train) %*% data$X.train
  # ) # sample covariance matrix
  survObj$n <- length(survObj$t)
  survObj$p <- ncol(survObj$X)

  model.type <- "Pooled"
  # time.train = sim_surv_data[[2]] #difference between this and previously defined time.train??
  # status.train = sim_surv_data[[1]]

  # log.like  = coxph( Surv(time.train, status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariates
  # initial$log.like.ini = log.like

  # check the formula
  cl <- match.call()

  # defining parameters for the main MCMC simulation function
  nu0 <- 0.05
  nu1 <- 5
  lambda <- 3
  S <- 1


  # estimate shape and scale parameter of Weibull distribution
  # (hyperparameters for prior of H* (mean function in Gamma process prior for baseline hazard))
  #surv_obj <- Surv(time = pmin(survObj$t, survObj$di), event = as.numeric(survObj$t <= survObj$di))
  surv_obj <- Surv(survObj$t, survObj$di)
  fit <- survreg(surv_obj ~ 1, dist = "weibull", x = TRUE, y = TRUE)
  kappa0 <- 1 / exp(fit$icoef["Log(scale)"])
  eta0 <- exp(fit$coefficients)^(-kappa0)
  priorPara$kappa0 <- kappa0
  priorPara$eta0 <- eta0

  ret <- list(input = list(), output = list(), call = cl)
  class(ret) <- "BayesSurv"

  # Copy the inputs
  ret$input["p"] <- survObj$p
  ret$input["nIter"] <- nIter
  ret$input["burnin"] <- 1
  ret$input["thin"] <- 1
  ret$input["rw"] <- NULL
  ret$input$priorPara <- priorPara

  ret$output <- func_MCMC(
    survObj = survObj,
    priorPara = priorPara,
    initial = initial,
    nIter = nIter,
    S = S,
    method = "Pooled", # "Pooled" , #"Pooled", "CoxBVSSL", "Sub-struct"
    MRF_2b = FALSE,
    seed = seed
  )

  return(ret)
}
