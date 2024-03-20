#' BayesSurv
#'
#' @title Time-dependent Brier scores (Not yet subgroups)
#'
#' @description
#' Predict time-dependent Brier scores based on Cox regression models
#'
#' @name plotBrier
#' 
#' @importFrom ggplot2 ggplot aes geom_step theme element_blank
#' 
#' @param object fitted object obtained with \code{BayesSurv}
#' @param survObj.new a list containing observed data from new subjects with
#' components \code{t}, \code{di}, \code{X}
#' @param times maximum time point to evaluate the prediction
#' @param method option to use the posterior mean (\code{"mean"}) of coefficients
#' for prediction or Bayesian model averaging (\code{"BMA"}) for prediction
#' @param subgroup index of the subgroup in \code{survObj.new} for prediction. 
#' Default value is 1 
#'
#' @keywords survival
#' @examples
#'
#' library("survival")
#' library("GGally")
#' library("BayesSurv")
#' set.seed(123)
#' 
#' # Load the example dataset
#' data("simData", package = "BayesSurv")
#' 
#' dataset = list("X" = simData[[1]]$X, 
#'                "t" = simData[[1]]$time,
#'                "di" = simData[[1]]$status)
#' 
#' # Initial value: null model without covariates
#' initial = list("gamma.ini" = rep(0, ncol(dataset$X)))
#' # Hyperparameters
#' priorParaPooled = list(
#'   "c0"     = 2,                      # prior of baseline hazard
#'   "tau"    = 0.0375,                 # sd for coefficient prior
#'   "cb"     = 20,                     # sd for coefficient prior
#'   "pi.ga"  = 0.02,                   # prior variable selection probability for standard Cox models
#'   "nu0"    = 0.05,                   # hyperparameter in graphical model
#'   "nu1"    = 5,                      # hyperparameter in graphical model
#'   "lambda" = 3,                      # hyperparameter in graphical model
#'   "a"      = -4,                     # hyperparameter in MRF prior
#'   "b"      = 0.1,                    # hyperparameter in MRF prior
#'   "G"      = simData$G               # hyperparameter in MRF prior
#' )    
#'
#' \donttest{
#' # run Bayesian Cox with graph-structured priors
#' fit = BayesSurv(Data=dataset, priorPara=priorParaPooled, 
#'                 initial=initial, nIter=100)
#' # predict survival probabilities of the train data
#' plotBrier(fit, survObj.new = dataset)
#' }
#'
#' @export
plotBrier <- function(object, survObj.new = NULL, 
                      method = "mean", times = NULL, subgroup = 1) {
  
  if (is.null(times)) 
    times <- sort(unique(survObj.new$t))
  capture.output(
    Brier_score <- predict.BayesSurv(object, 
                                     survObj.new = survObj.new, 
                                     method = method, 
                                     times = times,
                                     subgroup = subgroup),
    file = nullfile())
  
  Brier <- model <- NULL
  #Brier_score %>%
    ggplot2::ggplot(Brier_score, 
                    aes(times, Brier, group = model, color = model)) +
    xlab("Evaluation time points") +
    ylab("Brier score") +
    geom_step(direction = "vh") +
    theme(legend.title = element_blank())
}
