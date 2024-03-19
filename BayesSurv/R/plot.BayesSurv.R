#' @title create a plot of estimated coefficients
#' @description
#' Plot point estimates of regression coefficients and 95\% credible intervals
#'
#' @name plot.BayesSurv
#'
#' @importFrom GGally ggcoef
#' @importFrom ggplot2 xlab ylab
#' @importFrom stats quantile
#'
#' @param x an object of class \code{BayesSurv} or a matrix. If \code{x}
#' is a matrix, use \code{BayesSurv:::plot.BayesSurv(x)}
#' @param type type of point estimates of regression coefficients. One of
#' \code{c("mean", "median")}. Default is \code{mean}
#' @param interval logical argument to show 95\% credible intervals. Default
#' is \code{TRUE}
#' @param ... additional arguments sent to \code{ggplot2::geom_point()}
#'
#' @return ggplot object
#'
#' @examples
#' 
#' library("survival")
#' library("GGally")
#' library("BayesSurv")
#' 
#' # Load the example dataset
#' data("simData", package = "BayesSurv")
#' 
#' dataset = list("X" = simData[[1]]$X, 
#'                "t" = simData[[1]]$time,
#'                "di" = simData[[1]]$status)
#' 
#' # Initial value: null model without covariates
#' log.like  = coxph( Surv(dataset$t, dataset$di, type = c('right')) ~ 1 )$loglik 
#' initial = list("gamma.ini" = rep(0, ncol(dataset$X)), 
#'                "beta.ini" = rep(0, ncol(dataset$X)), 
#'                "log.like.ini" = log.like)
#'                
#' # Prior parameters
#' priorParaPooled = list(
#'   #"eta0"   = eta0,                   # prior of baseline hazard
#'   #"kappa0" = kappa0,                 # prior of baseline hazard
#'   "c0"     = 2,                      # prior of baseline hazard
#'   "tau"    = 0.0375,                 # sd for coefficient prior
#'   "cb"     = 20,                     # sd for coefficient prior
#'   "pi.ga"  = 0.02, #0.5, ga.pi,      # prior variable selection probability for standard Cox models
#'   "nu0"    = 0.05,                   # hyperparameter in graphical model
#'   "nu1"    = 5,                      # hyperparameter in graphical model
#'   "lambda" = 3,                      # hyperparameter in graphical model
#'   "a"      = -4, #a0,                # hyperparameter in MRF prior
#'   "b"      = 0.1, #b0,               # hyperparameter in MRF prior
#'   "G"      = simData$G               # hyperparameter in MRF prior
#' )    
#'
#' \donttest{
#' # run Bayesian Cox with graph-structured priors
#' fit = BayesSurv(survObj=dataset, priorPara=priorParaPooled, 
#'                 initial=initial, nIter=100, seed=123)
#'
#' # show posterior mean of coefficients and 95% credible intervals
#' plot(fit) + 
#'   coord_flip() + 
#'     theme(axis.text.x = element_text(angle = 90, size = 7))
#' }
#'
#' @export
plot.BayesSurv <- function(x, type = "mean", interval = TRUE, ...) {
  if (!(inherits(x, "BayesSurv") || is.matrix(x))) {
    stop("Use only with 'BayesSurv' object or a matrix!")
  }

  if (length(type) == 1) {
    if (!type %in% c("mean", "median")) {
      stop("'type' should be one of c('mean', 'median')!")
    }
  } else {
    stop("'type' should be one of c('mean', 'median')!")
  }

  if (!is.logical(interval)) {
    stop("Argument 'interval' must be a logical value!")
  }

  if (inherits(x, "BayesSurv")) {
    if (is.null(colnames(x$output$beta.p))) {
      x_names <- paste0("x", seq_len(ncol(x$output$beta.p)))
    } else {
      x_names <- colnames(x$output$beta.p)
    }
    beta_p <- x$output$beta.p[-(1:(x$input$burnin / x$input$thin + 1)), ]
  } else {
    if (is.null(colnames(x))) {
      x_names <- paste0("x", seq_len(ncol(x)))
    } else {
      x_names <- colnames(x)
    }
    beta_p <- x
  }

  # pdf("psbcBeta.pdf", height = 5, width = 3.5)
  beta_est <- apply(beta_p, 2, type)
  beta_L <- apply(beta_p, 2, quantile, 0.025)
  beta_U <- apply(beta_p, 2, quantile, 0.975)
  tbl <- data.frame(term = x_names, estimate = beta_est, conf.low = beta_L, conf.high = beta_U)
  tbl$term <- factor(tbl$term, levels = tbl$term)

  # Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
  pCoef <- ggcoef(tbl, ...) + xlab(expression(Posterior ~ ~beta)) + ylab("")
  pCoef
  # dev.off()
}
