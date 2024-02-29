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
#' # Load the example dataset
#' data("exampleData", package = "BayesSurv")
#' p <- exampleData$p
#' q <- exampleData$q
#' survObj <- exampleData[1:3]
#'
#' # Set hyperparameters
#' mypriorPara <- list(
#'   "groupInd" = 1:p, "eta0" = 0.02, "kappa0" = 1, "c0" = 2, "r" = 10 / 9,
#'   "delta" = 1e-05, "lambdaSq" = 1, "sigmaSq" = runif(1, 0.1, 10),
#'   "beta.prop.var" = 1, "beta.clin.var" = 1
#' )
#'
#' \donttest{
#' # run Bayesian Lasso Cox
#' library("BayesSurv")
#' set.seed(123)
#' fitBayesCox <- BayesSurv(survObj,
#'   p = p, q = q, hyperpar = mypriorPara,
#'   nIter = 10, burnin = 0, outFilePath = tempdir()
#' )
#' plot(fitBayesCox, color = "blue")
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
