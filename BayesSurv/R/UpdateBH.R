#' BayesSurv
#' @title Update increment in cumulative baseline hazard
#'
#' @description
#' This contains subfunctions to update increment in cumulative baseline hazard
#'
#' @name UpdateGamma
#'
#' @param sobj a list containing observed data
#' @param priorPara a list containing prior parameter values
#' @param ini a list containing prior parameters' initial values
#' @param S the number of subgroups
#' @param method a method option from 
#' \code{c("Pooled", "CoxBVSSL", "Sub-struct", "Subgroup")}
#'
#' @return An object of ...
#'
#' @examples
#'
#' # Load the example dataset
#'
#' @export
UpdateBH <- function(sobj, priorPara, ini, S, method) {
  p <- sobj$p
  c0 <- priorPara$c0

  if (method == "Pooled") {
    n <- sobj$n
    X <- sobj$X
    J <- priorPara$J
    ind.r_d <- priorPara$ind.r_d
    d <- priorPara$d
    hPriorSh <- priorPara$hPriorSh
    beta.ini <- ini$beta.ini

    xbeta <- apply(X * matrix(rep(beta.ini, n), n, p, byrow = T), 1, sum)
    xbeta[xbeta > 700] <- 700
    exp.xbeta <- exp(xbeta)
    exp.xbeta.mat <- matrix(rep(exp.xbeta, J), n, J)
    h.rate <- apply(exp.xbeta.mat * ind.r_d, 2, sum) + c0

    H <- rgamma(J, shape = hPriorSh + d, rate = h.rate)
  } else {
    H <- vector("list", S)

    for (g in 1:S) { # loop through subgroups

      n <- sobj$n[[g]]
      X <- sobj$X[[g]]
      J <- priorPara$J[[g]]
      ind.r_d <- priorPara$ind.r_d[[g]]
      d <- priorPara$d[[g]]
      hPriorSh <- priorPara$hPriorSh[[g]]
      beta.ini <- ini$beta.ini[[g]]

      xbeta <- apply(X * matrix(rep(beta.ini, n), n, p, byrow = T), 1, sum)
      xbeta[xbeta > 700] <- 700
      exp.xbeta <- exp(xbeta)
      exp.xbeta.mat <- matrix(rep(exp.xbeta, J), n, J)
      h.rate <- apply(exp.xbeta.mat * ind.r_d, 2, sum) + c0

      H[[g]] <- rgamma(J, shape = hPriorSh + d, rate = h.rate)
    }
  }

  return(H)
}
