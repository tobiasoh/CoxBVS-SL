#' BayesSurv
#' @title Function to run MCMC sampling
#'
#' @description
#' This an internal function for MCMC sampling
#'
#' @name func_MCMC
#' 
#' @import stats
#' @import utils
#'
#' @param survObj a list containing observed data from \code{n} subjects;
#' \code{t}, \code{di}, \code{X}. See details for more information
#' @param priorPara a list containing prior parameter values
#' @param initial a list containing prior parameters' initial values
#' @param nIter the number of iterations of the chain
#' @param burnin number of iterations to discard at the start of the chain.
#' Default is 0
#' @param thin thinning MCMC intermediate results to be stored
#' @param S the number of subgroups
#' @param method a method option from 
#' \code{c("Pooled", "CoxBVSSL", "Sub-struct", "Subgroup")}
#' @param MRF_2b two different b in MRF prior for subgraphs G_ss and G_rs
#' @param MRF_G logical value. \code{MRF_G = TRUE} is to fix the MRF graph which 
#' is provided in the argument \code{priorPara}, and \code{MRF_G = FALSE} is to 
#' use graphical model for leanring the MRF graph
#' @param output_graph_para allow (\code{TRUE}) or suppress (\code{FALSE}) the 
#' output for parameters 'G', 'V', 'C' and 'Sig' in the graphical model 
#' if \code{MRF.G = FALSE}
#'
#' @return A list object
#'
#'
#' @export
func_MCMC <- function(survObj, priorPara, initial, 
                      nIter, thin, burnin, 
                      S, method, MRF_2b, MRF_G, 
                      output_graph_para) {
  # prior parameters for grouped data likelihood of Cox model
  if (method != "Pooled") {
    s <- J <- intv <- vector("list", S)
    for (g in 1:S) {
      sg <- sort((survObj$t[[g]])[survObj$di[[g]] == 1])
      # s[[g]]    = unique(c(sg, 2 * max(survObj$t[[g]]) - max( (survObj$t[[g]])[-which(survObj$t[[g]] == max(survObj$t[[g]]))])))
      s[[g]] <- unique(c(sg, 2 * max(survObj$t[[g]]) - max((survObj$t[[g]])[-which(survObj$t[[g]] == max(survObj$t[[g]]))])))

      J[[g]] <- length(s[[g]])
      # intv[[g]] = setting.interval(survObj$t[[g]], survObj$di[[g]], s[[g]], J[[g]])
      intv[[g]] <- settingInterval_cpp(survObj$t[[g]], survObj$di[[g]], s[[g]], J[[g]])
    }
    priorPara$s <- s
    priorPara$J <- J
    priorPara$ind.r_d <- lapply(intv, function(x) x$ind.r_d)
    priorPara$ind.d <- lapply(intv, function(x) x$ind.d)
    priorPara$ind.r <- lapply(intv, function(x) x$ind.r)
    priorPara$d <- lapply(intv, function(x) x$d)

    H.star <- alpha0 <- vector("list", S)
    for (g in 1:S) {
      H.star.g <- alpha0.g <- numeric()
      for (j in 1:priorPara$J[[g]]) {
        H.star.g[j] <- priorPara$eta0[[g]] * (priorPara$s[[g]])[j]^priorPara$kappa0[[g]]
        alpha0.g[j] <- priorPara$c0 * H.star.g[j]
      }
      H.star[[g]] <- H.star.g
      alpha0[[g]] <- alpha0.g
    }
    priorPara$hPriorSh <- lapply(alpha0, function(x) diff(c(0, x)))

    h <- lapply(priorPara$J, function(x) rgamma(x, 1, 1))
  } else { # method = "Pooled"
    priorPara$s <- sort(survObj$t[survObj$di == 1])
    priorPara$s <- unique(c(priorPara$s, 2 * max(survObj$t) - max(survObj$t[-which(survObj$t == max(survObj$t))])))
    priorPara$J <- length(priorPara$s)

    # intv              = setting.interval(survObj$t, survObj$di, priorPara$s, priorPara$J)
    intv <- settingInterval_cpp(survObj$t, survObj$di, priorPara$s, priorPara$J)
    priorPara$ind.r_d <- intv$ind.r_d
    priorPara$ind.d <- intv$ind.d
    priorPara$ind.r <- intv$ind.r
    priorPara$d <- intv$d

    H.star <- alpha0 <- c()
    for (j in 1:priorPara$J) {
      H.star[j] <- priorPara$eta0 * priorPara$s[j]^priorPara$kappa0
      alpha0[j] <- priorPara$c0 * H.star[j]
    }
    priorPara$hPriorSh <- diff(c(0, alpha0))

    h <- rgamma(priorPara$J, 1, 1)
  }

  ini <- initial
  ini$h <- h

  # for posterior samples
  mcmcOutcome <- list()
  mcmcOutcome$gamma.p <- ini$gamma.ini
  mcmcOutcome$beta.p <- ini$beta.ini
  mcmcOutcome$h.p <- ini$h
  mcmcOutcome$log.jpost <- ini$log.like.ini <- NULL

  if (method == "Pooled") {
    mcmcOutcome$post.gamma <- NULL
    mcmcOutcome$gamma.margin <- ini$gamma.ini
    mcmcOutcome$beta.margin <- ini$beta.ini
  } else {
    mcmcOutcome$post.gamma <- vector("list", S)
    mcmcOutcome$gamma.margin <- rep(list(0), S)
    mcmcOutcome$beta.margin <- rep(list(0), S)
  }

  if (method %in% c("CoxBVSSL", "Sub-struct") && output_graph_para) {
    mcmcOutcome$G.p <- list(ini$G.ini)
    mcmcOutcome$V.p <- list(ini$V.ini)
    mcmcOutcome$C.p <- list(ini$C.ini)
    mcmcOutcome$Sig.p <- list(ini$Sig.ini)
  }

  # save prior parameters in MCMC output
  mcmcOutcome$s <- priorPara$s
  mcmcOutcome$eta0 <- priorPara$eta0
  mcmcOutcome$kappa0 <- priorPara$kappa0
  mcmcOutcome$c0 <- priorPara$c0
  mcmcOutcome$pi.ga <- priorPara$pi.ga
  mcmcOutcome$tau <- priorPara$tau
  mcmcOutcome$cb <- priorPara$cb

  if (method %in% c("CoxBVSSL", "Sub-struct")) {
    mcmcOutcome$V0 <- priorPara$V0
    mcmcOutcome$V1 <- priorPara$V1
    mcmcOutcome$lambda <- priorPara$lambda
    mcmcOutcome$a <- priorPara$a
    mcmcOutcome$b <- priorPara$b
    mcmcOutcome$pi.G <- priorPara$pi.G
  }

  if (method == "Pooled") {
    RW.accept <- numeric()
  } else {
    RW.accept <- vector("list", S)
  }

  # MCMC sampling
  
  # Initializes the progress bar  
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=")
  
  for (M in 1:nIter) {
    
    if (method %in% c("CoxBVSSL", "Sub-struct")) {
      # update graph and precision matrix
      network <- func_MCMC_graph(survObj, priorPara, ini, S, method, MRF_2b)

      Sig.ini <- ini$Sig.ini <- network$Sig.ini # precision matrix?
      C.ini <- ini$C.ini <- network$C.ini
      V.ini <- ini$V.ini <- network$V.ini # some sort of variance?
      G.ini <- ini$G.ini <- network$G.ini # graph
    }


    # update gamma (latent indicators of variable selection)
    # browser()
    sampleGam <- UpdateGamma(survObj, priorPara, ini, S, method, MRF_2b)
    gamma.ini <- ini$gamma.ini <- sampleGam$gamma.ini

    # update beta (regression parameters)
    # beta.tmp  = UpdateRP.lee11(survObj, priorPara, ini, S, method)
    beta.tmp <- UpdateRPlee11(survObj, priorPara, ini, S, method)
    beta.ini  = ini$beta.ini = beta.tmp$beta.ini

    # update increments in cumulative hazards
    #h <- ini$h <- UpdateBH(survObj, priorPara, ini, S, method)
    if (S == 1 && MRF_G) {
      h <- ini$h <- as.vector(updateBH_cpp(survObj$X, ini$beta.ini, 
                                 priorPara$J, priorPara$ind.r_d, 
                                 priorPara$hPriorSh, priorPara$d, priorPara$c0))
    } else {
      h <- updateBH_list_cpp(survObj$X, ini$beta.ini, 
                             priorPara$J, priorPara$ind.r_d, 
                             priorPara$hPriorSh, priorPara$d, priorPara$c0)
      h <- ini$h <- lapply(h, function(x) as.vector(x))
    }

    # profile joint posterior probability
    profJpost <- calJpost(survObj, priorPara, ini, S, method, MRF_2b)
    log.j <- profJpost$logjpost
    log.lh <- profJpost$loglike
    
    if (M > burnin) {
      if (S == 1  && MRF_G) {
        mcmcOutcome$gamma.margin <- mcmcOutcome$gamma.margin + gamma.ini
        mcmcOutcome$beta.margin <- mcmcOutcome$beta.margin + beta.ini
      } else {
        mcmcOutcome$gamma.margin <- Map("+", mcmcOutcome$gamma.margin, gamma.ini)
        mcmcOutcome$beta.margin <- Map("+", mcmcOutcome$beta.margin, beta.ini)
      }
    }
    
    # store thinning MCMC intermediate results
    if (M %% thin == 0) {
      if (method == "Pooled") {
        RW.accept = rbind(RW.accept, beta.tmp$acceptlee)
        # RW.accept <- rbind(RW.accept, as.vector(beta.tmp$acceptlee))
      } else {
        for (g in 1:S) {
          RW.accept[[g]] <- rbind(RW.accept[[g]], beta.tmp$acceptlee[[g]])
        }
      }
      mcmcOutcome$accept.RW <- RW.accept
      
      # store posterior samples
      if (method %in% c("CoxBVSSL", "Sub-struct")) {
        mcmcOutcome$log.jpost <- c(mcmcOutcome$log.jpost, log.j)
        mcmcOutcome$log.like <- rbind(mcmcOutcome$log.like, log.lh)
        
        if (output_graph_para) {
          mcmcOutcome$G.p <- c(mcmcOutcome$G.p, list(G.ini))
          mcmcOutcome$V.p <- c(mcmcOutcome$V.p, list(V.ini))
          mcmcOutcome$C.p <- c(mcmcOutcome$C.p, list(C.ini))
          mcmcOutcome$Sig.p <- c(mcmcOutcome$Sig.p, list(Sig.ini))
        }
      } else {
        if (method == "Subgroup") {
          mcmcOutcome$log.jpost <- rbind(mcmcOutcome$log.jpost, log.j)
          mcmcOutcome$log.like <- rbind(mcmcOutcome$log.like, log.lh)
        } else { # method == "Pooled"
          mcmcOutcome$log.jpost <- c(mcmcOutcome$log.jpost, log.j)
          mcmcOutcome$log.like <- c(mcmcOutcome$log.like, log.lh)
        }
      }
      
      
      if (method == "Pooled") {
        mcmcOutcome$gamma.p <- rbind(mcmcOutcome$gamma.p, gamma.ini, deparse.level = 0)
        mcmcOutcome$post.gamma <- rbind(mcmcOutcome$post.gamma, sampleGam$post.gamma, deparse.level = 0)
        mcmcOutcome$beta.p <- rbind(mcmcOutcome$beta.p, beta.ini, deparse.level = 0)
        mcmcOutcome$h.p <- rbind(mcmcOutcome$h.p, h, deparse.level = 0)
      } else {
        for (g in 1:S) {
          #browser()
          mcmcOutcome$gamma.p[[g]] <- rbind(mcmcOutcome$gamma.p[[g]], (gamma.ini)[[g]], deparse.level = 0)
          mcmcOutcome$post.gamma[[g]] <- rbind(mcmcOutcome$post.gamma[[g]], sampleGam$post.gamma[[g]], deparse.level = 0)
          mcmcOutcome$beta.p[[g]] <- rbind(mcmcOutcome$beta.p[[g]], (beta.ini)[[g]], deparse.level = 0)
          mcmcOutcome$h.p[[g]] <- rbind(mcmcOutcome$h.p[[g]], (h)[[g]], deparse.level = 0)
        }
      }
    }
    
    # if (M %% 1000 == 0) {
    #   print(M)
    # }
    
    # Sets the progress bar to the current state
    setTxtProgressBar(pb, M)
    
  } # the end of MCMC sampling
  close(pb) # Close the connection of progress bar
  
  if (S == 1 && MRF_G) {
    mcmcOutcome$gamma.margin <- mcmcOutcome$gamma.margin / (nIter - burnin)
    mcmcOutcome$beta.margin <- mcmcOutcome$beta.margin / (nIter - burnin)
  } else {
    mcmcOutcome$gamma.margin <- Map("/", mcmcOutcome$gamma.margin, nIter - burnin)
    mcmcOutcome$beta.margin <- Map("/", mcmcOutcome$beta.margin, nIter - burnin)
  }
  
  return(mcmcOutcome)
}
