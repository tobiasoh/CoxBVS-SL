
library(survival)
load("SimStudy/full_sim_thin6/empty1.RData")
load("SimStudy/full_sim_thin6/true1.RData")

# check paths of MCMC estimated baseline hazards
pdf("h_p.pdf", width = 8, height = 5)
matplot(simulation_result$result$h.p, type='l', xlab='MCMC iteration', ylab='h.p')
dev.off()

source("brier_calculations_zz.R")
warmup = ceiling( simulation_result$mcmcIterations / (2*simulation_result$thinning ) )

# It does not make sense to specify a new time partition, 
# because the increments of baseline hazards h.p were sampled (in MCMC) only based on the partition of observed events
#time_star = seq(0, 5, .1)
time_star = unique(sort(simulation_result$survival_data$time.train))
brier = BrierScoreVectorised( simulation_result$result$beta[-(1:warmup),], 
                              simulation_result$X.train, 
                              simulation_result$survival_data, 
                              time_star, 
                              simulation_result$result$h.p[-(1:warmup),], 
                              simulation_result$result$s )
layout(matrix(1:3, nrow = 1))
plot(colMeans(brier) ~ time_star, type = "S", lty = 1, ylim = c(0, max(brier)), 
     xlab = "time", ylab = "Brier score", main = "Tobias - modified")


#==============================================
## Compute Brier score based on 'riskRegression' package
## Note that the baseline hazard function (h.p) as estimated from the training data is not used for predictions, similar to Zucknick et al. (2015)
#==============================================
beta_mean <- apply(simulation_result$result$beta[-(1:warmup),], 2, mean)
lp <- as.vector(simulation_result$X.train %*% beta_mean)
times <- simulation_result$survival_data$time.train
status <- simulation_result$survival_data$status.train
data_train <- data.frame(time = times, status = status, lp = lp)
bayes_train <- coxph(Surv(time, status) ~ lp, 
                     data = data_train, y = TRUE, x = TRUE)
library(riskRegression)
Brier_train <- riskRegression::Score(list("Brier_train" = bayes_train), 
                                     formula = Surv(time, status) ~ 1, 
                                     data = data_train, conf.int = FALSE, 
                                     metrics = "brier", summary = "ibs", 
                                     times = sort(unique(data_train$time)))$Brier$score
                                     #times = time_star)$Brier$score
Brier_score <- Brier_train[Brier_train$model != "Null model", ]
plot(Brier_score$Brier ~ Brier_score$times, type = "S", lty = 1, 
     xlab = "time", ylab = "Brier score", main = "riskRegression")

(ibs_riskRegression <- Brier_score$IBS[which.max(Brier_score$times)])

#==============================================
## Compute Brier scores based on 'pec' package
#==============================================
library(pec)
coxBVS <- function(mcmcOutcome, ...){
  out <- list("mcmcOutcome"=mcmcOutcome)
  class(out) <- "coxBVS"
  out$call <- match.call()
  out
}
predictSurvProb.coxBVS <- function (object, newdata, times, ...) 
  #object = object of class coxBVS
  #newdata = new data set
  #times = times at which to predict survival probabilities
{
  x <- as.matrix(newdata[,-match(c("time","status"),colnames(newdata))])
  
      #use marginal posterior means:
      beta <- object$mcmcOutcome$beta.p #remove starting solution and burnin phase
      beta.est <- apply(beta, 2, mean)
  lp <- x %*% beta.est
  
  require(survival)
  survival.coxph <- getFromNamespace("coxph", ns = "survival")
  survival.survfit.coxph <- getFromNamespace("survfit.coxph", ns = "survival")
  survival.summary.survfit <- getFromNamespace("summary.survfit", ns = "survival")
  
  f <- survival.coxph(Surv(newdata$time, newdata$status) ~ lp)
  survfit.object <- survival.survfit.coxph(f, newdata=data.frame(lp=lp), se.fit=FALSE, conf.int=FALSE)
  inflated.pred <- survival.summary.survfit(survfit.object, times = times)
  p <- t(inflated.pred$surv)
  if ((miss.time <- (length(times) - NCOL(p))) > 0) 
    p <- cbind(p, matrix(rep(NA, miss.time * NROW(p)), nrow = NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop("Prediction failed")
  p
}
mcmcResults <- list()
mcmcResults$beta.p <- simulation_result$result$beta[-(1:warmup),]
coxBVS.mean.object <- coxBVS(mcmcResults)
pe <- pec(list("Bayes lasso (mean)"=coxBVS.mean.object), 
          data=data_train, formula=Surv(time,status)~1, times=max(time_star))

brier_pec <- pe$AppErr$`Bayes lasso (mean)`
plot(brier_pec ~ pe$time, type = "S", lty = 1, ylim = c(0, max(brier_pec)), 
     xlab = "time", ylab = "Brier score", main = "pec")

(ibs_pec <- crps(pe, times=max(pe$time)))
