

integratedBrierScore2 = function(brier_score, time_points) {
  ibs = brier_score[1]*time_points[1]
  
  for (t in 2:length(time_points)) {
    interval = time_points[t] - time_points[t-1]
    ibs = ibs + brier_score[t]
  }
  ibs = ibs / time_points[length(time_points)]
  
  return(ibs)
}

integratedBrierScore = function(brier_score, time_points) {
  mcmc_iter = dim(brier_score)[1]
  ibs = rep(0, mcmc_iter)
  for (i in 1:mcmc_iter) {
  
    ibs[i] = brier_score[i,1]*time_points[1]
    
    for (t in 2:length(time_points)) {
      interval = time_points[t] - time_points[t-1]
      ibs[i] = ibs[i] + brier_score[i,t]
    }
    ibs[i] = ibs[i] / time_points[length(time_points)]
    
  }
  
  return(ibs)
}



### OTHER WAY OF AVERAGING



survProb_MCMC = function(beta, X, h.g) {
  mcmc_iter = dim(beta)[1]
  n = dim(X)[1]
  time_points = dim(h.g)[2]
  survival_prob = array( 0, dim=c(mcmc_iter, time_points, n))
  #browser()
  for (m in 1:mcmc_iter) {
    linear_pred = exp(X %*% beta[m,])
    H0 = cumsum(h.g[m,])
    for (t in 1:time_points) { #can vectorise this loop?
      #H0 = sum(h.g[m, 1:t])
      Ht = H0[t] * linear_pred
      survival_prob[m,t,] = exp(-Ht)
      
    }
  }
  
  return(survival_prob)
}

#surv_MCMC = survProb_MCMC(result$beta.p[-(1:warmup),], dataset$X.train, result$h.p[-(1:warmup),])


BrierScore2 = function(beta, X, survival_data, time_points, h.g) {
  num_time_points = length(time_points)
  mcmc_iter = dim(beta)[1]
  brier_score = array(0, dim=c(mcmc_iter, num_time_points))
  #taking 1 - status.train to invert the event / censoring, to get km for censoring
  

  
  ordering = order(survival_data$time.train)
  surv_time = survival_data$time.train[ordering]
  surv_status = survival_data$status.train[ordering]
  
  event_times = surv_time[surv_status == 1]
  
  surv_object = Surv(time = surv_time, event = 1 - surv_status)
  km_censor = survfit(surv_object ~ 1)
  n = length(survival_data$time.train) #number of patients
  
  
  
  for (i in 1:mcmc_iter) {
    
    linear_pred = exp(X %*% beta[i,])
    H0 = cumsum(h.g[i,])
    
    for (t in 1:num_time_points) {
      #which time interval the current time point is in
      time_interval = sum(event_times <= time_points[t])
      if (time_interval == 0) {
        Ht = 0
      #} else if (time_interval == length(event_times)) { #in this case, time point is larger than final event time
        #stop(sprintf("Brier score: time point %f is after the final event time %f", time_points[t], event_times[length(event_times)]))
      } else {
        Ht = H0[time_interval] * linear_pred
      }
      
      #browser()
      
      #calculating S(t)
      surv_prob = exp(-Ht)
      
      
      #
      
      
      
      
      indicator = as.numeric(surv_time <= time_points[t])
      w_mt = rep(0, length(indicator))
      indices = which(indicator == 1)
      #browser()
      if (identical(indices, integer(0))) {
        w_mt = rep(1/summary(km_censor, time=time_points[t], extend=T)$surv, length(w_mt))
      } else {
        w_mt[indices] = (surv_status / summary(km_censor, time=surv_time, extend=T)$surv)[indices]
        w_mt[-indices] = rep(1/summary(km_censor, time=time_points[t], extend=T)$surv, length(w_mt) - length(indices))
        
      }
  
      
      
      add = w_mt * (indicator - surv_prob)^2
      isNaN = which(is.nan(add))
      #browser()
      if (identical(isNaN, integer(0))) {
        brier_score[i,t] = brier_score[i,t] + sum(add)/n
      } else {
        brier_score[i,t] = brier_score[i,t] + sum(add[-isNaN])/n
      }
      
        

      
    }
    #browser()
  }
  
  return(brier_score)
}


#brier2 = BrierScore2(result$beta[-(1:25000),], dataset$X.train, survival_data, result$s, result$h.p[-(1:25000),])








BrierScoreVectorised2 = function(beta, X, survival_data, time_points, h.g) {
  num_time_points = length(time_points)
  mcmc_iter = dim(beta)[1]
  brier_score = array(0, dim=c(mcmc_iter, num_time_points))
  #taking 1 - status.train to invert the event / censoring, to get km for censoring
  
  
  
  ordering = order(survival_data$time.train)
  surv_time = survival_data$time.train[ordering]
  surv_status = survival_data$status.train[ordering]
  
  event_times = surv_time[surv_status == 1]
  
  surv_object = Surv(time = surv_time, event = 1 - surv_status)
  km_censor = survfit(surv_object ~ 1)
  n = length(survival_data$time.train) #number of patients
  
  
  
  for (i in 1:mcmc_iter) {
    
    linear_pred = exp(X %*% beta[i,])
    H0 = cumsum(h.g[i,])
    
    time_interval = rep(0, num_time_points)
    indicator = array(0, dim=c(num_time_points, n))
    w_mt = array(0, dim=dim(indicator))
    Ht = array(dim=c(0,n))
    
    for (t in 1:num_time_points) {
      #which time interval the current time point is in
      time_interval[t] = sum(event_times <= time_points[t])
      if (time_interval[t] == 0) {
        Ht = rbind(Ht, rep(0, 100))
      } #else if (time_interval == length(event_times)) { #in this case, time point is larger than final event time
        #stop(sprintf("Brier score: time point %f is after the final event time %f", time_points[t], event_times[length(event_times)]))
      #} else {
      #  Ht[t] = H0[time_interval] * linear_pred
      #}
      
      indicator[t,] = as.numeric(surv_time <= time_points[t])
      
      
    }
    
    #browser()
    
    
    Ht =  rbind(Ht, H0[time_interval] %*% t(linear_pred))
    
      
    #calculating S(t)
    surv_prob = exp(-Ht) #why dim = c(50,100) instead of (51, 100)? there are 51 time points

    indices = which(indicator == 1, arr.ind=T)
    #browser()
    if (identical(indices, integer(0))) {
      w_mt = array(1/summary(km_censor, time=time_points, extend=T)$surv, dim=c(num_time_points,n))
    } else {
      w_mt[indices] = 1/summary(km_censor, time=surv_time[survival_data$status.train == 1], extend=T)$surv#(surv_status / summary(km_censor, time=surv_time[survival_data$status.train == 1], extend=T)$surv)#[indices[,1]]
      w_mt[-c(indices)] = 1/summary(km_censor, time=time_points, extend=T)$surv #rep(1/summary(km_censor, time=time_points, extend=T)$surv, n - length(indices[2]) - 1)
      
    }
    
    #browser()
    
    
    add = w_mt * (indicator - surv_prob)^2
    isNaN = which(is.nan(add))
    #browser()
    if (identical(isNaN, integer(0))) {
      brier_score[i,] = brier_score[i,] + rowSums(add)/n
    } else {
      brier_score[i,] = brier_score[i,] + rowSums(add[-isNaN])/n
    }
   browser()
    
    
    
    
  }
  #browser()
return(brier_score)
}
  



BrierScoreVectorised = function(beta, X, survival_data, time_points, h.g, time_int) {
  num_time_points = length(time_points)
  mcmc_iter = dim(beta)[1]
  brier_score = array(0, dim=c(mcmc_iter, num_time_points))
  #taking 1 - status.train to invert the event / censoring, to get km for censoring
  
  
  
  ordering = order(survival_data$time.train)
  surv_time = survival_data$time.train[ordering]
  surv_status = survival_data$status.train[ordering]
  
  event_times = surv_time[surv_status == 1]
  
  surv_object = Surv(time = surv_time, event = 1 - surv_status)
  km_censor = survfit(surv_object ~ 1)
  n = length(survival_data$time.train) #number of patients
  
  
  
  for (i in 1:mcmc_iter) {
    
    linear_pred = exp(X %*% beta[i,])
    H0 = cumsum(h.g[i,])
    
    time_interval = rep(0, num_time_points)
    indicator = array(0, dim=c(num_time_points, n))
    w_mt = array(0, dim=dim(indicator))
    Ht = array(dim=c(0,n))
    
    for (t in 1:num_time_points) {
      #which time interval the current time point is in
      #browser()
      time_interval[t] = sum(time_points[t] >= time_int)
      #if (is.nan(time_interval[t])) {time_interval[t] = length(time_int)}
      if (identical(time_interval[t], numeric(0))) {time_interval[t] = 0} 
      
      #time_interval[t] time_int[sum()]

      if (time_interval[t] == 0) {
        Ht = rbind(Ht, rep(0, 100))
      } #else if (time_interval == length(event_times)) { #in this case, time point is larger than final event time
      #stop(sprintf("Brier score: time point %f is after the final event time %f", time_points[t], event_times[length(event_times)]))
      #} else {
      #  Ht[t] = H0[time_interval] * linear_pred
      #}
      
      indicator[t,] = as.numeric(surv_time <= time_points[t])
      
      indices = which(indicator[t,] == 1)
      
      #browser()
      if (identical(indices, integer(0))) {
        w_mt[t,] = rep(1/summary(km_censor, time=time_points[t], extend=T)$surv, n)
      } else {
        w_mt[t,indices] = (surv_status/summary(km_censor, time=surv_time, extend=T)$surv)[indices]#(surv_status / summary(km_censor, time=surv_time[survival_data$status.train == 1], extend=T)$surv)#[indices[,1]]
        w_mt[t,-c(indices)] = rep(1/summary(km_censor, time=time_points[t], extend=T)$surv, n - length(indices)) #rep(1/summary(km_censor, time=time_points, extend=T)$surv, n - length(indices[2]) - 1)
        
      }
      
      #browser()
      
      
    }
    
    #browser()
    
    
    Ht =  rbind(Ht, H0[time_interval] %*% t(linear_pred))
    
    
    #calculating S(t)
    surv_prob = exp(-Ht) #why dim = c(50,100) instead of (51, 100)? there are 51 time points
    
    
    
    #browser()
    
    
    
    add = w_mt * ((1-indicator) - surv_prob)^2
    isNaN = which(is.nan(add), arr.ind=T)
    add = replace(add, is.nan(add), 0)
    brier_score[i,] = brier_score[i,] + rowSums(add)/n
    #browser()
    #if (identical(isNaN, integer(0))) {
    #  brier_score[i,] = brier_score[i,] + rowSums(add)/n
    #} else {
    #  brier_score[i,] = brier_score[i,] + rowSums(add[-c(isNaN)])/n
    #}
    #browser()
    
    
    
    
  }
  #browser()
  return(brier_score)
}


brier_score_mpm = function(beta, X, survival_data, h.g) {
  num_time_points = length(h.g)
  brier_score = array(0, dim=c(mcmc_iter, num_time_points))
  #taking 1 - status.train to invert the event / censoring, to get km for censoring
  
  
  
  ordering = order(survival_data$time.train)
  surv_time = survival_data$time.train[ordering]
  surv_status = survival_data$status.train[ordering]
  
  event_times = surv_time[surv_status == 1]
  
  surv_object = Surv(time = surv_time, event = 1 - surv_status)
  km_censor = survfit(surv_object ~ 1)
  n = length(survival_data$time.train) #number of patients
  
  
  

  
  linear_pred = exp(X %*% beta[i,])
  H0 = cumsum(h.g[i,])
  
  time_interval = rep(0, num_time_points)
  indicator = array(0, dim=c(num_time_points, n))
  w_mt = array(0, dim=dim(indicator))
  Ht = array(dim=c(0,n))
  
  for (t in 1:num_time_points) {
    #which time interval the current time point is in
    #browser()
    time_interval[t] = sum(time_points[t] >= time_int)
    #if (is.nan(time_interval[t])) {time_interval[t] = length(time_int)}
    if (identical(time_interval[t], numeric(0))) {time_interval[t] = 0} 
    
    #time_interval[t] time_int[sum()]
    
    if (time_interval[t] == 0) {
      Ht = rbind(Ht, rep(0, 100))
    } #else if (time_interval == length(event_times)) { #in this case, time point is larger than final event time
    #stop(sprintf("Brier score: time point %f is after the final event time %f", time_points[t], event_times[length(event_times)]))
    #} else {
    #  Ht[t] = H0[time_interval] * linear_pred
    #}
    
    indicator[t,] = as.numeric(surv_time <= time_points[t])
    
    indices = which(indicator[t,] == 1)
    
    #browser()
    if (identical(indices, integer(0))) {
      w_mt[t,] = rep(1/summary(km_censor, time=time_points[t], extend=T)$surv, n)
    } else {
      w_mt[t,indices] = (surv_status/summary(km_censor, time=surv_time, extend=T)$surv)[indices]#(surv_status / summary(km_censor, time=surv_time[survival_data$status.train == 1], extend=T)$surv)#[indices[,1]]
      w_mt[t,-c(indices)] = rep(1/summary(km_censor, time=time_points[t], extend=T)$surv, n - length(indices)) #rep(1/summary(km_censor, time=time_points, extend=T)$surv, n - length(indices[2]) - 1)
      
    }
    
    #browser()
    
    
  }
  
  #browser()
  
  
  Ht =  rbind(Ht, H0[time_interval] %*% t(linear_pred))
  
  
  #calculating S(t)
  surv_prob = exp(-Ht) #why dim = c(50,100) instead of (51, 100)? there are 51 time points
  
  
  
  #browser()
  
  
  
  add = w_mt * ((1-indicator) - surv_prob)^2
  isNaN = which(is.nan(add), arr.ind=T)
  add = replace(add, is.nan(add), 0)
  brier_score[i,] = brier_score[i,] + rowSums(add)/n
  #browser()
  #if (identical(isNaN, integer(0))) {
  #  brier_score[i,] = brier_score[i,] + rowSums(add)/n
  #} else {
  #  brier_score[i,] = brier_score[i,] + rowSums(add[-c(isNaN)])/n
  #}
  #browser()
  
  
    
    
  #browser()
  return(brier_score)
}
}





#brier2 = BrierScoreVectorised(result$beta[-(1:29000),], dataset$X.train, survival_data, seq(5,8,.1), result$h.p[-(1:29000),], result$s)
#brier_old = BrierScore2(result$beta[-(1:29000),], dataset$X.train, survival_data, seq(0,5,.1), result$h.p[-(1:29000),])

