
BrierScoreVectorised = function(beta, X, survival_data, time_points, h.g, time_int) {
  num_time_points = length(time_points)
  mcmc_iter = dim(beta)[1]
  brier_score = array(0, dim=c(mcmc_iter, num_time_points))
  #taking 1 - status.train to invert the event / censoring, to get km for censoring
  
  
  
  ordering = order(survival_data$time.train)
  surv_time = survival_data$time.train[ordering]
  surv_status = survival_data$status.train[ordering]
  
  event_times = surv_time[surv_status == 1]
  
  surv_object = Surv(time = surv_time, event = surv_status) # be sure event is surv_status or 1-surv_status=?
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
        #w_mt[t,-c(indices)] = rep(1/summary(km_censor, time=time_points[t], extend=T)$surv, n - length(indices)) #rep(1/summary(km_censor, time=time_points, extend=T)$surv, n - length(indices[2]) - 1)
        w_mt[t,-c(indices)] = 0
        
      }
      
      #browser()
      
      
    }
    
    #browser()
    
    #Ht =  rbind(Ht, H0[time_interval] %*% t(linear_pred)) ## why rbind Ht from different MCMC iterations, and add 0
    if(time_interval[1] == 0){
      Ht =  c(0, H0[time_interval]) %*% t(linear_pred)
    } else {
      Ht =  H0[time_interval] %*% t(linear_pred)
    }
    
    
    #calculating S(t)
    surv_prob = exp(-Ht) #why dim = c(50,100) instead of (51, 100)? there are 51 time points
    
    
    
    #browser()
    
    
    
    add = w_mt * ((1-indicator) - surv_prob)^2
    if(any(is.nan(add))) browser()
    isNaN = which(is.nan(add), arr.ind=T)
    add = replace(add, is.nan(add), 0)
    #brier_score[i,] = brier_score[i,] + rowSums(add)/n ## why accumulate 'brier_score[i,]', i is for MCMC iteration?
    brier_score[i,] = rowSums(add)/n 
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
