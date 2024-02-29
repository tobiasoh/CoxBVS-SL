#' BayesSurv
#' @title Update coefficients of the Bayesian Cox Lasso Model (information below to be updated)
#'
#' @description
#' This a 
#'
#' @name UpdateRPlee11
#'
#' @param sobj object of ...
#'
#' @return An object of ...
#'
#' @examples
#'
#' # Load the example dataset
#' 
#'
#' @export
UpdateRPlee11 <- function(sobj, priorPara, ini, S, method)
{
  p   = sobj$p
  tau = priorPara$tau
  cb  = priorPara$cb
  
  if(method == "Pooled"){
    n       = sobj$n
    x       = sobj$X
    J       = priorPara$J
    ind.r   = priorPara$ind.r
    ind.d   = priorPara$ind.d
    ind.r_d = priorPara$ind.r_d	  
    be.ini  = ini$beta.ini
    ga.ini  = ini$gamma.ini
    h       = ini$h
    
    #erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
    erg = updateRP_genomic_cpp(p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
    
    beta.ini  = erg$be.ini
    acceptlee = erg$acceptl
    
  }else{
    beta.ini = acceptlee = vector("list", S )  
    for(g in 1:S){  # loop through subgroups
      
      n       = sobj$n[[g]]
      x       = sobj$X[[g]]
      J       = priorPara$J[[g]]
      ind.r   = priorPara$ind.r[[g]]
      ind.d   = priorPara$ind.d[[g]]
      ind.r_d = priorPara$ind.r_d[[g]]		
      be.ini  = ini$beta.ini[[g]]
      ga.ini  = ini$gamma.ini[[g]]
      h       = ini$h[[g]]
      
      #erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
      erg = updateRP_genomic_cpp(p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
      
      beta.ini[[g]] = erg$be.ini
      acceptlee[[g]] = erg$acceptl
    }
  }
  return(list(beta.ini = beta.ini, acceptlee = acceptlee)) 
} 
# the end of "UpdateRP.lee11" function
