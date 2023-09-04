source("func_MCMC.R")
source("func_MCMC_cox.R")
source("func_MCMC_graph.R")

source("Data_Simulation.R")



#prepare data
beta1.p = c(1,1,1,-1,-1,-1,0,0,0, 1)
beta2.p = c(1,1,1,-1,-1,-1,0,0,0,1)
surv.e = Surv.e
surv.c = Surv.c

n = 1000
p = 10

#simulated gene expression data, for two subgroups, split into test and training data
sim_data = sim_data_fun(n=n, p=p, surv.e=surv.e, surv.c=surv.c, beta1.p = beta1.p, beta2.p=beta2.p )


beta1 = beta1.p


sim_surv_data = sim.surv(X = sim_data$Traindata[[2]]$X.train, beta = beta1, surv.e=Surv.e[[1]], surv.c = Surv.c[[1]], n=n)

surv_obj = Surv(sim_surv_data[[2]], sim_surv_data[[1]])


X.train = sim_data$Traindata[[1]]$X.train
time.train = sim_data$Traindata[[1]]$time.train
status.train = sim_data$Traindata[[1]]$status.train



survObj   = list("t" = time.train, "c" = status.train, "X" = X.train,
                 "SSig" = t(X.train)%*%X.train ) # sample covariance matrix
survObj$n = length(survObj$t) 
survObj$p = ncol(survObj$X)

model.type = "Pooled"
time.train = sim_surv_data[[2]]
status.train = sim_surv_data[[1]]


#defining parameters for the main MCMC simulation function
nu0 = 0.05
nu1 = 5
lambda = 3
S = 1


# estimate shape and scale parameter of Weibull distribution
# (hyperparameters for prior of H* (mean function in Gamma process prior for baseline hazard))
fit     = survreg(surv_obj ~ 1, dist = "weibull", x = TRUE, y = TRUE)
kappa0  = 1/exp(fit$icoef["Log(scale)"]) 
eta0    = exp(fit$coefficients)^(-kappa0)  





priorPara = list("eta0"   = eta0,                   # prior of baseline hazard
                 "kappa0" = kappa0,                 # prior of baseline hazard
                 "c0"     = 2,                      # prior of baseline hazard
                 "tau"    = 0.0375,                 # standard deviation for prior of regression coefficients
                 "cb"     = 20,                     # standard deviation for prior of regression coefficients
                 "pi.ga"  = 0.5, #ga.pi,                  # prior variable selection probability for standard Cox models
                 "a"      = -4, #a0,                     # hyperparameter in MRF prior
                 "b"      = 0, #b0,                     # hyperparameter in MRF prior
                 "pi.G"   = 2/(p-1),                # prior edge selection probability for all edges
                 "V0"     = nu0^2 * matrix(1,p,p),  # matrix with small variances for prior of precision matrices
                 "V1"     = nu1^2 * matrix(1,p,p),  # matrix with large variances for prior of precision matrices
                 "lambda" = lambda                  # for prior of diagonal elements in precision matrices
) 


ga.ini = 0
# set starting values of parameters to be updated
if( 0 < ga.ini & ga.ini < 1 ){
  Gamma.ini = rep(0, p)
  Gamma.ini[ sample(1:p, round(p*ga.ini)) ] = 1
}else{
  Gamma.ini = rep(ga.ini, p)
}


Beta.ini = numeric(p)
Beta.ini[ which(Gamma.ini == 0) ] = runif( sum(Gamma.ini == 0), -0.02, 0.02) 
Beta.ini[ which(Gamma.ini == 1) ] = runif( sum(Gamma.ini == 1), -2, 2)

if(model.type == "Pooled"){
  beta.ini  = Beta.ini 
  gamma.ini = Gamma.ini
  log.like  = coxph( Surv(time.train, status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariables

}
  
  
  
  

G_rs = diag(1, p, p) # subgraph between subgroups (since all subgroups have the same starting values for gamma)
G_ss                  = matrix(0, p, p) # subgraph within each subgroup 
#updiag                = sample( which(upper.tri(G_ss)), round(((p^2-p)/2)*g.ini) )   
#G_ss[updiag]          = 1
G_ss[lower.tri(G_ss)] = t(G_ss)[lower.tri(G_ss)]

G.ini = do.call("cbind", rep(list(do.call("rbind", rep(list(G_rs), S))), S) )
for(g in 1:S){
  G.ini[(g-1)*p+(1:p),(g-1)*p+(1:p)] = G_ss
}

# variances for prior of precision matrices
v_ini                     = priorPara$V1
v_ini[ which(G_ss == 0) ] = nu0^2
V.ini                     = rep( list(v_ini), S ) 


C.ini = Sig.ini = rep( list(diag(1, p, p)), S) 


initial = list("G.ini" = G.ini, "V.ini" = V.ini, "C.ini" = C.ini,  "Sig.ini" = Sig.ini, 
               "gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = log.like)


mcmc_iterations = 100
#run main MCMC simulation:
fit = func_MCMC(survObj = survObj,
          priorPara = priorPara, 
          initial = initial, 
          num.reps = mcmc_iterations, 
          S = 1, 
          method = "Pooled" , #"Pooled", "CoxBVSSL", "Sub-struct" 
          MRF_2b = FALSE, 
          seed=2036237
          )



avg_gamma1 = rep(0, p)
avg_beta.p = rep(0,p)
for (i in 1:p) {
  avg_gamma1[i] = mean( fit$gamma.p[2:(mcmc_iterations+1),i] )
  avg_beta.p[i] = mean( fit$beta.p[2:(mcmc_iterations+1),i] )
  
}


