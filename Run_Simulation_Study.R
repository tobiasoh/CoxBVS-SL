library(survival)
source("/data/tobiasoh/Simulation_Study.R")


#setting the true parameters, that we use to simulate the data set
seed = sample(1:10e7, 1)
set.seed(seed)

args = commandArgs(trailingOnly = T)
graph = args[1]
mcmc_iterations = as.integer(args[2])

p = 200
trueBeta = c(1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1)
#trueBeta = c(1,1,1,1,-1,-1,-1)
trueBeta = c(trueBeta, rep(0, p - length(trueBeta)))

#sigma = diag(p)
#block = matrix(rep(.5,9), nrow=3); diag(block) = 1
#sigma[1:3, 1:3] = block#sigma[4:6, 4:6] = sigma[7:9, 7:9] = block


sigma = diag(p)
block = matrix(rep(.5,7*7), nrow=7); diag(block) = 1
sigma[1:7, 1:7] = block

truePara = list("beta" = trueBeta, "sigma" = sigma)


truePara$gamma = as.numeric(truePara$beta != 0)




G = diag(p)
#block = matrix(rep(1, 4*4), nrow=4); diag(block) = 1
#G[2:5, 2:5] = block
#G[6,7] = G[7,6] = 1
#G[6:7,50] = G[50,6:7] = 1
#G[2:5,25] = G[25,2:5] = 1
#G[1,12] = G[12,1] = 1

if (graph == "true") {
  G = matrix(data = as.numeric( sigma != 0 ), nrow=p, ncol=p)
  diag(G) = 0  
}

if (graph == "empty") {
  G = matrix(0, nrow=p, ncol=p)  
}


truePara = list("beta" = trueBeta, "sigma" = sigma)





S = 1

priorParaPooled = list(#"eta0"   = eta0,                   # prior of baseline hazard
                       #"kappa0" = kappa0,                 # prior of baseline hazard
                       "c0"     = 2,                      # prior of baseline hazard
                       "tau"    = 0.0375,                 # standard deviation for prior of regression coefficients
                       "cb"     = 20,                     # standard deviation for prior of regression coefficients
                       "pi.ga"  = 0.02, #0.5, ga.pi,                   # prior variable selection probability for standard Cox models
                       "a"      = -4, #a0,                     # hyperparameter in MRF prior
                       "b"      = 0.1, #b0,                     # hyperparameter in MRF prior
                       "G"       = G
) 






Beta.ini = rep(0, p)


beta.ini  = Beta.ini 
gamma.ini = rep(0,p)

initial = list("gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = 0)




n = 100
#mcmc_iterations = 30000

#run main MCMC simulation:
# fit = run_simulation(n=n,
#                      truePara=truePara, #true parameters
#                      priorPara = priorParaPooled, #prior parameters
#                      initial = initial, #initial values of the parameters
#                      num.reps = mcmc_iterations,
#                      seed=seed)

#dataset = make_dataset(n, truePara)
#load(file="SimulationStudy/sparse_largeN/dataset_n100_p200.RData")
#load(file="SimulationStudy/sparse_smallN/dataset_sparse_smallN.RData")
load(file="dataset_n100_p200.RData")
log.like  = coxph( Surv(dataset$time.train, dataset$status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariates
initial$log.like.ini = log.like

simulation = simulate(data=dataset,
                      priorPara = priorParaPooled,
                      initial=initial,
                      num.reps=mcmc_iterations,
                      seed=seed)

result = simulation[c("gamma.p", "beta.p", "h.p", "s", "log.jpost", "log.like", "post.gamma", "accept.RW")]


simulation_result = list("truePara"=truePara,
                         "priorPara" = priorParaPooled,
                         "initial" = initial,
                         "result" = result,
                         "mcmcIterations" = mcmc_iterations,
                         "n" = n,
                         "seed" = seed
                         )

#save(dataset, file="SimulationStudy/sparse_largeN/dataset_n100_p200.RData")

#save(simulation_result, file=sprintf("SimulationStudy/sparse_smallN/N=100_P=100/simulation_results%s_empty_N100_P100_MCMClong.RDdata", seed) )
save(simulation_result, file=sprintf("simulation_results/simulation_results%s_%s_N100_P200.RDdata", seed, graph) )

