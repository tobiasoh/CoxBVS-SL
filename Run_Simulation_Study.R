library(survival)
library(RhpcBLASctl)


args = commandArgs(trailingOnly = T)

server = "server" %in% args
path = "./"
if (server) {
  path = "/data/tobiasoh/"
}
source(sprintf("%sSimulation_Study.R", path))
source(sprintf("%sbrier_calculations.R", path))
source(sprintf("%splot_helper.R", path))

load(args)
graph = 

RhpcBLASctl::blas_set_num_threads(20)
#setting the true parameters, that we use to simulate the data set

graph = args[1]
mcmc_iterations = as.integer(args[2])
if (any( lapply("t[0-9]+", grepl, args)[[1]] ) ) {
  index = which( lapply("t[0-9]+", grepl, args)[[1]] ) 
  #thinning = as.integer(gregexpr("[0-9]+", args[index]))
  thinning = as.integer(gsub("\\D", "", args[index]))
} else {
  thinning = 1
}


if (any( lapply("d[0-9]+", base::grepl, args)[[1]] ) ) {
  index = which( lapply("d[0-9]+", grepl, args)[[1]] ) 
  num_dataset = as.integer(gsub("\\D", "", args[index]))
} else {
  num_dataset = 1
}




#trueBeta = runif(20, min=-1, max=1)
#trueBeta = c(1,1,1,1,-1,-1,-1)
#trueBeta = c(trueBeta, rep(0, p - length(trueBeta)))

#sigma = diag(p)
#block = matrix(rep(.5,9), nrow=3); diag(block) = 1
#sigma[1:3, 1:3] = block#sigma[4:6, 4:6] = sigma[7:9, 7:9] = block


#sigma = diag(p)
#block = matrix(rep(.5,15*15), nrow=15); diag(block) = 1
#sigma[1:15, 1:15] = block

#truePara = list("beta" = trueBeta, "sigma" = sigma)
#mcmc_iterations = 1000
#graph="noise1"
seed = sample(1:10e7, 1)
set.seed(seed)

p = 200

if (server) {
  load("/data/tobiasoh/truePara.RData")
} else {
  load("./SimStudy/truePara.RData")
}

truePara$gamma = as.numeric(truePara$beta != 0)


G = matrix(data = as.numeric( truePara$sigma != 0 ), nrow=p, ncol=p)
diag(G) = 0 
#block = matrix(rep(1, 4*4), nrow=4); diag(block) = 1
#G[2:5, 2:5] = block
#G[6,7] = G[7,6] = 1
#G[6:7,50] = G[50,6:7] = 1
#G[2:5,25] = G[25,2:5] = 1
#G[1,12] = G[12,1] = 1
sigma = truePara$sigma

if (grepl("true", graph, fixed=T)) {
  G = matrix(data = as.numeric( sigma != 0 ), nrow=p, ncol=p)
  diag(G) = 0  
}

if (grepl("empty", graph, fixed=T)) {
  G = matrix(0, nrow=p, ncol=p)  
}
if (grepl("partial_uniform", graph, fixed=T)) {
  G = matrix(data = as.numeric( sigma != 0 ), nrow=p, ncol=p)
  diag(G) = 0  
  vectorised = c(G)
  removed_last = FALSE
  
  for (i in 1:length(vectorised)) {
    if (vectorised[i] == 1) {
      if (!removed_last) {
        vectorised[i] = 0
        removed_last = TRUE
      }
      else {
        removed_last = FALSE
      }
    }
  }
  
  G = matrix(vectorised, nrow=p, ncol=p)
  
  
}


if (grepl("partial_non_uniform", graph, fixed=T)) {
  #vectorised = c(G)
  #edges = which(G == 1, arr.ind=T)
  #edges_remove = sample(unique(edges[,1]), 5)
  #G[edges_remove,] = 0
  #G[,edges_remove] = 0
  num_vars = as.integer(gsub("\\D", "", graph))
  
  
  G[1:num_vars, 1:num_vars] = 0
  
}

if (grepl("noise", graph, fixed=T)) {

  num_noise = as.integer(gsub("\\D", "", graph))
  if (num_noise == 100) {
    load(sprintf("%s/SimStudy/real_sim/noise_graph.RData", path))
  }
  if (num_noise == 50) {
    load(sprintf("%s/SimStudy/real_sim/noise_graph50.RData", path))
  }
  if (num_noise == 200) {
    load(sprintf("%s/SimStudy/real_sim/noise_graph200.RData", path))
  }
  
  
  # num_edges = sum(G == 1)
  # edges_added = 0
  # 
  # while (edges_added < num_edges) {
  #   indices = sample.int(p, 2)
  # 
  # 
  #   if (G[indices[1], indices[2]] == 0) {
  #     G[indices[1], indices[2]] = G[indices[2], indices[1]] = 1
  #     edges_added = edges_added + 1
  #   }
  # }

}
if (grepl("line", graph, fixed=T)) {
  degree = as.integer(gsub("\\D", "", graph)) #9 gives 10%, 3 ish 25%. 2 gives 1/3. 15 ish 6%
  G = matrix(data = as.numeric( sigma != 0 ), nrow=p, ncol=p)
  diag(G) = 0  
  vectorised = c(G)
  num_removed = 1
  
  for (i in 1:length(vectorised)) {
    if (vectorised[i] == 1) {
      if (num_removed == degree) {
        #vectorised[i] = 0
        num_removed = 0
      }
      else {
        vectorised[i] = 0
        num_removed = num_removed + 1
      }
    }
  }
  
  G = matrix(vectorised, nrow=p, ncol=p)
  #G[1,15] = G[15,1] = 1
}
if (grepl("non_uniform_max", graph, fixed=T)) {
  num_vars = as.integer(gsub("\\D", "", graph))
  G = matrix(data = as.numeric( sigma != 0 ), nrow=p, ncol=p)
  diag(G) = 0  
  G[1:num_vars,] = 0
  G[,1:num_vars] = 0
  
}

#truePara = list("beta" = trueBeta, "sigma" = sigma)





S = 1

priorParaPooled = list(#"eta0"   = eta0,                   # prior of baseline hazard
                       #"kappa0" = kappa0,                 # prior of baseline hazard
                       "c0"     = 2,                      # prior of baseline hazard
                       "tau"    = 0.0375,                 # standard deviation for prior of regression coefficients
                       "cb"     = 20,                     # standard deviation for prior of regression coefficients
                       "pi.ga"  = 0.02, #0.5, ga.pi,                   # prior variable selection probability for standard Cox models
                       "a"      = -1, #a0,                     # hyperparameter in MRF prior
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

####dataset = make_dataset(n, truePara)

#num_dataset = gsub("\\D", "", graph)
#num_dataset= as.integer(gsub("\\D", "", args[3]))
load(sprintf("%s/SimStudy/datasets/dataset%s.RData", path, num_dataset))


#load(file="SimulationStudy/sparse_largeN/dataset_n100_p200.RData")
#load(file="SimulationStudy/sparse_smallN/dataset_sparse_smallN.RData")
#load(file="dataset_n100_p200.RData")
log.like  = coxph( Surv(dataset$time.train, dataset$status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariates
initial$log.like.ini = log.like

simulation = simulate(data=dataset,
                      priorPara = priorParaPooled,
                      initial=initial,
                      num.reps=mcmc_iterations,
                      thinning = thinning,
                      seed=seed)

result = simulation[c("gamma.p", "beta.p", "h.p", "s", "log.jpost", "log.like", "post.gamma", "accept.RW")]


#result=simulation_result$result
#load(file="SimulationStudy/sparse_largeN/dataset_n100_p200.RData")




#surv_prob_MCMC = survProb_MCMC(result$beta.p[-(1:warmup),], dataset$X.train, result$h.p[-(1:warmup),])

#brier_mcmc = BrierScore2(result$beta.p[-(1:warmup),], dataset$X.train, survival_data, surv_prob_MCMC, result$s)




# if (thinning > 1) {
#   result$gamma.p = result$gamma.p[seq(2,mcmc_iterations+1, thinning),]
#   result$beta.p = result$beta.p[seq(2,mcmc_iterations+1, thinning),]
#   result$h.p = result$h.p[seq(2,mcmc_iterations+1, thinning),]
#   result$post.gamma = result$post.gamma[seq(1,mcmc_iterations, thinning),]
#   result$accept.RW = result$accept.RW[seq(2,mcmc_iterations+1, thinning),]
#   
#   result$log.jpost = result$log.jpost[seq(2,mcmc_iterations+1, thinning)]
#   result$log.like = result$log.like[seq(2,mcmc_iterations+1, thinning)]
# }

#need to adapt warmup to the thinnning parameter. Wants half of iterations kept as warmup
warmup = ceiling( length(seq(1,mcmc_iterations, thinning))/2 )#ceiling(mcmc_iterations/2)

#metrics = sensitivity_and_specificity(result, p, warmup)


  
  #beta_estimate = variableSelection(simulation_result, warmup)
  
  #survival_prob = calculateSurvivalProb(beta_estimate, dataset$X.train, result$h.p)
  
  survival_data = list("time.train"=dataset$time.train, "status.train"=dataset$status.train)


time_points =  seq(0,5,0.1)
#brier_score = BrierScoreVectorised(result$beta.p[-(1:warmup),], dataset$X.train, survival_data, time_points, result$h.p[-(1:warmup),], result$s)

#ibs = integratedBrierScore(brier_score, time_points)


simulation_result = list("truePara"=truePara,
                         "priorPara" = priorParaPooled,
                         "initial" = initial,
                         "result" = result,
                         "mcmcIterations" = mcmc_iterations,
                         "n" = n,
                         "seed" = seed,
                         #"brier_score" = brier_score,
                         #"ibs" = ibs,
                         "thinning" = thinning,
                         #"metrics" = metrics,
                         "survival_data" = survival_data,
                         "X.train" = dataset$X.train
)
#save(simulation_result, file=sprintf("%sSimStudy/real_sim/%s.RData", path, graph))

p=200
simulation_result$metrics = sensitivity_and_specificity(simulation_result, p, warmup)


#save(dataset, file="SimulationStudy/sparse_largeN/dataset_n100_p200.RData")

#save(simulation_result, file=sprintf("SimulationStudy/sparse_smallN/N=100_P=100/simulation_results%s_empty_N100_P100_MCMClong.RDdata", seed) )

#save(simulation_result, file=sprintf("/data/tobiasoh/simulation_results/%s.RData", graph) )

save(simulation_result, file=sprintf("%s/SimStudy/real_sim_thin6/%s_%d.RData", path, graph, num_dataset))
