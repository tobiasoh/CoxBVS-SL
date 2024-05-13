time_start = Sys.time()
library(survival)
library(RhpcBLASctl)

seed = sample(1:10e7, 1)
set.seed(seed)
args = commandArgs(trailingOnly = T)

server = "server" %in% args
path = "./"
if (server) {
  path = "/data/tobiasoh/"
}
source(sprintf("%sSimulation_Study.R", path))
source(sprintf("%sbrier_calculations.R", path))
source(sprintf("%splot_helper.R", path))


RhpcBLASctl::blas_set_num_threads(20)
#setting the true parameters, that we use to simulate the data set


load(args[1]) #dataset
load(args[2]) #graph G

filename = args[3]
last_part = substr(filename, nchar(filename) - 1, nchar(filename))

digits = gsub("\\D", "", last_part)

num_dataset = as.integer(digits)

mcmc_iterations = as.integer(args[4])

if (any( lapply("-t[0-9]+", grepl, args)[[1]] ) ) {
  index = which( lapply("-t[0-9]+", grepl, args)[[1]] ) 
  #thinning = as.integer(gregexpr("[0-9]+", args[index]))
  thinning = as.integer(gsub("\\D", "", args[index]))
} else {
  thinning = 1
}

if (any( lapply("-b[0-9]+", grepl, args)[[1]] ) ) {
  index = which( lapply("-b[0-9]+", grepl, args)[[1]] ) 
  #thinning = as.integer(gregexpr("[0-9]+", args[index]))
  #b = as.integer(gsub("\\D", "", args[index]))
  splt = strsplit(args[index], "b")
  b = as.numeric(splt[[1]][2])
} else {
  b = 0.1
}

if (any( lapply("-a[0-9]+", grepl, args)[[1]] ) ) {
  index = which( lapply("-a[0-9]+", grepl, args)[[1]] ) 
  #thinning = as.integer(gregexpr("[0-9]+", args[index]))
  #b = as.integer(gsub("\\D", "", args[index]))
  splt = strsplit(args[index], "a")
  a = -1*abs( as.numeric(splt[[1]][2]) )
} else {
  a = -3
}


#80-20 training split
#training_percentage = 0.8
#num_patients = nrow(X)
#test_indices = list()


# for (i in 1:20) {
#   train_indices = sample(1:num_patients, size=training_percentage*num_patients)
#   
#   
#   test_indices = append( test_indices, list( dplyr::setdiff(1:nrow(X), train_indices) ) )
#   
#   
# }

#train_indices = sample(1:num_patients, size=training_percentage*num_patients)
#test_indices = dplyr::setdiff(1:nrow(X), train_indices) 
#load( sprintf("%sRealData/datasplits%s.RData") ) #load list of test indices

load(sprintf("./RealData/datasplits/test_indices%i.RData", num_dataset)) # loading dataset (object test_indices)
train_indices = dplyr::setdiff(1:nrow(X), test_indices) 
X.train = X[train_indices,]
X.test = X[test_indices,]

#save(test_indices, file="./RealData/datasplits.RData")
#meta$time = meta[which(meta$vital_status=="Dead"),]$days_to_death
#meta[which(is.na(meta$time)), "time"] = meta[which(is.na(meta$time)), "days_to_last_follow_up"]



time.train = meta[train_indices, "time"]
time.test = meta[test_indices, "time"]

status.train = as.integer(meta[train_indices,]$vital_status == "Dead")
status.test = as.integer(meta[test_indices,]$vital_status == "Dead")

dataset = list("X.train" = X.train,
               "X.test" = X.test,
               "time.train" = time.train,
               "status.train" = status.train,
               "status.test" = status.test)

#maybe change a to -3, 

priorPara = list(#"eta0"   = eta0,                   # prior of baseline hazard
  #"kappa0" = kappa0,                 # prior of baseline hazard
  "c0"     = 2,                      # prior of baseline hazard
  "tau"    = 0.0375,                 # standard deviation for prior of regression coefficients
  "cb"     = 20,                     # standard deviation for prior of regression coefficients
  "pi.ga"  = 0.02, #0.5, ga.pi,                   # prior variable selection probability for standard Cox models
  "a"      = a,#-4, #a0,                     # hyperparameter in MRF prior
  "b"      = b, #b0,                     # hyperparameter in MRF prior
  "G"       = G
) 

print(cat("a: ", a))
print(cat("b: ",b))

p = ncol(X.train)
beta.ini = rep(0, p)
gamma.ini = c(rep(0,p-2), 1, 1) #for clin. Age and treatment included in mod from initial
#gamma.ini[vars_in_graph] = 1
initial = list("gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = 0)



log.like  = coxph( Surv(dataset$time.train, dataset$status.train, type = c('right')) ~ 1 )$loglik # initial value: null model without covariates
initial$log.like.ini = log.like

library(Rcpp)
library(RcppArmadillo)
sourceCpp("./BayesSurv/src/updateRP_genomic_cpp.cpp")

UpdateRP.lee11.helper_c = updateRP_genomic_cpp
#print(dim(X.train))

print(thinning)

simulation = simulate(data=dataset,
                      priorPara = priorPara,
                      initial=initial,
                      num.reps=mcmc_iterations,
                      thinning=thinning,
                      seed=seed) 

result = simulation[c("gamma.p", "beta.p", "h.p", "s", "log.jpost", "log.like", "post.gamma", "accept.RW")]

# 
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



model_output = list(
                         "priorPara" = priorPara,
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
                         "test_indices" = test_indices
)



#if (grepl("full", args[2]) ) {
#  filename = "full"
  
#}else if (grepl("empty", args[2])) {
#  filename = "empty"
  
#}else {
#  filename="kegg"
#}


save(model_output, file=sprintf("%s/RealData/%s.RData", path, filename))
print(time_start)
print(Sys.time())
