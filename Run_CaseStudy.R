# Run case study with all settings described in the paper

library("batchtools")

# Batchtools registry to save all information regarding computational jobs and all results
reg =  makeExperimentRegistry(file.dir = "CoxBVSSL_CaseStudy")

# We have one data set and sample different training and test indices (subsampling).
# "Problem" function in batchtools that creates problem instance (here: performs 
# subsampling) used in the "algorithmic" function below. 
subsampling.fun <- function(data, ratio, ...) {   
  study = levels(data$group)  
  train = list()
  test  = list()
  for(i in 1:length(study)) {
    index  = which( data$group == study[i] ) # samples belonging to i-th subgroup/ study
    status = data$surv[index, 2L]  
    # Stratified sampling of training indices, stratified according to subgroup and event
    # indicator
    train[[i]] = unlist(lapply(split(index, status),  
                                function(id) sample(id, floor(length(id) * ratio))), 
                         use.names = FALSE )
    test[[i]]  = setdiff(index, train[[i]])
  }
  return(list(train = unlist(train), # training indices across all subgroups
              test = unlist(test)
  ))
}

# Parameters for subsampling
prob.designs = vector("list", 1)
prob.designs[[1]] = data.frame(ratio = .8) 
names(prob.designs)[1] = "GBdata"

load("GBM_Data.RData") # Glioblastoma data

addProblem(name  = "GBdata",
           data  = Data,
           fun   = subsampling.fun,
           seed  = 293 
)

source("WrapperCaseCoxBVSSL.R")  # Algorithm function 
addAlgorithm(name  = "CaseCoxBVSSL", fun = WrapperCaseCoxBVSSL ) 

m.seed = 2932000 

# Parameters for algorithm function: different Cox model types

design = CJ( model.type = c("CoxBVSSL", "Pooled", "Subgroup"),
             MRF2b = FALSE, b01 = 0.5, b02 = 0.5, m.seed = m.seed
)
design = rbind(design, 
               CJ(model.type = "CoxBVSSL", 
                  MRF2b = TRUE, b01 = 0.5, b02 = c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3), m.seed = m.seed
) )
design = rbind(design, 
               CJ(model.type = "CoxBVSSL", 
                  MRF2b = TRUE, b01 = c(1,1.5,2,2.5,3), b02 = 0.5, m.seed = m.seed
) )

# Add experiments to batchtools registry (combinations of problems and algorithms) to define jobs
addExperiments(prob.designs  = prob.designs,
               algo.designs  = list(CaseCoxBVSSL = design ),
               repls = 10L  
) 


###
# Run all jobs:

ids1 = findExperiments(algo.pars = (model.type != "CoxBVSSL"))
ids2 = findExperiments(algo.pars = (model.type == "CoxBVSSL"))

submitJobs( ids1, resources = list(walltime = 2*60*60, memory = 8192L, measure.memory = TRUE) )
submitJobs( ids2, resources = list(walltime = 2*60*60, memory = 20000L, measure.memory = TRUE) )
