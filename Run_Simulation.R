# Run both simulations with all settings described in the paper

library("batchtools")

# Batchtools registry to save all information regarding computational jobs and all results
reg =  makeExperimentRegistry(file.dir = "CoxBVSSL_Sim") 

# "Problem" function in batchtools that creates problem instance (here: generation of training
# and test data) used in the "algorithmic" function below. 
prob.fun = function( p, n, ...){
  source("Data_Simulation.R")
  Data = sim_data_fun(n = n, p = p,
                      # Weibull parameters for simulation of event and censoring times in both subgroups
                      # (we use the same parameters for event and censoring times to obtain about 
                      # 50% censoring rates):
                      surv.e = Surv.e, surv.c = Surv.e, 
                      # Set true effects of the first prognostic genes in both subgroups:
                      beta1.p = c(1,1,1,-1,-1,-1,0,0,0), beta2.p = c(0,0,0,-1,-1,-1,1,1,1)
  )
}

###
# Simulation I

# Parameter combinations for data simulation
prob.designs = vector("list", 1)
vp = rbind(  CJ(n = c(50,75,100,125,150), p = 100),
             CJ(n = c(50,75,100), p = 20) 
)
prob.designs[[1]] = vp
names(prob.designs)[1] = "Simdata"

addProblem(name  = "Simdata",
           data  = NULL,
           fun   = prob.fun,
           seed  = 283 
)

source("WrapperSimCoxBVSSL.R")  # Algorithm function 
addAlgorithm(name  = "SimCoxBVSSL", fun = WrapperSimCoxBVSSL ) 

m.seed = 2832000 

# Parameters for algorithm function: different Cox model types
design = CJ(model.type = c("CoxBVSSL", "Pooled", "Subgroup"),
            MRF2b = FALSE, b01 = 1, b02 = 1, m.seed = m.seed )  

# Add experiments to batchtools registry (combinations of problems and algorithms) to define jobs
addExperiments(prob.designs  = prob.designs,
               algo.designs  = list(SimCoxBVSSL = design ),
               repls = 10L  
)   
# Adding 240 experiments


###
# Add further simulations to registry for Simulation II

# Parameter combinations for data simulation
prob.designs = vector("list", 1)
vp = CJ(n = c(50,75,100,125,150), p = 100)
prob.designs[[1]] = vp
names(prob.designs)[1] = "Simdata"

design = rbind(
  CJ(model.type = "Sub-struct", 
     MRF2b = FALSE, b01 = 1, b02 = 1, m.seed = m.seed ),  
  CJ(model.type = "CoxBVSSL", 
     MRF2b = TRUE, b01 = 1, b02 = c(1,1.5,2,2.5,3), m.seed = m.seed)  
)

addExperiments(prob.designs  = prob.designs,
               algo.designs  = list(SimCoxBVSSL = design ),
               repls = 10L  
)
# Adding 300 experiments


###

submitJobs(ids = findNotSubmitted())  # Run (all) jobs or run single job with job ID 1 (example):
submitJobs(1)
