# Cox model with Bayesian Variable Selection 

This code is based on the R source code and data associated with the publication *Madjar K, Zucknick M, Ickstadt K, and Rahnenführer J (2020): Combining heterogeneous subgroups with graph-structured variable selection priors for Cox regression. arXiv: 2004.07542*.

This method is focused on the situation of predefined, possibly heterogenous subgroups of patients with available survival endpoint and continuous molecular measurements such as gene expression
data, with the aim of obtaining a separate risk prediction model for each subgroup.
For this purpose, we propose a Cox regression model with a graph-structured variable selection prior that incorporates information on the relationships among the covariates and encourages the joint selection of linked variables.
These links correspond to variables either being conditionally dependent within each subgroup (e.g. functional or regulatory pathways) or being simultaneously prognostic across different subgroups.

We compare our approach to Cox regression models with an independent Bernoulli prior for variable selection as proposed by Treppmann et al. (2017).
We evaluate all models through simulations and a case study with Glioblastoma protein expression data.


The difference between the method in Madjar et. al (2021) and my master thesis, will be that we will not focus on heterogenous subgroups, and we will not do structure learning on the graph.
We will use Stochastic Search Variable Selection (SSVS), just as in Madjar et. al (2021), with a Markov Random Field prior, where we keep the graph fixed. This will allow us to use the
graph to incorporate biological knowledge about funcitional and regulatory molecular pathways, and encourage the joint inclusion of genes that are in the same pathway.


The main file to run the MCMC simulation is **Run_func_MCMC.R**

The main files to run the simulation study from Madjar et. al (2021) is **Run_Simulation.R** and the main file for the case study is **Run_CaseStudy.R**. Files from
this paper might be removed if they are not necessary to my thesis.

## Overview of R files:

#### Run_func_MCMC.R
Running the main MCMC simulation. Currently using simulated data

#### Run_Simulation.R

Setup and running of the simulations.
All simulations are run using the R package batchtools for parallelization.


#### Run_CaseStudy.R

Setup and running of the case study.
All settings are run using the R package batchtools for parallelization.

#### Data_Simulation.R
Helper functions for generation of simulated training and test data sets, including simulated gene expression data and survival outcome. 
Only needed in the simulation study, not in the case study.

#### WrapperSimCoxBVSSL.R
Wrapper function for the inference of a specific type of Cox model with MCMC sampling based on simulated data.
First, the training data is prepared (standarize covariates for Pooled model), then the hyperparameters for all priors are defined and starting values 
for the MCMC sampler are set.

#### WrapperCaseCoxBVSSL.R
Analogous to WrapperSimCoxBVSSL.R but for case study data, wrapper function for the inference of a specific type of Cox model with MCMC sampling.

#### func_MCMC.R
Main function for performing MCMC sampling of all parameters and storage of posterior results.

#### func_MCMC_cox.R
Helper functions for MCMC sampling of Cox model parameters and joint posterior distribution.

#### func_MCMC_graph.R
Helper function for MCMC sampling of graph and precision matrix (structure learning, adapted code from Wang (2015).


#### Weibull_param.RData
RData object with parameters of the Weibull distribution used for data simulation (see file "Data_Simulation.R")

#### GBM_Data.RData

RData object with preprocessed Glioblastoma protein expression data for case study, including original expression data of 20 selected proteins, 
simulated survival endpoint, subgroup membership for each patient.


Start with file Run_Simulation or Run_CaseStudy.R.
