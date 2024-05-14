# Discovery of prognostic biomarkers using a Bayesian Cox model with molecular pathway information

This code is based on the R source code and data associated with the publication *Madjar K, Zucknick M, Ickstadt K, and Rahnenführer J (2020): Combining heterogeneous subgroups with graph-structured variable selection priors for Cox regression. arXiv: 2004.07542*.

This method is focused on the situation of predefined, possibly heterogenous subgroups of patients with available survival endpoint and continuous molecular measurements such as gene expression
data, with the aim of obtaining a separate risk prediction model for each subgroup.
For this purpose, we propose a Cox regression model with a graph-structured variable selection prior that incorporates information on the relationships among the covariates and encourages the joint selection of linked variables.
These links correspond to variables either being conditionally dependent within each subgroup (e.g. functional or regulatory pathways) or being simultaneously prognostic across different subgroups.

The difference between the method in Madjar et. al (2021) and my master thesis, will be that we will not focus on heterogenous subgroups, and we will not do structure learning on the graph.
We will use Stochastic Search Variable Selection (SSVS), just as in Madjar et. al (2021), with a Markov Random Field prior, where we keep the graph fixed. This will allow us to use the
graph to incorporate biological knowledge about functional and regulatory molecular pathways, and encourage the joint inclusion of genes that are in the same pathway.

The main files to run the simulation study from Madjar et. al (2021) is **Run_Simulation.R** and the main file for the case study is **Run_CaseStudy.R**. Files from
this paper might be removed if they are not necessary to my thesis.

## Overview of R files:

### Code for the MCMC algorithm
Files that are used for the MCMC simulation are:


#### func_MCMC.R
Main function for performing MCMC sampling of all parameters and storage of posterior results.

#### func_MCMC_cox.R
Helper functions for MCMC sampling of Cox model parameters and joint posterior distribution.


### The simulation studies


#### Data_Simulation.R
Helper functions for generation of simulated training and test data sets, including simulated gene expression data and survival outcome. 
Only needed in the simulation study, not for real data analysis.

#### Weibull_param.RData
RData object with parameters of the Weibull distribution used for data simulation (see file "Data_Simulation.R")

#### Run_Simulation_Study.R
The file used to run the MCMC algorithm for the simulation studies. 

#### generate_plot_info.R
File that generates different values needed for the plots in the thesis from the MCMC samples, such as the MPM coefficient estimates, MPM model size, IBS and more. The plot_info objects for Simulation study 1 and simulation study 2 (partial uniform, partial non-uniform and noise) is in the SimStudy folder. As the MCMC samples of all runs for the simulation chapter went over the file limit of GitHub, they are not included.

#### final_plots.R
Code for the plots in the Simulations chapters.


### Application: The TCGA breast cancer data set

#### real_data_download.R
this file downloads the gene expression and mutation data from the TCGA-BRCA data set.

#### real_data_preprocessing.R
In this file, the preprocessing of the gene expression and mutation features is done, as well as the pre-filtering of features. The construction of the graph in the MRF prior is also done here, from the KEGG breast cancer pathways.

#### real_data_descriptive.R
File generating plots and tables for the TCGA breast cancer data set, such as estimated survival probabilities, and Kaplan-Meier estimates of survival curves of treatment/control group, different breast cancer subtypes etc.

#### run_real_data.R
File that runs the MCMC algorithm for the Application chapter of the thesis.

#### generate_plot_info_real_data.R
file that generates a plot_info object for the real data analysis, containing information such as the MPM coefficient estimates, MPM model size and the IBS. The plot_info object for the real data analysis is included in the RealData folder in this repo. All MCMC samples (after thinning, to reduce file size) used in the Application chapter is also inlcuded in the RealData folder. The objects simulation_result in those output-files contain all prior parameters and the random seed, so the results can be reproduced.


#### real_data_plots.R
File generating plots from the Application chapter of the thesis.



