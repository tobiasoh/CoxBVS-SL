# BayesSurv


This is a R/Rcpp package **BayesSurv** for a Bayesian Cox model with graph-structured selection priors for sparse identification of omics features predictive of survival ([Madjar et al., 2021](https://doi.org/10.1186/s12859‐021‐04483‐z)) and its extension with the use of a fixed graph via a Markov Random Field (MRF) prior for capturing known structure of omics features, e.g. disease-specific pathways from the Kyoto Encyclopedia of Genes and Genomes (KEGG) database.

## Installation

Install the latest development version from GitHub

```r
#options(timeout = max(300, getOption("timeout")))
#install.packages("remotes")
remotes::install_github("tobiasoh/master_thesis/BayesSurv")
```

## Examples

### Simulate data

```r
# Load the example dataset
data("simData", package = "BayesSurv")
dataset = list("X" = simData[[1]]$X, 
               "t" = simData[[1]]$time,
               "di" = simData[[1]]$status)
```


### Run a Bayesian Cox model

```r
## Initial value: null model without covariates
library("survival")
log.like  = coxph( Surv(dataset$t, dataset$di, type = c('right')) ~ 1 )$loglik 
initial = list("gamma.ini" = rep(0, ncol(dataset$X)), 
               "beta.ini" = rep(0, ncol(dataset$X)), 
               "log.like.ini" = log.like)
# Prior parameters
priorParaPooled = list(
  #"eta0"   = eta0,                   # prior of baseline hazard
  #"kappa0" = kappa0,                 # prior of baseline hazard
  "c0"     = 2,                      # prior of baseline hazard
  "tau"    = 0.0375,                 # sd for coefficient prior
  "cb"     = 20,                     # sd for coefficient prior
  "pi.ga"  = 0.02, #0.5, ga.pi,      # prior variable selection probability for standard Cox models
  "nu0"    = 0.05,                   # hyperparameter in graphical model
  "nu1"    = 5,                      # hyperparameter in graphical model
  "lambda" = 3,                      # hyperparameter in graphical model
  "a"      = -4, #a0,                # hyperparameter in MRF prior
  "b"      = 0.1, #b0,               # hyperparameter in MRF prior
  "G"      = simData$G               # hyperparameter in MRF prior
)   

## run Bayesian Cox with graph-structured priors
library("BayesSurv")
fit = BayesSurv(survObj=dataset, priorPara=priorParaPooled, 
                initial=initial, nIter=100, seed=123)

## show posterior mean of coefficients and 95% credible intervals
library(ggplot2)
plot(fit) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 7))

#plot(fit$output$beta.p[,1], type="l")
#colMeans(fit$output$beta.p[-c(1:(nrow(fit$output$beta.p)/2)), ])
#colMeans(fit$output$gamma.p[-c(1:(nrow(fit$output$gamma.p)/2)), ])
#trueBeta
```

<img src="man/figures/README_plot_beta.png" width="100%" />


## References

> Katrin Madjar, Manuela Zucknick, Katja Ickstadt, Jörg Rahnenführer (2021).
> Combining heterogeneous subgroups with graph‐structured variable selection priors for Cox regression.
> _BMC Bioinformatics_, 22(1):586. DOI:[10.1186/s12859‐021‐04483‐z](https://doi.org/10.1186/s12859‐021‐04483‐z).
