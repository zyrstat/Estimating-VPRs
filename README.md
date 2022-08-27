# Estimating_VPRs
We used the publicly available epidemiological data to estimate the vaccine protection rates with the methods provided in "Real-World COVID-19 Vaccine Protection Rates against Infection in the Delta and Omicron Eras" R version 4.2.0 was used for data analysis. 
R packages Rcpp and RcppArmadillo were used in the estimation of model parameters.

First, use Rcpp and RcppArmadillo to build the R packages named "objfunction" and "objectfunctionsbooster", which contain functions for the estimation. 
Then, use "real data Omicron-reinfection.R" to implement the estimation.
