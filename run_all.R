## Convenience driver: runs the full analysis end to end.
## Open kisumu_variability.Rproj first (so paths resolve), then:
##   source("run_all.R")
##
## Step 1 fits all four insecticides (SLOW - Bayesian MCMC, may take
## a while and uses multiple cores). Step 2 needs Step 1's outputs.

source(here::here("scripts", "1_variability_analysis.R"))
source(here::here("scripts", "2_spearman_correlation.R"))
