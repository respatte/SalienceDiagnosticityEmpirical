# LIBRARY IMPORTS ==================================================================================
library(brms)
library(tidyverse)

# BRM FIXEF BAYES FACTOR
# Function computing all nested models in formulas nd computing Bayes factors for pairs of models
# Formulas must be given with no-intercept formulas first, then intercept-only formulas, and
# finally other formulas.
# no_intercept is the number of formulas with no intercept (default = 0)
# intercept_only is the number of formulas with intercept-only (default = 1)
# priors is a list with the prior for intercept-only models and priors for other models
# (no prior for no-intercept models)
bayes_factor.brm_fixef <- function(formulas, df, priors,
                                   no_intercept = 0, intercept_only = 1,
                                   controls = NULL, iter = 2000, family = gaussian()){
  ## Run all Bayesian models
  models <- lapply(seq_along(formulas),
                   function(i){
                     brm(formulas[[i]],
                         data = df,
                         prior = if(i <= no_intercept){
                           NULL
                         }else{
                           if(i <= no_intercept + intercept_only){
                             # Intercept-only model, no slope priors
                             priors[[1]]
                           }else{
                             # Slope priors
                             priors[[2]]
                           }
                         },
                         family = family,
                         chains = 4, cores = 4, iter = iter,
                         control = controls,
                         save_all_pars = T)
                   })
  ## Compute all Bayes fators
  n_models <- length(formulas)
  bayes_factors <- lapply(2:n_models,
                          function(i){
                            bayes_factor(models[[i]],
                                         models[[i-1]])
                          })
  ## Return both lists as a list
  return(list("models" = models, "BF" = bayes_factors))
}
