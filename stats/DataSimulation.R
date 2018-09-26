# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(lmerTest)
library(emmeans)
library(brms)
library(coda)
library(tidyverse)

# BAYESIAN POINT NULL HYPOTHESIS TESTING ===========================================================
# Generate random data
# n=42 * 2trials from real data, mean to 0, sd=.5 from real data
new_old.sims <- rep(list(tibble(ChanceArcsin = rnorm(42*2, 0, .5),
                                Participant = rep(1:42, each = 2),
                                Condition = factor(rep(c("Label", "No Label"), times = 21, each = 2)),
                                ContrastType = factor(rep(c("Head", "Tail"), times = 42)))),
                    times = 10)

# Define STB analysis function, returning emmeans analysis
stb.analysis <- function(df){
  m <- lmer(ChanceArcsin ~ ContrastType*Condition +
              (1 | Participant),
            data = df)
  t <- emmeans(m, ~ ContrastType | Condition,
               options = list(infer = c(T, T), null = 0,
                              level = .89))
  return(as_tibble(t))
}

# Define Bayesian analysis function, returning emmeans bf
bayesian.analysis <- function(df){
  p <- c(set_prior("uniform(-.8,.8)",
                   class = "Intercept"),
         set_prior("normal(0,.5)", class = "b"))
  m <- brm(ChanceArcsin ~ ContrastType + Condition +
             ContrastType:Condition +
             (1 | Participant),
           data = df, prior = p, family = gaussian(),
           chains = 4, cores = 4, iter = 2000,
           save_all_pars = T)
  bf <- hypothesis(m,
                   c("Intercept > 0",
                     "Intercept + ContrastTypeTail > 0",
                     "Intercept + ConditionNo Label > 0",
                     paste("Intercept +",
                           "ConditionNo Label +",
                           "ContrastTypeTail +",
                           "ContrastTypeTail:ConditionNo Label",
                           "> 0")),
                   alpha = .11)
  return(as_tibble(bf$hypothesis))
}

# Define summary functions for stats results
stb.summary <- function(l){
  e <- lapply(l,
              function(t){sum(t$p.value < .05)})
  return(unlist(e))
}
bayesian.summary <- function(l){
  e <- lapply(l,
              function(t){sum(t$Evid.Ratio > 3)})
  return(unlist(e))
}

# Get statistics and evidence summary from simulations for STB
t <- proc.time()
new_old.sims.stb <- lapply(new_old.sims,
                           stb.analysis)
new_old.sims.stb.evid <- stb.summary(new_old.sims.stb)
stb.time <- proc.time() - t

# Get statistics and evidence summary from simulations for Bayesian
t <- proc.time()
new_old.sims.bayesian <- lapply(new_old.sims,
                                bayesian.analysis)
new_old.sims.bayesian.evid <- bayesian.summary(new_old.sims.bayesian)
bayesian.time <- proc.time() - t