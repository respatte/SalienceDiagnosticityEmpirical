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
                    times = 100)

# Define STB analysis function, returning emmeans analysis
stb.analysis <- function(df.list){
  e <- lapply(df.list,
              function(df){
                m <- lmer(ChanceArcsin ~ ContrastType*Condition +
                       (1 | Participant),
                     data = df)
                t <- emmeans(m, ~ ContrastType | Condition,
                          options = list(infer = c(T, T), null = 0,
                                         level = .89)) %>%
                  as_tibble()
                #return(sum(t$p.value < .05))
                return(t)
              })
  return(bind_rows(e))
}

# Define Bayesian analysis function, returning emmeans bf
bayesian.analysis <- function(df.list){
  p <- c(set_prior("uniform(-.8,.8)",
                   class = "Intercept"),
         set_prior("normal(0,.5)", class = "b"))
  m.list <- brm_multiple(ChanceArcsin ~ ContrastType + Condition +
                      ContrastType:Condition +
                      (1 | Participant),
                    data = df.list, prior = p, family = gaussian(),
                    chains = 4, cores = 4, iter = 2000,
                    save_all_pars = T, combine = F)
  bf.list <- lapply(m.list,
                    function(m){
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
                    })
  #e <- lapply(bf.list,
  #            function(bf){sum(bf$Evid.Ratio > 3)})
  bf <- bind_rows(bf.list)
  return(bf)
}

# Get evidence summary from simulations
t <- proc.time()
new_old.sims.stb.evid <- stb.analysis(new_old.sims)
stb.time <- proc.time() - t
t <- proc.time()
new_old.sims.bayesian.evid <- bayesian.analysis(new_old.sims)
bayesian.time <- proc.time() - t

# Plot p-values
new_old.sims.stb.plot <- ggplot(new_old.sims.stb.evid,
                                aes(y = p.value,
                                    x = ContrastType:Condition,
                                    colour = ContrastType:Condition)) +
  geom_violin(width=2) +
  geom_boxplot(width=.1) +
  geom_jitter(width = .1, height = 0, alpha = .2) +
  theme(legend.position = "top")
# Plot b-values
new_old.sims.bayesian.plot <- ggplot(new_old.sims.bayesian.evid,
                                     aes(y = Evid.Ratio,
                                         x = Hypothesis,
                                         colour = Hypothesis)) +
  geom_violin(width=2) +
  geom_boxplot(width=.1) +
  geom_jitter(width = .1, height = 0, alpha = .2) +
  theme(legend.position = "top")
