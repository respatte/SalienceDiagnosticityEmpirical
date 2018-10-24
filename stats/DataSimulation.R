# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(lmerTest)
library(emmeans)
library(brms)
library(coda)
library(tidyverse)

# BAYESIAN POINT NULL HYPOTHESIS TESTING ===========================================================
H_naught.test <- function(N){
  # Generate random data
  # n=42 * 2trials from real data, mean to 0, sd=.5 from real data
  new_old.sims <- replicate(25,
                            list(tibble(ChanceArcsin = rnorm(N*2, 0, .5),
                                        Participant = rep(1:N, each = 2),
                                        Condition = factor(rep(c("Label", "No Label"),
                                                               times = N/2, each = 2)),
                                        ContrastType = factor(rep(c("Head", "Tail"), times = N)))))

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
                           chains = 4, cores = 4, iter = 20000, warmup = 2000,
                           control = list(adapt_delta = .99, max_treedepth = 20),
                           save_all_pars = T, combine = F, sample_prior = "yes")
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
  new_old.sims.bayesian.evid <- bayesian.analysis(new_old.sims) %>%
    mutate(Condition = factor(ifelse(grepl("NoLabel", Hypothesis), "No Label", "Label")),
           ContrastType = factor(ifelse(grepl("Tail", Hypothesis), "Tail", "Head")))
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
                                           x = ContrastType:Condition,
                                           colour = ContrastType:Condition)) +
    geom_violin(width=2) +
    geom_boxplot(width=.1) +
    geom_jitter(width = .1, height = 0, alpha = .2) +
    theme(legend.position = "top")
  return(list(BayesianEvidence = new_old.sims.bayesian.evid,
              #BayesianPlot = new_old.sims.bayesian.plot,
              SampleTheoryEvidence = new_old.sims.stb.evid))
              #SampleTheoryPlot = new_old.sims.stb.plot))
}

#new_old.sims.results.42 <- H_naught.test(N = 42)
#new_old.sims.results.200 <- H_naught.test(N = 200)
new_old.sims.results.42 <- H_naught.test(N = 42)
l <- names(new_old.sims.results.42)
for(i in 1:39){
  tmp <- H_naught.test(N = 42)
  new_old.sims.results.42 <- Map(bind_rows, new_old.sims.results.42[l], tmp)
}
