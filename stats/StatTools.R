# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(doSNOW)
library(tidyverse)

lmer_test.SS2.par <- function(formula, data){
  terms <- tibble(label = attr(terms(f), "term.labels"),
                  order = attr(terms(f), "order"))
}
