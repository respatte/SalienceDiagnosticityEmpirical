# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(doSNOW)
library(tidyverse)

lmer_test.SS2.par <- function(formula, data){
  terms <- tibble(label = attr(terms(f), "term.labels"),
                  order = attr(terms(f), "order"))
  n.order <- max(terms$order)
  for(i in 1:n.order){
    n.order_i <- sum(terms$order == i)
    for(j in 1:n.order_i){
      #TODO
    }
  }
}
