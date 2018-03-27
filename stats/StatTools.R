# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(doSNOW)
library(tidyverse)

lmer_test.SS2.par <- function(formula, data = NULL){
  # Generate all formulae
  terms <- tibble(label = attr(terms(formula), "term.labels"),
                  order = attr(terms(formula), "order")) %>%
    subset(!grepl("\\|", label))
  n.order <- max(terms$order)
  formulae <- list("Intercept" = update.formula(formula,
                                                str_c(". ~ . - ",
                                                      str_c(terms$label, collapse = " - "))))
  # for(i in 1:n.order){
  #   n.order_i <- sum(terms$order == i)
  #   order_i <- terms %>% subset(order <= i)
  #   for(j in 1:n.order_i){
  # 
  #   }
  # }
}
