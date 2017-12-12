library(eyetrackingR)
library(Matrix)
library(lme4)
library(tidyverse)

source("Routines.R")

# GATHER DATA
# Define AOIs
# AOIs.plot <- ggplot(NULL, aes(xmin = L, xmax = R, ymin = T, ymax = B)) +
#   xlim(c(0,1920)) + scale_y_reverse(limits = c(1080,0)) +
#   geom_rect(data = AOIs.infants.head, aes(fill = AOI_type)) +
#   geom_rect(data = AOIs.infants.tail, aes(fill = AOI_type))
# ggsave("../results/infants/data_cleaning_graphs/AOIs.png",
#        plot = AOIs.plot , width = 6.05, height = 2.95)

d <- LT_data.gather("infants")