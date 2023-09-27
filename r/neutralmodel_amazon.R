# amazon case https://www.science.org/doi/10.1126/science.1243092 

# load necessary packages

library(entropart)
library(dbmss)
library(vegan)
library(tidyverse)
# for progress bar
library(cli)

source("r/functions.R")

# create metacommunity ----------------------------------------------------

# species name code generation (combination of 26^3 -> more than the number of species in amazonian number of species)

species_names <- outer(LETTERS, LETTERS, FUN = paste0) %>% 
  # add a third letter
  outer(LETTERS, paste0)

# metacommunity size
# number of trees and the alpha value in the amzonian (ter steege 2013)
mc_size <- 3.9E11 #number of individuals
alpha <- 754
# Number of species (Fisher, 1943)
n_species <- -alpha * log(alpha / (mc_size + alpha)) %>% 
  # integer needed
  as.integer()

# load the preloaded result
if (file.exists("data/amazonia.Rda")) {
  load("data/amazonia.Rda")
} else {
  # Apply LSabundance to each species
  #  replicate(n_species, expr = LSabundance(mc_size, alpha)) does not have a progress bar
  # Loop instead
  # Prepare a progress bar
  if (interactive()) cli_progress_bar("Log-series", total = n_species)
  # And a vector of abundances
  amz_abundances <- numeric(n_species)
  # Draw each species abundance
  for (s in 1:n_species) {
    amz_abundances[s] <- LSabundance(mc_size, alpha)
    if (interactive()) cli_progress_update()
  }
  save(amz_abundances, file="data/amazonia.Rda")
}

# Set the species names
names(amz_abundances) <- species_names[seq_along(amz_abundances)]

# Declare it is a vector of abundances (package entropart)
amz_abundances %>% 
  as.AbdVector() ->
  amz_abundances

# Range
summary(amz_abundances)
# Rank-Abundance Curve
autoplot(amz_abundances, distribution = "lseries")
