# barro colorado example

library(entropart)
library(dbmss)
library(vegan)
library(tidyverse)
# for progress bar
library(cli)

source("r/functions.R")


# speceis names generation ------------------------------------------------

# LETTERS is the vector of capital letters A:Z
# outer() combines the first two args by FUN
outer(LETTERS, LETTERS, FUN = paste0) %>% 
  # add a third letter
  outer(LETTERS, paste0) ->
  species_names

# create the metacommunity ------------------------------------------------

mc_size <- 1E9
# theta is also Fisher's alpha
alpha <- 50
# Number of species (Fisher, 1943)
n_species <- -alpha * log(alpha / (mc_size + alpha)) %>% 
  # integer needed
  as.integer()


# calculate abundance -----------------------------------------------------

# Prepare a progress bar
if (interactive()) cli_progress_bar("Log-series", total = n_species)
# And a vector of abundances
mc_abundances <- numeric(n_species)
# Draw each species abundance
for (s in 1:n_species) {
  mc_abundances[s] <- LSabundance(mc_size, alpha)
  if (interactive()) cli_progress_update()
}

# Set the species names
names(mc_abundances) <- species_names[seq_along(mc_abundances)]

# Declare it is a vector of abundances (package entropart)
mc_abundances %>% 
  as.AbdVector() ->
  mc_abundances

# Range
summary(mc_abundances)
# Rank-Abundance Curve
autoplot(mc_abundances, Distribution = "lseries")


# draw a local community sample -------------------------------------------

# Choose the size of the community
c_size <- 256
# Draw a random community according to the abundances of the metacommunity
c_abundances <- rmultinom(1, size = c_size, prob = mc_abundances) %>% 
  # We need a vector, not a matrix
  as.vector()
# Give the species a name
names(c_abundances) <- names(mc_abundances)

# Eliminate species with 0 individual
c_abundances[c_abundances > 0] %>% 
  # Declare a vector of abundances
  as.AbdVector() ->
  c_abundances
# Show the RAC
autoplot(c_abundances, Distribution = "lnorm")


# demographic drift -------------------------------------------------------

# Build a dataframe for the marks of a point pattern
marks <- data.frame(
  # Species. Each line of the dataframe is an individual
  PointType = rep(names(c_abundances), times = c_abundances), 
  # Size of the points
  PointWeight = rep(1, c_size)
)
# Draw random locations in a square window to obtain a point pattern
c_size %>% 
  runifpoint() ->
  c_ppp
# Add the marks to the point pattern
marks(c_ppp) <- marks
# Make it a wmppp object (package dbmss) to make nice plots
c_wmppp <- as.wmppp(c_ppp)
# Set the factor levels of the point types identical to those of the metacommunity
# so that species of the metacommunity can be used in the local community
c_wmppp$marks$PointType <- factor(
  c_wmppp$marks$PointType, 
  levels = names(mc_abundances)
)
# Map the point pattern
autoplot(c_wmppp)

# test evolve_drift function
c_wmppp %>% 
  evolve_drift(show = TRUE)

# kill all the species until it reduces to one species

# Copy the community to modify it
the_community <- c_wmppp
# Number of species
richness <- length(unique(the_community$marks$PointType))
# Prepare a long vector to save the richness along time
richness_time <- numeric(1E6)
# Count the steps
step <- 0
# Prepare a progress bar
if (interactive()) cli_progress_bar("Drifting")
while (richness > 1) {
  step <- step + 1
  # Replace a tree
  the_community <- evolve_drift(the_community)
  # Update the richness
  richness <- length(unique(the_community$marks$PointType))
  # Save it
  richness_time[step] <- richness
  # Show the progress
  if (interactive()) cli_progress_update()
}
# Finalise the progress bar
if (interactive()) cli_progress_update(force = TRUE)
# Save the necessary time to reach a single species
drift_time <- step
# Plot the community
autoplot(the_community)
# Plot the decreasing richness
data.frame(Time = seq_len(step), Richness = richness_time[seq_len(step)]) %>% 
  ggplot() +
  geom_line(aes(x = Time, y = Richness))


# imigration drift --------------------------------------------------------

# Choose the migration probability
m <- 0.5

# Copy the community to modify it
the_community <- c_wmppp
# Number of species
richness <- length(unique(the_community$marks$PointType))
# Prepare a long vector to save the richness along time
richness_time <- numeric(1E6)
# Prepare a progress bar
if (interactive()) cli_progress_bar(paste("Migration", m, ":"), total = drift_time)
for (step in seq_len(drift_time)) {
  # Replace a tree
  the_community <- evolve_migrate(the_community, mc_abundances, m)
  # Save the richness
  richness_time[step] <- length(unique(the_community$marks$PointType))
  # Show the progress
  if (interactive()) cli_progress_update()
}
autoplot(the_community)
# Plot richness
data.frame(Time = seq_len(drift_time), Richness = richness_time[seq_len(drift_time)]) %>% 
  ggplot() +
  geom_line(aes(x = Time, y = Richness))
# Plot the rank-abundance curve
the_community %>% 
  # Consider the dataframe of marks
  marks() %>% 
  # Group by species
  group_by(PointType) %>% 
  # Count the numnber of individuals by species
  summarise(Abundance = n()) %>% 
  # Extract the column 
  pull(Abundance) %>% 
  # Make it an abundance vector (package entropart)
  as.AbdVector() %>% 
  # Plot the rank-abundance curve and try to fit a log-normal distribution
  autoplot(Distribution = "lnorm")
