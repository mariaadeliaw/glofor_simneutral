#' Stochastic abundance of a (single) species in a log-series distribution
#'
#' @param N The size of the metacommunity
#' @param alpha Fisher's alpha, that links the number od species to that of individuals
#'
#' @return The number of individuals of a species
LSabundance <- function(N, alpha) {
  # adapted from Dan Lunn, http://www.stats.ox.ac.uk/~dlunn/BS1_05/BS1_Rcode.pdf
  # Fisher's x is log-series 1-theta
  x <- N / (N + alpha)
  # Draw a random number between 0 and 1
  u <- stats::runif(1)
  # k is the number of individuals to draw
  k <- 1
  # Calculate the probability at k=1
  P <- -x / log(1 - x)
  # Store it in the distribution function
  F <- P
  # Repeat while the cumulated probability is below u
  while (u >= F) {
    # Probability at k+1 obtained from that at k
    P <- P * k * x / (k + 1)
    # Increment k
    k <- k + 1
    # Increment the cumulated probability
    F <- F + P
  }
  return(k)
}

#' One-step evolution of a community without migration
#'
#' @param community A wmppp object (package dbmss) where points are trees.
#' @param show If TRUE, a map of the evolved community is produced.
#'
#' @return A wmppp object with the modified community.
evolve_drift <- function(community, show = FALSE) {
  # Choose a tree to be killed
  focus_tree <- sample(seq_len(community$n), size = 1)
  # Species of the survivors
  survivors_sp <- community$marks$PointType[-focus_tree]
  # Draw the species of the new tree
  community$marks$PointType[focus_tree] <- 
    survivors_sp[sample(seq_along(survivors_sp), size = 1)]
  if (show) {
    # Copy the community
    community_to_plot <- community
    # Increase the size of the modified tree
    community_to_plot$marks$PointWeight[focus_tree] <- 2
    # Plot the community (print is mandatory inside a function)
    print(autoplot(community_to_plot))
  }
  # return
  return(community)
}

#' One-step evolution of a community with migration
#'
#' @param community A wmppp object (package dbmss) where points are trees.
#' @param metacommunity A named vector of abundances or probability that describes the metacommunity.
#' @param m The migration parameter, that is the probability for a dead tree to be replaced by an immigrant species.
#'
#' @return A wmppp object with the modified community.
evolve_migrate <- function(community, metacommunity, m = 0) {
  # Choose a tree
  focus_tree <- sample(seq_len(community$n), size = 1)
  # Migration or local parent? Draw in a uniform distribution.
  u <- runif(1)
  if (u > m) {
    # Local parent: Species of the survivors
    survivors_sp <- community$marks$PointType[-focus_tree]
    # Draw the species of the new tree
    community$marks$PointType[focus_tree] <- 
      survivors_sp[sample(seq_along(survivors_sp), size = 1)]
  } else {
    # Migration : Draw the species of the new tree in the metacommunity
    community$marks$PointType[focus_tree] <- 
      names(metacommunity)[which(rmultinom(1, size = 1, prob = metacommunity) == 1)]
  }
  # return
  return(community)
}