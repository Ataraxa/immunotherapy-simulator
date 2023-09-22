source("simulator.R")
require(FamilyRank)
require(MASS)
set.seed(123)

# Settings of the samples generator
NUM_TRAJ <- 10

for (traj in 1:NUM_TRAJ) {
  # Sample parameter values using hyperparameters
  µ1 <- 1; µ2 <- 1; s1 <- 1; s2 <- 1; prop <- 0.5
  k6 <- rbinorm(1, µ1, µ2, s1, s2, prop)

  µ1 <- 1; µ2 <- 1; s1 <- 1; s2 <- 1; prop <- 0.5
  d1 <- rbinorm(1, µ1, µ2, s1, s2, prop)

  µ1 <- 1; µ2 <- 1; s1 <- 1; s2 <- 1; prop <- 0.5
  s2 <- rbinorm(1, µ1, µ2, s1, s2, prop)

  # Default parameters
  td <- 1.556
  t_delay  <-  0.4151
  t_last  <-  5.2695
  t_delay12  <-  3.9
  t_last12  <-  1.3903

  k1 <- 0.19; k2 <- 6.04; k3 <- 79.566
  k4 <- 1054; k5 <- 5.44 

  d2 <- 9.72; d3 <- 1.26
  d4 <- 4.85; d5 <- 0.02; d6 <- 0.04
  d7 <- 51.4; d8 <- 0.56

  s1 <- 14.5 
  v_max <- 600

  # Simulate trajectory
  simulation <- simulate_traj()
  temp <- simulation[4,]+simulation[5,]
  plot(simulation[4,], type = "l", lwd = 2, col='blue', main = "test")
  lines(simulation[5,], type = "l", lwd = 2, col='orange', main = "test")
  lines(temp, type = "l", lwd = 2, col='red', main = "test")

  # Save trajectories 
  write.matrix(t(simulation[,]), file=sprintf("Data/traj-%d.csv", traj), sep=",")
  
  # Save parameters
  params_mat <- matrix(c(k6, d1, s2), byrow=TRUE)
  dim(params_mat) <- c(1,3)
  colnames(params_mat) <- c("k6", "d1", "s2")
  write.matrix(params_mat, file=sprintf("Data/params-%d.csv", traj), sep=",")
}

