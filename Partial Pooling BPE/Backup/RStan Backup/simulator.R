require(deSolve)
source("treatment_doser.R")
require(colorout)

simulate_traj <- function() {
  ##-----------------------------
  ## the derivative function
  ##-----------------------------
  derivs  <-  function(t, y, params) {
    if (t < td)
      lag  <-  0
    else 
      lag  <-  lagvalue(t - td)

    d_12  <- 0; d_cpi <- 0
    d_cbd <- treatment_doser(t, c(7), t_delay, t_last)

    dy1  <-  k1 + k2*(d_cbd + d_12) - d1 * y[1]
    dy2  <-  k3 + k4*lag[1] - d2*y[2]
    dy3  <-  k5 - (d3 + d4*y[1])*y[3]
    dy4  <-  k6*(1-(y[4]+y[5])/v_max)*y[4] - (d5 + (d6*y[2]/(1+s1*y[3]*(1-d_cpi)) + d7*y[1])/(1+s2*(y[4]+y[5])))*y[4]
    dy5  <-  (d5 + (d6*y[2]/(1+s1*y[3]*(1-d_cpi)) + d7*y[1])/(1+s2*(y[4]+y[5])))*y[4] - d8*y[5]

    list(c(dy1, dy2, dy3, dy4, dy5))
  }

  ## ----------------------------
  ##  initial values and times
  ## ------------------------
  yinit  <-  c(0.0084, 9.56, 4.96, 6.65, 0)
  times  <- seq(0, 30)

  # Solver
  yout  <-  dede(y = yinit, times = times, func = derivs, parms = NULL)

  return(yout)
}
