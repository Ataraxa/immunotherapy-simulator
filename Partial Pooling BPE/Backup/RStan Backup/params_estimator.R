# This is the main file where the parameter estimation is done (by calling the STAN model)
library(rstan)

fit <- stan(model_code=model.stan, data=list(), iter=5000)

