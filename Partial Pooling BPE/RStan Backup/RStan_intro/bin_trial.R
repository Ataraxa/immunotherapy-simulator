bern.stan =
"
data {
  int<lower=0> N;               // number of trials
  int<lower=0, upper=1> y[N];   // success on trial n
}

parameters {
  real<lower=0, upper=1> theta; // chance of success
}

model {
  theta ~ uniform(0, 1);        // prior
  y ~ bernoulli(theta);         // likelihood
}
"

theta = 0.3
N = 20
y = rbinom(N, 1, 0.3)

library(rstan)

fit = stan(model_code=bern.stan, data=list(y=y, N=N), iter=5000)

print(fit, probs=c(0.1, 0.9))

theta_draws = extract(fit)$theta

mean(theta_draws)
