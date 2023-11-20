Hpothesis for weird results in individual parameters
1) Need to logtransform the param space?
2) due to upper and lower bounds of model?

1) try indivudal
2) check nuts initial leap

8473263        Short           fdc3_n10_cauchy      Queued   no start time estimate yet
8473265        Short           fdc3_n1_cauchy       Queued   no start time estimate yet

Do fake data check while comparing following paramters:
- 1 vs 3 parameters
- lognormal vs normal 
- informative vs non informative priors
- cauchy vs normal priors
- noise scale