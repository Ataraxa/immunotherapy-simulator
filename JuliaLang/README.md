<h1>Welcome to the folder JuliaLang</h1>

<h2>Analysis Scr</h2>
This folder contains all the code necessary to perfrom a Bayesian inference on the data. The full theoretical procedure is outlined in the <b>Documents</b> folder.

The general flow is as follows: the `main_inference.jl` is called, which is the heart of the code. It combines a `DDEProblem` object with the data vector by passing them to the Bayesian model. The resulting fitted model is sampled by a MCMC chain, stored in a `.h5` file for further analysis.

<h2>Backup & RStan Backup</h2>
Contains legacy functional code that is now outdated.

<h2>Data</h2>
Contains sample (generated) data for validation as well as the true experimental data for the final inference.

<h2>Model</h2>
Where most of the mathematical code is. Contains the DDE model, the statistical model as well as a custom statistical distribution implemented in Turing.jl - the Binormal distribution.

<h2>Res</h2>
Folder that stores the MCMC chains calculated by the HPC cluster. It acts as an interface between the local machine and the remote servers, to share the results.

<h2>Test</h2>
Collection of tests to check that each individual component of the project behaves as expected. 

<h1>Useful commnand<h1>
`; tput rmam` to artificially reduce stacktrace