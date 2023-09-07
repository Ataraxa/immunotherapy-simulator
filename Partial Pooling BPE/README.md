<h1>Bayesian Parameter Estimation</h1>

This folder contains all the code necessary to perfrom a Bayesian inference on the data. The full theoretical procedure is outlined in the <b>Documents</b> folder.

<h2>Backup</h2>
Contains legacy functional code that is now outdated

<h2>Data</h2>
Contains sample (generated) data for validation as well as the true experimental data for the final inference.

<h2>Model</h2>
Where most of the mathematical code is. Contains the DDE model, the statistical model as well as a custom statistical distribution implemented in Turing.jl - the Binormal distribution.

<h2>Res</h2>
Folder that stores the MCMC chains calculated by the HPC cluster. t acts as an interface between the local machien and the remote servers, to share the results.