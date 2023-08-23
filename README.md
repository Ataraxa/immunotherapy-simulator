This repository contains all the code and data used for Cancer Immunotherapy UROP project.

# How to use the repo
There is one folder for each major component of the project. 

## Complete Poolig GA
Code and tools used to perform parameter fitting through a Genetic Algorithm. The parameters of the ODE model were fitted to the average tumour volume data. This design choice is questionnable and discussed in the various report. (Please note that the code in mainly written in MATLAB)

## Partial Pooling BPE
Code to implement Bayesian Inference on the individual data. A hierarchical model was used. Written in RSTAN.

## Data
Folder that contains all the data, and some data processing tools to convert the data into more usesable formats.

## Data Visualisation
Folder that contains various code to plot the processed data.