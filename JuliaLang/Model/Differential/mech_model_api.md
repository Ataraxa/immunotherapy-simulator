# Mechanistic Model - API Guide

## Folder Structure
There are three files to make the API work.
- `ode_core.jl`: contains the function `create_problem`. This function returns a parameterised model that can be solved using the `solve` function.

- `ode_restricted.jl`: contains an auxiliary function used to convert a restricted-space parameter vector ($\vec{\theta}\in\mathbb{R}^3$ for example) into full-space vector compatible with the `ode_core` set of functions.

- `ode_params`: self explanatory.