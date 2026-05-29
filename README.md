# Dormancy stabilizes non-transitive competitive dynamics

This repository contains the simulation code used in the manuscript:

> *Dormancy stabilizes non-transitive competitive dynamics*

## Repository Contents

### `SDELIMITSIMULATION.R`

Contains the simulations used for the *Rock--Paper--Scissors* and *Tournaments* sections.

Main functions include:

* `simulate_SDE`: Euler--Maruyama simulation of the SDE in Equation (SDEXY)
* `createSigmaOrg`: constructs the function (\zeta)
* `createSigmaAlt`: numerically stable implementation with equivalent law
* `simulate_SDE_Times_rcpp` and `run_SDE_Times_vec`: computation of fixation times under varying parameter configurations

### `SDELIMITSIMULATION.cpp`

C++ implementation of the Euler--Maruyama scheme used for faster computation of fixation times.

### `ACSIMULATION.R`

Contains the Euler--Maruyama implementation used for the *Majority Voting* section and the SDE in Equation (AC).


## Citation

If you use this code, please cite the associated PNAS article.
