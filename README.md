# FAI-capstone-GAVI
Reproduction of Kanoh et al. 1997, Solving CSPs by a Genetic Algorithm Adopting Viral Infection. CS5100 Spring 2026.

# FAI Capstone: Reproducing GAVI for CSPs

This repository contains the implementation for my CS5100 final project, 
reproducing "Solving Constraint Satisfaction Problems by a Genetic Algorithm 
Adopting Viral Infection" (Kanoh et al., 1997).

## Contents
- `FAI-capstone-GAVI.py`  — Implementation of IHC, GA and GAVI, along with the 
  random CSP generator and the experimental loop.

Running the script sweeps over d in {0.98, 2, 3, 4, 5} and R_mean in {2, 4, 6} 
at n = 50, |L| = 4, and prints the success rate and mean solution time 
for each of IHC, GA and GAVI.

## Parameters
- Population size: 500
- Max generations: 200
- GA mutation rate: 2%
- GAVI infection rate: 40%
- Instances per cell: configurable at the top of the script

## Report
See the final report PDF submitted alongside this repository.
