# FAI Capstone: Reproducing GAVI for CSPs

This repository contains the implementation for my CS5100 Foundations of 
Artificial Intelligence final project, a reproduction of 
Kanoh, Matsumoto, Hasegawa, Kato and Nishihara (1997), *Solving Constraint 
Satisfaction Problems by a Genetic Algorithm Adopting Viral Infection*, 
Engineering Applications of Artificial Intelligence 10(6), 531 to 537.

## What this project does

The project reproduces the paper's core claim that replacing mutation in 
a genetic algorithm with *viral infection* of partial solutions produces 
faster and more reliable search on binary constraint satisfaction problems. 
Three solvers are implemented and compared head to head on the same 
randomly generated instances:

- **IHC**, iterated hillclimbing with the min-conflicts heuristic and 
  random restarts, matching the GSAT-style baseline used in the paper.
- **GA**, a standard genetic algorithm with roulette selection, 
  uniform crossover and 2% per-locus mutation.
- **GAVI**, the paper's proposed method: the same GA scaffolding but 
  with mutation replaced by an infection operator that injects 
  fitness-weighted partial solutions (viruses) into individuals.

## No dataset required

The project **does not use any external dataset**. CSP instances are 
generated on the fly inside the script itself, so you do not need to 
download anything. Given a target number of variables `n`, domain size 
`|L|`, constraint density `d` and mean partial solution set size 
`R_mean`, the generator:

1. Samples a random hidden solution.
2. Samples `m = d * n` distinct variable pairs as constraint scopes.
3. For each constraint, builds an allowed value pair set of size near 
   `R_mean`, with the hidden solution's induced pair always included 
   so the instance is guaranteed to be satisfiable.

This matches the experimental setup of the original paper and means the 
repository is fully self-contained.

## Repository contents

- **`FAI-capstone-GAVI.py`**  The only source file. It contains:
  - A `CSP` class and the random CSP generator `generate_csp`.
  - Shared utilities: `fitness`, `num_conflicts`, `is_solution`, 
    `roulette`, `uniform_cx`.
  - Implementations of the three solvers: `ihc`, `ga`, `gavi`.
  - `run_experiment`, which runs all three solvers on the same set of 
    generated instances and aggregates success rate and mean solution 
    time per algorithm.
  - A `__main__` block that sweeps the full grid of 5 densities x 3 
    R_mean values x 3 instance counts (10, 30 and 50 CSPs per cell), 
    prints the per-cell results, and prints two comparative tables 
    showing how the estimates stabilise as the number of CSPs per cell 
    grows.

## How to run

No external libraries are needed, only the Python standard library.

```bash
python FAI-capstone-GAVI.py
```

The script will sweep all configurations and print results to the 
console. A full run with 50 CSPs per cell takes on the order of a 
few hours on a modern laptop; to do a quick sanity check, edit 
`CSP_COUNTS` at the bottom of the file to `[10]` first.

## Experimental parameters

| Parameter | Value |
|---|---|
| Number of variables `n` | 50 |
| Domain size `|L|` | 4 |
| Constraint density `d` | 0.98, 2, 3, 4, 5 |
| R_mean | 2, 4, 6 |
| Instances per cell | 10, 30, 50 |
| GA / GAVI population | 500 |
| GA / GAVI max generations | 200 |
| GA mutation rate | 2% per locus |
| GAVI infection rate | 40% |
| IHC hillclimb steps / restarts | 500 / 200 |

All of these are defined at the bottom of `FAI-capstone-GAVI.py` and can 
be changed directly in the script.

## Output

For every combination of `N_CSPS`, `R_mean` and `d`, the script prints 
success rate and mean solution time for IHC, GA and GAVI. After the 
sweep finishes, it also prints two comparative tables showing how 
success rate and mean solve time evolve across the three instance 
counts, which is useful for checking that the reported estimates are 
stable.

## Report

A written report accompanying this implementation is submitted 
alongside the repository. See the final report PDF for the full 
description of the method, experiments, analysis and reproducibility 
verdict.
