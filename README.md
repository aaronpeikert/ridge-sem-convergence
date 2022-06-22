# Simulation for Convergence of Ridge SEM

Benedikt Friemelt, Christian Bloszies & Tobias Koch found that `RegSEM` and `lslx` had many convergence issues when plain ML did not (for a particular class of models with high multicollinearity).
This result is very suprising since Ridge should be at least as well behaved for normal optimizers as ML (Lasso is another story).
We (Aaron Peikert, Maximilian S. Ernst, Andreas M. Brandmaier) hence repeated the simulation with `StructuralEquationsModell.jl` (without `Lasso` as this *was* beyond the ability of the package at that time).
Our implementation of Lasso differs notably from the others.
It is a minimal implementation of Ridge, with standard optimizers (L-BFGS) used for ML in `lavaan` and `OpenMX`.
`RegSEM` and `lslx`, on the other hand, use optimizers and algorithms that can handle the nondifferentiable Lasso.
Ridge is differentiable, so these are not quite necessary, and the hypothesis was that these add-ons lead to convergence issues.

Indeed we found that the simulation is considerably faster (2h on a laptop, so 200ms per condition (ca. 200) per repetition (200) or 18ms per single model (50 hyperparameters)) and had, at most, a handful of convergence problems (<5 for the 200000 models (200x200x50)).

This simulation uses an old experimental version of `StructuralEquationsModell.jl` with a custom interface to lavaan partables, so nothing to be used in production.

# Reproduction

The commands (paste into terminal):

```
git clone https://github.com/aaronpeikert/ridge-sem-convergence.git
cd ridge-sem-convergence
make
```

reproduces the whole project if the software requirements are met.
It might be that we overlooked to set a seed (not easy across different programming languages); however, if the seed would influence the results considerably, we would be worried.
The simulation results are saved to the file `simulation_results.csv` and can be analyzed further.

## Dependencies

This project needs `julia` and `R` and many packages (we are lazy, we know).
Exact versions of the used julia packages are recorded in `julia/StructuralEquationModels.jl/Manifest.toml`, and the old experimental version of `StructuralEquationsModell.jl` is saved as a complete source.
For R we dumped the `sessionInfo()` in the file `RsessionInfo.txt`.
