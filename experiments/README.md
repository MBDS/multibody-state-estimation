# Paper experiments

This directory holds scripts to replicate the results reported in the library
accompanying papers.

## Instructions

- Build the library locally, or open the GitPod online environment.

- From within this `MBSE/experiments` directory, run:

```bash
./run-experiments-forward-dynamics-dcoord.sh
./run-experiments-forward-dynamics-icoord.sh
```

- Then, from MATLAB/Octave, generate the graphics with:

```octave
% See results for dependent coordinates FG forward dynamics:
plots_forward_dynamics_comparison('bars4_dc');

% See results for independent coordinates FG forward dynamics:
plots_forward_dynamics_comparison('bars4_ic');
```
