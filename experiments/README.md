# Paper experiments

This directory holds scripts to replicate the results reported in the library
accompanying papers.

## Before running experiments

Build the library locally, or open the GitPod online environment.

## Forward dynamics

From within this `MBSE/experiments` directory, run:

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


## Inverse dynamics

From the `MBSE/experiments` directory, run:

```bash
../build/bin/mbse-fg-inverse-dynamics   \
  --mechanism ../config/mechanisms/fourbars1-with-rel-angle.yaml \
  --desired-trajectory ../config/trajectories/fourbars1-with-rel-angle-trajectory.txt \
  --imposed-coordinates "[ 4 ]"
  # --verbose

../build/bin/mbse-fg-inverse-dynamics  \
  --mechanism ../config/mechanisms/pick-and-place-robot.yaml \
  --desired-trajectory ../config/trajectories/pick-and-place-robot-trajectory.txt \
  --imposed-coordinates "[ 20 ; 21 ]"
  #--verbose
```  

And see results with:

```octave
plot_inverse_dynamics_results_4bars
```
