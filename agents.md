# agents.md — context for AI coding agents

This file gives AI agents the context needed to work effectively in this
repository. Keep it up to date (see "Keeping this file in sync" below).

## What this project is

`libmbse` (MultiBody State Estimation) is a C++ library for multibody
forward/inverse dynamics, state estimation, and parameter identification.
Two paradigms are implemented:

1. **Factor-graph based** (GTSAM): kinematic/dynamic constraints of
   multibody systems expressed as factors in a factor graph. This is the
   actively developed paradigm.
2. **Particle filter estimator** (`MultiBodyParticleFilter`): a Bayesian
   filter, tested only on planar dynamics models.

Related theory is in the papers linked in `README.md`.

## Repository layout

- `libmbse/include/mbse/`, `libmbse/src/` — the library itself.
  - `ModelDefinition.h/.cpp` — user-facing description of a mechanism
    (points, bodies, constraints, relative DOFs). Can be loaded from YAML
    (`ModelDefinition::FromYAML`, format documented in
    `docs/mechanism-definition-yaml.md`).
  - `AssembledRigidModel.*` — the assembled/preprocessed model used to run
    kinematic/dynamic simulations (mass matrix, generalized forces, initial
    position problem).
  - `constraints/` — one class per joint/constraint type (constant
    distance, fixed/mobile slider, relative angle, relative position, etc.),
    all deriving from `ConstraintBase`.
  - `dynamics/` — dynamic simulators (dense/sparse, several solvers:
    Lagrange multipliers via CHOLMOD/KLU/UMFPACK/dense LU, augmented
    Lagrangian, independent-coordinates formulations).
  - `factors/` — GTSAM factors implementing the factor-graph paradigm
    (dynamics, constraints — position/velocity/acceleration,
    independent-coordinate variants, integrators, sensor factors).
- `apps/` — example/demo executables, each its own CMake subproject:
  `mbse-dynamic-simulation`, `mbse-viewer`, `mbse-fg-smoother-forward-dynamics`
  (+ `-icoords` variant), `mbse-fg-inverse-dynamics`, `mbse-pf-demo`.
- `config/mechanisms/*.yaml` — example mechanism definitions consumed by the
  apps and tests. `config/trajectories/` — reference trajectory data.
- `libmbse/tests/` — GoogleTest-based unit tests (bundled gtest 1.6.0
  sources; treated as external code, excluded from clang-format checks).
- `docs/` — Sphinx/Doxygen documentation sources (published to
  readthedocs.io); `doc/` is generated Doxygen output, not source.
- `experiments/` — MATLAB/Octave scripts and reference data used for the
  factor-graph paper's experiments; not part of the build.

## Build

```
mkdir build && cd build
cmake ..
make                # library + apps
ctest --verbose     # or: make test_legacy
```

Dependencies: CMake, a C++ compiler, SuiteSparse, MRPT (>=2.0.0), GTSAM
(>=4.2), GTSAM_UNSTABLE.

## GTSAM version branches

- `master` branch targets GTSAM <=4.2.
- Branch `newer-gtsam` targets GTSAM >=4.3 (or GTSAM's `develop` branch).

Check which branch you're on before assuming GTSAM API availability.

## Code style

- Formatted with clang-format (config in `.clang-format`) using
  `clang-format-14`; run `./formatter.sh` to reformat in place, or
  `./formatter.sh --check` to verify. `libmbse/tests/gtest-1.6.0/` is
  external/vendored and excluded from formatting.
- No one-line `if`/`for` bodies — always brace on its own line.
- One variable declaration per line.
- No em/en dashes in code or comments; American spelling.
- Use anonymous namespaces instead of `static` for internal linkage.
- Comments should give short, generic reasoning, not "fixes bug X in
  dataset Y" style detail.

## CI

- `.github/workflows/linux-build.yml`: builds + runs `ctest` on a matrix of
  Ubuntu 20.04/22.04 x gcc/clang, with and without the MRPT PPA.
- `.github/workflows/check-clang-format.yml`: runs `formatter.sh --check`.

## Keeping this file in sync

When you make a **significant** change to the project (new module/paradigm,
new app, changed build/dependency requirements, changed directory layout,
new supported GTSAM/MRPT version range, etc.), update this file accordingly.
Routine bug fixes, small refactors, or one-off experiment tweaks do not need
an entry here.
