# DAMASK 3.1.0 with extensions

This branch is DAMASK 3.1.0 (tag `v3.1.0`) plus the changes below. Everything else is upstream DAMASK; see `README.md` and https://damask-multiphysics.org for documentation, and the upstream licence (AGPL-3.0-or-later) for terms.

## Additions

| Feature | Files | Configuration |
|---|---|---|
| Dislocation-density plasticity with twinning for hP lattices | `src/phase_mechanical_plastic_dislotwinhcp.f90` | `plastic: {type: dislotwinhcp, ...}` |
| Isotropic ductile damage source (DAMASK 2.0.3 `isoDuctile`) | `src/phase_damage_isoductile.f90` | `damage: {type: isoductile, gamma_crit: ..., N: ...}` |
| Local (per-cell) damage solver, no gradient term | `src/grid/grid_damage_local.f90`, `src/grid/DAMASK_grid.f90` | load case `solver: {damage: local}` |

The local solver reproduces DAMASK 2.0.3's local damage homogenisation and is mesh dependent. For predictive damage studies use the stock nonlocal solver (`damage: spectral`) with a physical length scale and a mesh-sensitivity check.

## Fixes

- **phenopowerlaw twinning:** the untwinned-fraction prefactor is bounded at zero, so a point can no longer exceed 100 % twinned and un-twin. Two optional parameters, `f_sat_tw` (attainable twin fraction, default 1) and `m_tw` (exponent, default 1), reproduce the released model at their defaults while the twin fraction stays at or below 1.
- **material.yaml with mixed constituent counts:** stock 3.1.0 crashes at start-up when homogenizations with different `N_constituents` are used in one file. Backport of upstream f19b1138e9 (`material.f90`, `result.f90`, `python/damask/_result.py`); constituents a cell does not have are written with an empty phase label in `cell_to/phase`. The DADF5 format version is unchanged.

## Build

Build as stock DAMASK 3.1 (CMake with `-DGRID=ON`); the new sources are picked up automatically. Requirements are those of DAMASK 3.1: PETSc 3.19 to 3.25, HDF5 with Fortran bindings, FFTW, libfyaml and zlib.
