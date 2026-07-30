# Examples

The standalone Julia scripts in this directory demonstrate grids,
fields, time stepping, diffusion, KPP, TKE closures, and diagnostic
plumes. The two subdirectories below provide more complete,
paper-based examples:

- [`free_convection_les/`](free_convection_les/):
  compares
  OceanTurb closures with the official LES profiles underlying Figure
  5 of Wagner et al. (2025).
- [`rotating_strong_wind_les/`](rotating_strong_wind_les/):
  compares the closures with the official rotating strong-wind LES
  profiles underlying Figure 7.

These are physical-case validation examples rather than exact
reproductions of every manuscript curve. They compare CATKE, KPP,
Mellor--Yamada 2.5, k--omega, and k--epsilon from common LES initial
conditions. The rotating strong-wind LES includes surface-wave effects
that are not yet represented by the OceanTurb columns.

Each directory contains:

- compact, preprocessed LES profiles needed by the runner;
- a Julia runner for all five closures;
- an optional preprocessing script documenting how the LES profiles
  were produced from the official archive;
- a plotting notebook with outputs cleared; and
- one representative comparison figure.

Large source archives and generated model output are intentionally not
version controlled. See each case README for commands, data provenance,
and interpretation.
