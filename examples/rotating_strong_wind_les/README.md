# Rotating strong-wind LES example

This experiment compares OceanTurb CATKE, KPP, Mellor--Yamada 2.5,
k--omega, and k--epsilon with rotating strong-wind LES profiles. These
cases underlie Figure 7 of Wagner et al. (2025). Every closure starts
from the same official LES velocity, buoyancy, and passive-tracer
profiles at 600 seconds.

This is a physical-case comparison, not an exact reproduction of every
manuscript curve. The manuscript used CATKE, KPP, and SMC-LT, whereas
this workflow compares five OceanTurb closures. The official
wind-driven LES also includes surface-wave/Stokes effects. OceanTurb
does not yet include explicit Langmuir or Stokes-drift enhancement, so
LES--OceanTurb differences include both closure behavior and that
missing physics.

## Configuration

| Setting | Value |
|---|---:|
| Domain depth | 256 m |
| OceanTurb cells | 128 |
| Start time | 600 s |
| Default time step | 10 s |
| Final times | 6, 12, 24, 48, 72 h |
| Surface buoyancy flux | 0 |
| Coriolis parameter | \(10^{-4}\ {\rm s^{-1}}\) |
| Deep buoyancy gradient | \(2\times10^{-6}\ {\rm s^{-2}}\) |
| Initial TKE for prognostic closures | \(10^{-9}\ {\rm m^2\,s^{-2}}\) |

The five surface \(u\)-momentum fluxes are
\(-1.4\times10^{-3}\), \(-9.0\times10^{-4}\),
\(-6.8\times10^{-4}\), \(-4.5\times10^{-4}\), and
\(-4.1\times10^{-4}\ {\rm m^2\,s^{-2}}\), respectively.

## Run

The processed LES inputs required by the runner are committed in
`data/`. From the OceanTurb repository root, run:

```bash
julia --project=. examples/rotating_strong_wind_les/run.jl
```

For a shorter single-suite check:

```bash
OCEANTURB_SUITES=6 julia --project=. \
  examples/rotating_strong_wind_les/run.jl
```

Set `OCEANTURB_STRONG_WIND_DT` to override the default 10-second time
step.

The runner writes profiles and a summary to `output/dt10s/`. Generated
output is ignored by Git.

To regenerate the comparison figures after running the cases:

```bash
cd examples/rotating_strong_wind_les
jupyter execute --inplace --timeout=600 plot_profiles.ipynb
```

The committed representative result is
`figures/rotating_strong_wind_closures.png`.

## Data provenance

- CATKE manuscript: <https://glwagner.github.io/assets/pdf/CATKE.pdf>
- Archived comparison code, commit
  `76d7dc54dbd081e5822b17d9592ce3c0df9a9141`:
  <https://github.com/glwagner/SingleColumnModelCalibration.jl/blob/76d7dc54dbd081e5822b17d9592ce3c0df9a9141/comparison/les_catke_gotm_with_tracers.jl>
- Official LES archive SHA-256:
  `167307a30e161ca98d0f9ded1e248329339c33877174e32ca77d5358ce35e9f8`

The large LES archive is deliberately excluded from the repository.
`prepare_les_profiles.jl` records the interpolation procedure used to
produce the committed CSV files. It reuses the extracted one-metre
archive under `../free_convection_les/input/les_1m/`; the optional
script requires JLD2.jl.
