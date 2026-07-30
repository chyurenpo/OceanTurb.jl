# Free-convection LES example

This experiment compares OceanTurb CATKE, KPP, Mellor--Yamada 2.5,
k--omega, and k--epsilon with free-convection LES profiles. These cases
underlie Figure 5 of Wagner et al. (2025). Every closure starts from
the same official LES buoyancy and passive-tracer profiles at 600
seconds and uses the same 128-cell grid and 10-second time step.

This is a physical-case comparison, not an exact reproduction of every
manuscript curve. KPP and the additional two-equation closures are the
OceanTurb implementations.

## Configuration

| Setting | Value |
|---|---:|
| Domain depth | 256 m |
| OceanTurb cells | 128 |
| Start time | 600 s |
| Default time step | 10 s |
| Final times | 6, 12, 24, 48, 72 h |
| Surface buoyancy fluxes | \(9.6, 4.8, 2.4, 1.2, 0.87\times10^{-7}\ {\rm m^2\,s^{-3}}\) |
| Coriolis parameter | \(10^{-4}\ {\rm s^{-1}}\) |
| Deep buoyancy gradient | \(2\times10^{-6}\ {\rm s^{-2}}\) |
| Initial TKE for prognostic closures | \(10^{-9}\ {\rm m^2\,s^{-2}}\) |

The LES archive has no corresponding \(k\), \(q^2\ell\), \(\omega\),
or \(\epsilon\) profiles. CATKE therefore starts from the paper's TKE
floor. The two-equation closures use that TKE together with their
canonical native second-variable floors.

## Run

The processed LES inputs required by the runner are committed in
`data/`. From the OceanTurb repository root, run:

```bash
julia --project=. examples/free_convection_les/run.jl
```

For a shorter single-suite check:

```bash
OCEANTURB_SUITES=6 julia --project=. \
  examples/free_convection_les/run.jl
```

Set `OCEANTURB_FREE_CONVECTION_DT` to override the default 10-second
time step.

The runner writes profiles and a summary to `output/dt10s/`. Generated
output is ignored by Git.

To regenerate the comparison figures after running the cases:

```bash
cd examples/free_convection_les
jupyter execute --inplace --timeout=600 plot_profiles.ipynb
```

The committed representative result is
`figures/free_convection_closures.png`.

## Data provenance

- CATKE manuscript: <https://glwagner.github.io/assets/pdf/CATKE.pdf>
- Archived comparison code, commit
  `76d7dc54dbd081e5822b17d9592ce3c0df9a9141`:
  <https://github.com/glwagner/SingleColumnModelCalibration.jl/blob/76d7dc54dbd081e5822b17d9592ce3c0df9a9141/comparison/les_catke_gotm_with_tracers.jl>
- Official LES archive SHA-256:
  `167307a30e161ca98d0f9ded1e248329339c33877174e32ca77d5358ce35e9f8`

The 94 MB LES archive is deliberately excluded from the repository.
`prepare_les_profiles.jl` records the interpolation procedure used to
produce the committed CSV files. To repeat that preprocessing, provide
the extracted one-metre LES suites under `input/les_1m/`; the optional
script requires JLD2.jl.
