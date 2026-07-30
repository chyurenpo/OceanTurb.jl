# Convective Adjustment and Turbulent Kinetic Energy

`OceanTurb.CATKE` implements the one-equation Convective Adjustment
and Turbulent Kinetic Energy closure described by Wagner et al. (2025).
The prognostic fields are

\[
(U, V, T, S, e),
\]

where \(e\) is turbulent kinetic energy. The TKE equation is

\[
\partial_t e =
\partial_z(K_e\partial_z e)
+K_u M^2-K_cN^2-\frac{e^{3/2}}{\ell_D}.
\]

CATKE diagnoses distinct momentum, tracer, TKE-transport, and
dissipation length scales. The corresponding diffusivities are

\[
K_u=\ell_u\sqrt{e},\qquad
K_c=\ell_c\sqrt{e},\qquad
K_e=\ell_e\sqrt{e}.
\]

Each mixing length chooses the larger of a shear/stability component
and a dynamic convective component:

\[
\ell_\psi=\max(\ell_\psi^{\rm shear},
               \ell_\psi^{\rm conv}).
\]

The shear component is limited by stratification and distance from
the surface and bottom, then multiplied by a calibrated, piecewise
Richardson-number stability function. The convective component uses
the filtered destabilizing surface buoyancy flux,

\[
\ell_\psi^{\rm conv}
\sim C_\psi^{\rm conv}
\frac{e^{3/2}}{\widetilde J_b+J_b^{\min}},
\]

with separate penetrative tracer mixing immediately beneath a
convecting layer and a shear-suppression factor.

## Basic use

```julia
using OceanTurb

const CA = OceanTurb.CATKE

grid = UniformGrid(N=128, H=100.0)
constants = Constants(f=1e-4)
Jb = 1e-7

bcs = CA.ModelBoundaryConditions(
    T = BoundaryConditions(
        top = FluxBoundaryCondition(
            Jb / (constants.g * constants.α),
        ),
    ),
)

model = CA.Model(
    grid = grid,
    constants = constants,
    bcs = bcs,
    stepper = :BackwardEuler,
    initial_tke = 1e-6,
)

set!(
    model.solution,
    U = 0.0,
    V = 0.0,
    T = z -> 10.0 + 0.01*z,
    S = 35.0,
)

time_step!(model, 1.0, 100)
CA.update_state!(model)
```

The prognostic turbulence field is `model.solution.e`. Diagnostics
include `model.state.ell_u`, `ell_c`, `ell_e`, `ell_D`, `nu_t`,
`kappa_t`, `K_e`, `P`, `B`, and `epsilon`.

## Surface forcing

CATKE's default surface TKE flux is

\[
J_e=-C_{Wu}u_\tau^3
    -C_{Wb}\Delta z\max(J_b,0).
\]

It is diagnosed from the surface boundary conditions for momentum,
temperature, and salinity. Positive \(J_b\) denotes buoyancy loss and
drives convection. The destabilizing surface buoyancy flux is filtered
over a diagnosed convective mixing time before it enters CATKE's
convective length scales.

## Parameters and current scope

`CATKEMixingLength` and `CATKEEquation` contain the calibrated
coefficients. Their defaults follow the maintained Oceananigans
implementation, including bottom-distance limitation and the optional
near-bottom TKE sink. Set
`CATKEEquation(C_bottom_dissipation=0)` to recover a no-flux bottom TKE
budget closer to the surface-focused formulation in the 2025 paper.

This OceanTurb implementation assumes a linear temperature-salinity
equation of state and a one-dimensional column. Stokes-shear
production and explicit nonlocal mass-flux transport are outside its
current scope.

## Reference

G. L. Wagner et al. (2025), “Formulation and calibration of CATKE, a
one-equation parameterization for microscale ocean mixing,” *Journal
of Advances in Modeling Earth Systems*, 17, e2024MS004522,
[doi:10.1029/2024MS004522](https://doi.org/10.1029/2024MS004522).

## API

```@docs
OceanTurb.CATKE.Model
OceanTurb.CATKE.ModelBoundaryConditions
OceanTurb.CATKE.Forcing
OceanTurb.CATKE.CATKEParameters
OceanTurb.CATKE.CATKEMixingLength
OceanTurb.CATKE.CATKEEquation
OceanTurb.CATKE.BackgroundDiffusivities
OceanTurb.CATKE.CATKETKEBoundaryConditions
OceanTurb.CATKE.friction_velocity
OceanTurb.CATKE.surface_buoyancy_flux
OceanTurb.CATKE.surface_tke_flux
OceanTurb.CATKE.stability_scale
```
