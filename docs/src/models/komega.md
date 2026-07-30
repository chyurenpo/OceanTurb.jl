# Oceanic k--omega closure

`OceanTurb.KOmega` implements the Wilcox \(k\)-\(\omega\) closure as
extended to stratified geophysical flows by Umlauf, Burchard, and
Hutter (2003). This is an ocean-column closure, not the engineering
\(k\)-\(\omega\) SST model.

The six prognostic fields are horizontal velocity, temperature,
salinity, turbulent kinetic energy, and inverse turbulent time scale:

\[
(U, V, T, S, k, \omega).
\]

The turbulence equations are

\[
\partial_t k =
\partial_z \left[
  \left(K_{k0} + \frac{\nu_t}{\sigma_k}\right)
  \partial_z k
\right]
+ P + B - \epsilon
\]

and

\[
\partial_t \omega =
\partial_z \left[
  \left(K_{\omega0} + \frac{\nu_t}{\sigma_\omega}\right)
  \partial_z \omega
\right]
+ \frac{\omega}{k}
  \left(c_{\omega1} P + c_{\omega3} B
  - c_{\omega2}\epsilon\right).
\]

With the Canuto-A weak-equilibrium stability functions,

\[
\epsilon = (c_\mu^0)^4 k \omega, \qquad
\ell = \frac{\sqrt{k}}{c_\mu^0\omega}, \qquad
\nu_t = c_\mu \sqrt{k}\ell, \qquad
\kappa_t = c_\mu' \sqrt{k}\ell.
\]

Stable buoyancy destruction and dissipation are treated as positive
implicit sinks. The optional Galperin limit prevents the turbulent
length scale from exceeding

\[
\ell_{\max} = c_G \sqrt{\frac{2k}{N^2}}
\]

under stable stratification.

## Basic use

```julia
using OceanTurb

const KW = OceanTurb.KOmega

grid = UniformGrid(N=128, H=60.0)
parameters = KW.KOmegaParameters(
    z0_bottom = 1e-3,
    length_limit = true,
)

model = KW.Model(
    grid = grid,
    parameters = parameters,
    constants = Constants(f=1e-4),
    stepper = :BackwardEuler,
    initial_tke = 1e-6,
    initial_omega = 0.05,
)

set!(
    model.solution,
    U = z -> 0.1 * (z + grid.H) / grid.H,
    V = 0.0,
    T = z -> 10.0 + 0.01*z,
    S = 35.0,
)

time_step!(model, 1.0, 100)
KW.update_state!(model)
```

The prognostic turbulence fields are `model.solution.k` and
`model.solution.omega`. Useful diagnostic fields include
`model.state.epsilon`, `model.state.ell`, `model.state.nu_t`,
`model.state.kappa_t`, `model.state.c_mu`,
`model.state.c_mu_prime`, `model.state.P`, and `model.state.B`.

## GOTM formulation

`KOmegaParameters()` and `KOmega.Model()` use the GOTM/Canuto-A
formulation by default:

```julia
parameters = KW.KOmegaParameters()

bcs = KW.ModelBoundaryConditions(
    Float64;
    T = BoundaryConditions(
        bottom = FluxBoundaryCondition(0.0),
        top = FluxBoundaryCondition(surface_temperature_flux),
    ),
    parameters = parameters,
)

model = KW.Model(
    grid = grid,
    parameters = parameters,
    bcs = bcs,
    stepper = :BackwardEuler,
)
```

The canonical configuration reproduces GOTM's:

- weak-equilibrium Canuto-A momentum and tracer stability functions;
- derived \(c_\mu^0\), shear-free \(c_\mu\), von Kármán constant, and
  stable \(c_{\omega3}\);
- \(k_{\min}=10^{-10}\) and
  \(\epsilon_{\min}=10^{-12}\) native initialization;
- molecular momentum and temperature diffusivities used by the
  reference free-convection cases; and
- logarithmic Neumann conditions: zero \(k\) flux and inward
  \(\omega\) flux
  \(k/[\sigma_\omega(\Delta z/2+z_0)]\) at both boundaries.

For the pinned GOTM configuration, the derived coefficients are

\[
c_\mu^0 = 0.5264647108353244,\qquad
\kappa = 0.3925605571319366,\qquad
c_{\omega3}^{stable} = -0.6386116064525458.
\]

Matching the closure equations does not make OceanTurb and GOTM
bitwise-identical column solvers. GOTM places prognostic turbulence
variables on interfaces, whereas OceanTurb places them at cell
centers; the update order and finite-volume operators also differ.

## Logarithmic wall conditions

```julia
u_tau = 0.005
k_bottom, omega_bottom =
    KW.BottomLogarithmicWallConditions(
        u_tau;
        parameters = parameters,
    )

turbulence_bcs = KW.ModelBoundaryConditions(
    k = BoundaryConditions(
        bottom = k_bottom,
        top = FluxBoundaryCondition(0.0),
    ),
    omega = BoundaryConditions(
        bottom = omega_bottom,
        top = FluxBoundaryCondition(0.0),
    ),
)

model = KW.Model(
    grid = grid,
    parameters = parameters,
    bcs = turbulence_bcs,
)
```

`u_tau` may be a number or a function `u_tau(model)`, which permits
coupling the turbulence wall conditions to a dynamically diagnosed
bottom stress. Because boundary-condition types are part of
OceanTurb's model type, value wall conditions must be supplied when
the model is constructed.

The basic-use example above shows a complete runnable column. The
free-convection and rotating strong-wind workflows under `examples/`
exercise this closure alongside the other OceanTurb models.

## Current scope

The former constant-stability implementation has been removed. The
canonical implementation does not yet include surface-wave TKE
injection, Stokes-shear production, or GOTM's full
interface-staggered numerical layout.

## Reference

L. Umlauf, H. Burchard, and K. Hutter (2003), “Extending the
\(k\)-\(\omega\) turbulence model towards oceanic applications,”
*Ocean Modelling*, 5, 195–218,
[doi:10.1016/S1463-5003(02)00039-2](https://doi.org/10.1016/S1463-5003(02)00039-2).

## API

```@docs
OceanTurb.KOmega.Model
OceanTurb.KOmega.ModelBoundaryConditions
OceanTurb.KOmega.Forcing
OceanTurb.KOmega.KOmegaParameters
OceanTurb.KOmega.BackgroundDiffusivities
OceanTurb.KOmega.LogarithmicNeumannBoundaryConditions
OceanTurb.KOmega.logarithmic_omega_flux
OceanTurb.KOmega.BottomLogarithmicWallConditions
OceanTurb.KOmega.TopLogarithmicWallConditions
OceanTurb.KOmega.friction_velocity
OceanTurb.KOmega.logarithmic_tke_value
OceanTurb.KOmega.logarithmic_omega_value
OceanTurb.KOmega.canuto_a_coefficients
OceanTurb.KOmega.canuto_a_reference_values
OceanTurb.KOmega.canuto_a_weak_equilibrium
OceanTurb.KOmega.canuto_a_quasi_equilibrium
OceanTurb.KOmega.canuto_a_compute_c3_stable
OceanTurb.KOmega.canuto_a_komega_coefficients
```
