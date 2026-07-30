# Oceanic k--epsilon closure

`OceanTurb.KEpsilon` implements GOTM's Rodi \(k\)-\(\epsilon\)
closure with the weak-equilibrium Canuto-A second-moment model.

The six prognostic fields are

\[
(U, V, T, S, k, \epsilon),
\]

where \(k\) is turbulent kinetic energy and \(\epsilon\) is its
dissipation rate. The turbulence equations are

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
\partial_t \epsilon =
\partial_z \left[
  \left(K_{\epsilon0} + \frac{\nu_t}{\sigma_\epsilon}\right)
  \partial_z \epsilon
\right]
+ \frac{\epsilon}{k}
  \left(c_{\epsilon1} P + c_{\epsilon3} B
  - c_{\epsilon2}\epsilon\right).
\]

The dissipative length is

\[
\ell = (c_\mu^0)^3\frac{k^{3/2}}{\epsilon}.
\]

GOTM's Canuto-A stability functions depend on
\(\alpha_M=(k/\epsilon)^2M^2\) and
\(\alpha_N=(k/\epsilon)^2N^2\):

\[
\nu_t = c_\mu(\alpha_M,\alpha_N)\sqrt{k}\ell,\qquad
\kappa_t =
c_\mu'(\alpha_M,\alpha_N)\sqrt{k}\ell.
\]

Stable buoyancy destruction and dissipation are treated as positive
implicit sinks. The optional Galperin limit increases \(\epsilon\)
when necessary to enforce

\[
\ell \le c_G \sqrt{\frac{2k}{N^2}}
\]

under stable stratification.

## Basic use

```julia
using OceanTurb

const KE = OceanTurb.KEpsilon

grid = UniformGrid(N=128, H=60.0)
parameters = KE.KEpsilonParameters(
    z0_bottom = 1e-3,
    length_limit = true,
)

model = KE.Model(
    grid = grid,
    parameters = parameters,
    constants = Constants(f=1e-4),
    stepper = :BackwardEuler,
)

set!(
    model.solution,
    U = z -> 0.1 * (z + grid.H) / grid.H,
    V = 0.0,
    T = z -> 10.0 + 0.01*z,
    S = 35.0,
)

time_step!(model, 1.0, 100)
KE.update_state!(model)
```

The prognostic turbulence fields are `model.solution.k` and
`model.solution.epsilon`. Useful diagnostic fields include
`model.state.ell`, `model.state.nu_t`, `model.state.kappa_t`,
`model.state.c_mu`, `model.state.c_mu_prime`, `model.state.P`, and
`model.state.B`.

## GOTM formulation

`KEpsilonParameters()` and `KE.Model()` use the GOTM/Canuto-A
formulation by default. The canonical configuration reproduces
GOTM's:

- weak-equilibrium Canuto-A momentum and tracer stability functions;
- derived \(c_\mu^0\), shear-free \(c_\mu\), von Kármán constant, and
  stable \(c_{\epsilon3}\);
- \(k_{\min}=10^{-10}\) and
  \(\epsilon_{\min}=10^{-12}\) native initialization;
- molecular momentum and temperature diffusivities used by the
  reference free-convection cases; and
- logarithmic Neumann conditions: zero \(k\) flux and inward
  \(\epsilon\) flux
  \[
  \frac{(c_\mu^0)^4 k^2}
       {\sigma_\epsilon(\Delta z/2+z_0)}
  \]
  at both boundaries.

For the pinned GOTM configuration, the derived coefficients are

\[
c_\mu^0 = 0.5264647108353244,\qquad
\kappa = 0.4158737996737366,\qquad
c_{\epsilon3}^{stable} = -0.6209121262490007.
\]

Matching the closure equations does not make OceanTurb and GOTM
bitwise-identical column solvers. GOTM places prognostic turbulence
variables on interfaces, whereas OceanTurb places them at cell
centers; the update order and finite-volume operators also differ.

## Logarithmic wall conditions

```julia
u_tau = 0.005
k_bottom, epsilon_bottom =
    KE.BottomLogarithmicWallConditions(
        u_tau;
        parameters = parameters,
    )

turbulence_bcs = KE.ModelBoundaryConditions(
    k = BoundaryConditions(
        bottom = k_bottom,
        top = FluxBoundaryCondition(0.0),
    ),
    epsilon = BoundaryConditions(
        bottom = epsilon_bottom,
        top = FluxBoundaryCondition(0.0),
    ),
)

model = KE.Model(
    grid = grid,
    parameters = parameters,
    bcs = turbulence_bcs,
)
```

`u_tau` may be a number or a function `u_tau(model)`. Because
boundary-condition types are part of OceanTurb's model type, value
wall conditions must be supplied when the model is constructed.

The basic-use example above shows a complete runnable column. The
free-convection and rotating strong-wind workflows under `examples/`
exercise this closure alongside the other OceanTurb models.

## Current scope

The former constant-stability implementation has been removed. The
canonical implementation does not yet include surface-wave TKE
injection, Stokes-shear production, or GOTM's full
interface-staggered numerical layout.

## API

```@docs
OceanTurb.KEpsilon.Model
OceanTurb.KEpsilon.ModelBoundaryConditions
OceanTurb.KEpsilon.Forcing
OceanTurb.KEpsilon.KEpsilonParameters
OceanTurb.KEpsilon.BackgroundDiffusivities
OceanTurb.KEpsilon.LogarithmicNeumannBoundaryConditions
OceanTurb.KEpsilon.logarithmic_epsilon_flux
OceanTurb.KEpsilon.BottomLogarithmicWallConditions
OceanTurb.KEpsilon.TopLogarithmicWallConditions
OceanTurb.KEpsilon.friction_velocity
OceanTurb.KEpsilon.logarithmic_tke_value
OceanTurb.KEpsilon.logarithmic_epsilon_value
OceanTurb.KEpsilon.canuto_a_coefficients
OceanTurb.KEpsilon.canuto_a_reference_values
OceanTurb.KEpsilon.canuto_a_weak_equilibrium
OceanTurb.KEpsilon.canuto_a_quasi_equilibrium
OceanTurb.KEpsilon.canuto_a_compute_c3_stable
OceanTurb.KEpsilon.canuto_a_kepsilon_coefficients
```
