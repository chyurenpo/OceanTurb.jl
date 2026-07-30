# Mellor--Yamada level 2.5 closure

`OceanTurb.MellorYamada25` implements the GOTM ocean configuration of
the Mellor--Yamada level 2.5 closure. The prognostic turbulence
variables are turbulent kinetic energy \(e=q^2/2\) and the
Mellor--Yamada length variable \(q^2\ell=2e\ell\):

\[
(U,V,T,S,e,q^2\ell).
\]

The turbulence equations are

\[
\partial_t e =
\partial_z\left(S_q q\ell\,\partial_z e\right)
+P+B-\epsilon
\]

and

\[
\partial_t(q^2\ell) =
\partial_z\left(S_\ell q\ell\,
\partial_z(q^2\ell)\right)
+\ell(E_1P+E_3B)
-\frac{q^3}{B_1}
\left[1+E_2\left(\frac{\ell}{L_z}\right)^2\right].
\]

The dissipation rate follows from

\[
\epsilon=(c_\mu^0)^3\frac{e^{3/2}}{\ell}.
\]

Momentum and tracer diffusivities use the local weak-equilibrium
Kantha--Clayson stability functions,

\[
\nu_t=c_\mu(\alpha_M,\alpha_N)\sqrt{e}\ell,\qquad
\kappa_t=c_\mu'(\alpha_M,\alpha_N)\sqrt{e}\ell,
\]

where

\[
\alpha_M=(e/\epsilon)^2M^2,\qquad
\alpha_N=(e/\epsilon)^2N^2.
\]

## Basic use

```julia
using OceanTurb

const MY = OceanTurb.MellorYamada25

grid = UniformGrid(N=128, H=60.0)
parameters = MY.MY25Parameters()

bcs = MY.ModelBoundaryConditions(
    T = BoundaryConditions(
        top = FluxBoundaryCondition(temperature_flux),
    ),
    parameters = parameters,
)

model = MY.Model(
    grid = grid,
    parameters = parameters,
    bcs = bcs,
    constants = Constants(f=1e-4),
    stepper = :BackwardEuler,
)

set!(
    model.solution,
    U = 0.0,
    V = 0.0,
    T = z -> 10.0 + 0.01z,
    S = 35.0,
)

time_step!(model, 1.0, 100)
MY.update_state!(model)
```

Useful diagnostic fields include `model.state.ℓ`, `model.state.ε`,
`model.state.νt`, `model.state.κt`, `model.state.c_mu`,
`model.state.c_mu_p`, `model.state.P`, and `model.state.B`.

## Canonical GOTM configuration

The default model reproduces the closure choices in the pinned GOTM
free-convection cases:

- weak-equilibrium second-moment closure;
- Kantha--Clayson coefficients;
- \(Ri_{st}=0.23\) and a Galperin coefficient of \(0.53\);
- native \(k_{\min}=10^{-10}\) and
  \(\epsilon_{\min}=10^{-12}\) initialization;
- molecular momentum and temperature diffusivities of
  \(1.3\times10^{-6}\) and \(1.4\times10^{-7}\ \mathrm{m^2\,s^{-1}}\);
  and
- logarithmic Neumann conditions with zero TKE flux and GOTM's
  \(q^2\ell\) wall contribution at both boundaries.

Dependent constants are derived with the same arithmetic used by the
native GOTM executable:

\[
c_\mu^0=0.5549055381391544,\qquad
B_1=16.55342282307772,
\]

\[
\kappa=0.4001096562684119,\qquad
E_3=4.825855926750471.
\]

The former MY82 quasi-equilibrium implementation has been removed from
the canonical model.

Matching the closure equations does not make OceanTurb and GOTM
bitwise-identical column solvers. GOTM places its turbulence variables
on interfaces, whereas OceanTurb places them at cell centers. The
finite-volume operators and split update order also differ. In the
five 128-level, 1-second free-convection validations, the diagnosed
mixed-layer base differs by one to two grid cells.

## API

```@docs
OceanTurb.MellorYamada25.Model
OceanTurb.MellorYamada25.ModelBoundaryConditions
OceanTurb.MellorYamada25.Forcing
OceanTurb.MellorYamada25.MY25Parameters
OceanTurb.MellorYamada25.BackgroundDiffusivities
OceanTurb.MellorYamada25.KanthaClaysonCoefficients
OceanTurb.MellorYamada25.kantha_clayson_reference_values
OceanTurb.MellorYamada25.kantha_clayson_quasi_equilibrium
OceanTurb.MellorYamada25.kantha_clayson_weak_equilibrium
OceanTurb.MellorYamada25.kantha_clayson_compute_my25_E3
OceanTurb.MellorYamada25.kantha_clayson_my25_coefficients
OceanTurb.MellorYamada25.LogarithmicNeumannBoundaryConditions
OceanTurb.MellorYamada25.logarithmic_q2l_flux
OceanTurb.MellorYamada25.BottomLogarithmicWallConditions
OceanTurb.MellorYamada25.TopLogarithmicWallConditions
OceanTurb.MellorYamada25.friction_velocity
OceanTurb.MellorYamada25.logarithmic_tke_value
OceanTurb.MellorYamada25.logarithmic_q2l_value
```
