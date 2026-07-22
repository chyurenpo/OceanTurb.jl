module MellorYamada25

using OceanTurb

using OceanTurb: minuszero
import OceanTurb: oncell, onface

const nsol = 6
@solution U V T S e q2l

"""
    Model

A first-pass OceanTurb implementation of the classical Mellor–Yamada
level-2.5 closure, following GOTM's `q2over2eq`, `lengthscaleeq`, and
quasi-equilibrium MY82 stability-function organization.

Prognostic fields:
- `U`, `V`: horizontal velocity components
- `T`, `S`: temperature and salinity
- `e = q²/2`: turbulent kinetic energy
- `q2l = q² ℓ = 2eℓ`: Mellor–Yamada length-scale variable
"""
mutable struct Model{P, SC, K0, C, ST, G, TS, S, BC, T} <: AbstractModel{TS, G, T}
    clock                    :: Clock{T}
    grid                     :: G
    timestepper              :: TS
    solution                 :: S
    bcs                      :: BC
    parameters               :: P
    second_moment_coefficients :: SC
    background_diffusivities :: K0
    constants                :: C
    state                    :: ST
end

include("parameters.jl")
include("stability_functions.jl")
include("state.jl")
include("diffusivities.jl")
include("turbulence_equations.jl")
include("wall_models.jl")

"""
    ModelBoundaryConditions([FT=Float64]; U, V, T, S, e, q2l)

Construct boundary conditions for the six prognostic MY2.5 fields.
The turbulence variables default to zero-flux conditions. For a
wall-bounded calculation, replace the lower `e` and `q2l` conditions
with the logarithmic-wall helpers in `wall_models.jl`.
"""
function ModelBoundaryConditions(FT=Float64;
    U    = DefaultBoundaryConditions(FT),
    V    = DefaultBoundaryConditions(FT),
    T    = DefaultBoundaryConditions(FT),
    S    = DefaultBoundaryConditions(FT),
    e    = ZeroFluxBoundaryConditions(FT),
    q2l  = ZeroFluxBoundaryConditions(FT),
)
    return (U=U, V=V, T=T, S=S, e=e, q2l=q2l)
end

"""
    Model(; grid, parameters, background_diffusivities, constants,
            stepper=:BackwardEuler, bcs, initial_tke)

Construct a Mellor–Yamada 2.5 model.

This implementation uses OceanTurb's backward-Euler diffusion and
linear-sink machinery. Destructive TKE and q²ℓ terms are represented
as positive implicit sink coefficients.
"""
function Model(;
    grid = UniformGrid(N=64, H=50),
    parameters = MY25Parameters(),
    second_moment_coefficients = MY82Coefficients(),
    background_diffusivities = BackgroundDiffusivities(),
    constants = Constants(),
    stepper = :BackwardEuler,
    bcs = ModelBoundaryConditions(eltype(grid)),
    initial_tke = max(parameters.k_min, 1e-8),
)
    @assert grid.N ≥ 2 "MellorYamada25.Model requires at least two vertical cells."

    solution = Solution((CellField(grid) for _ = 1:nsol)...)
    initialize_turbulence!(solution, grid, parameters, initial_tke)

    Kϕ = (
        U    = KU,
        V    = KV,
        T    = KT,
        S    = KS,
        e    = Ke,
        q2l  = Kq2l,
    )

    Rϕ = (
        U    = RU,
        V    = RV,
        T    = RT,
        S    = RS,
        e    = Re,
        q2l  = Rq2l,
    )

    Lϕ = (
        U    = minuszero,
        V    = minuszero,
        T    = minuszero,
        S    = minuszero,
        e    = Le,
        q2l  = Lq2l,
    )

    equation = Equation(K=Kϕ, R=Rϕ, L=Lϕ, update=update_state!)
    lhs = OceanTurb.build_lhs(solution)
    timestepper = Timestepper(stepper, equation, solution, lhs)

    Tclock = eltype(grid)
    clock = Clock(zero(Tclock), 0)
    state = State(grid)

    model = Model(
        clock,
        grid,
        timestepper,
        solution,
        bcs,
        parameters,
        second_moment_coefficients,
        background_diffusivities,
        constants,
        state,
    )

    update_state!(model)
    return model
end

export Model,
       ModelBoundaryConditions,
       MY25Parameters,
       BackgroundDiffusivities,
       MY82Coefficients,
       BottomLogarithmicWallConditions,
       TopLogarithmicWallConditions,
       friction_velocity,
       logarithmic_tke_value,
       logarithmic_q2l_value,
       update_state!,
       shear_production,
       buoyancy_production,
       dissipation

end # module MellorYamada25
