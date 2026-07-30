module CATKE

using OceanTurb

using OceanTurb: minuszero

const nsol = 5
@solution U V T S e

"""
    Model

The Convective Adjustment and Turbulent Kinetic Energy (CATKE)
one-equation closure of Wagner et al. (2025).

Prognostic fields:
- `U`, `V`: horizontal velocity components
- `T`, `S`: temperature and salinity
- `e`: turbulent kinetic energy

Momentum, tracer, TKE-transport, and dissipation length scales are
diagnosed from `e`, stratification, shear, distance from the surface
and bottom, and the filtered destabilizing surface buoyancy flux.
"""
mutable struct Model{P, ML, TE, K0, C, ST, G, TS, S, BC, T, F} <: AbstractModel{TS, G, T}
    clock                    :: Clock{T}
    grid                     :: G
    timestepper              :: TS
    solution                 :: S
    bcs                      :: BC
    parameters               :: P
    mixing_length            :: ML
    tke_equation             :: TE
    background_diffusivities :: K0
    constants                :: C
    state                    :: ST
    forcing                  :: F
end

include("parameters.jl")
include("mixing_lengths.jl")
include("state.jl")
include("diffusivities.jl")
include("tke_equation.jl")
include("boundary_conditions.jl")

addzero(args...) = 0

"""
    Forcing(; U=addzero, V=addzero, T=addzero, S=addzero, e=addzero)

Construct forcing functions for a CATKE `Model`. Each function must
have the signature `forcing(model, i)`, where `i` is a cell index.
"""
Forcing(; U=addzero, V=addzero, T=addzero, S=addzero, e=addzero) =
    (U=U, V=V, T=T, S=S, e=e)

"""
    ModelBoundaryConditions([FT=Float64]; U, V, T, S, e)

Construct boundary conditions for CATKE's five prognostic fields.
The default TKE condition is no flux at the bottom and the calibrated
wind- and convection-driven TKE flux at the surface.
"""
function ModelBoundaryConditions(FT=Float64;
    U = DefaultBoundaryConditions(FT),
    V = DefaultBoundaryConditions(FT),
    T = DefaultBoundaryConditions(FT),
    S = DefaultBoundaryConditions(FT),
    e = CATKETKEBoundaryConditions(FT),
)
    return (U=U, V=V, T=T, S=S, e=e)
end

"""
    Model(; grid, parameters, mixing_length, tke_equation,
            background_diffusivities, constants, stepper, bcs,
            initial_tke, forcing)

Construct an ocean-column CATKE model.

Destructive buoyancy flux, TKE dissipation, and the optional
near-bottom TKE sink are treated as positive implicit sinks.
"""
function Model(;
    grid = UniformGrid(N=64, H=50),
    parameters = CATKEParameters(),
    mixing_length = CATKEMixingLength(),
    tke_equation = CATKEEquation(),
    background_diffusivities = BackgroundDiffusivities(),
    constants = Constants(),
    stepper = :BackwardEuler,
    bcs = ModelBoundaryConditions(eltype(grid)),
    initial_tke = parameters.minimum_tke,
    forcing = Forcing(),
)
    @assert grid.N ≥ 2 "CATKE.Model requires at least two vertical cells."
    validate_parameters(parameters)
    validate_mixing_length(mixing_length)
    validate_tke_equation(tke_equation)
    validate_background_diffusivities(background_diffusivities)

    solution = Solution((CellField(grid) for _ = 1:nsol)...)
    initialize_tke!(solution, grid, parameters, initial_tke)

    Kphi = (
        U = KU,
        V = KV,
        T = KT,
        S = KS,
        e = Ke,
    )

    Rphi = (
        U = RU,
        V = RV,
        T = RT,
        S = RS,
        e = Re,
    )

    Lphi = (
        U = minuszero,
        V = minuszero,
        T = minuszero,
        S = minuszero,
        e = Le,
    )

    equation = Equation(K=Kphi, R=Rphi, L=Lphi, update=update_state!)
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
        mixing_length,
        tke_equation,
        background_diffusivities,
        constants,
        state,
        forcing,
    )

    OceanTurb.update!(model)
    return model
end

export Model,
       ModelBoundaryConditions,
       Forcing,
       CATKEParameters,
       CATKEMixingLength,
       CATKEEquation,
       BackgroundDiffusivities,
       CATKETKEBoundaryConditions,
       friction_velocity,
       surface_buoyancy_flux,
       surface_tke_flux,
       stability_scale,
       update_state!,
       shear_production,
       buoyancy_production,
       dissipation,
       momentum_mixing_length,
       tracer_mixing_length,
       tke_mixing_length,
       dissipation_length

end # module CATKE
