module KOmega

using OceanTurb

using OceanTurb: minuszero
import OceanTurb: onface

const nsol = 6
@solution U V T S k omega

"""
    Model

An oceanic `k`-`omega` two-equation turbulence closure following the
Wilcox (1988) model as extended to stratified geophysical flows by
Umlauf, Burchard, and Hutter (2003).

Prognostic fields:
- `U`, `V`: horizontal velocity components
- `T`, `S`: temperature and salinity
- `k`: turbulent kinetic energy
- `omega`: inverse turbulent time scale
"""
mutable struct Model{P, K0, C, ST, G, TS, S, BC, T, F} <: AbstractModel{TS, G, T}
    clock                    :: Clock{T}
    grid                     :: G
    timestepper              :: TS
    solution                 :: S
    bcs                      :: BC
    parameters               :: P
    background_diffusivities :: K0
    constants                :: C
    state                    :: ST
    forcing                  :: F
end

include("../canuto_a_stability_functions.jl")
include("parameters.jl")
include("state.jl")
include("diffusivities.jl")
include("turbulence_equations.jl")
include("wall_models.jl")

addzero(args...) = 0

"""
    Forcing(; U=addzero, V=addzero, T=addzero, S=addzero)

Construct mean-field forcing functions for a k-omega `Model`. Each
function must have the signature `forcing(model, i)`, where `i` is a
cell index. The default forcing is identically zero.
"""
Forcing(; U=addzero, V=addzero, T=addzero, S=addzero) =
    (U=U, V=V, T=T, S=S)

"""
    ModelBoundaryConditions([FT=Float64];
        U, V, T, S, k, omega, parameters)

Construct boundary conditions for the six prognostic k-omega fields.
The turbulence variables default to GOTM's logarithmic Neumann
conditions.
"""
function ModelBoundaryConditions(FT=Float64;
    U     = DefaultBoundaryConditions(FT),
    V     = DefaultBoundaryConditions(FT),
    T     = DefaultBoundaryConditions(FT),
    S     = DefaultBoundaryConditions(FT),
    k     = nothing,
    omega = nothing,
    parameters = KOmegaParameters(),
)
    default_k, default_omega =
        LogarithmicNeumannBoundaryConditions(
            FT;
            parameters=parameters,
        )
    k === nothing && (k = default_k)
    omega === nothing && (omega = default_omega)
    return (U=U, V=V, T=T, S=S, k=k, omega=omega)
end

"""
    Model(; grid, parameters, background_diffusivities, constants,
            stepper=:BackwardEuler, bcs, initial_tke, initial_omega,
            initial_epsilon, forcing)

Construct an oceanic k-omega model.

Diffusion and destructive turbulence terms are treated with OceanTurb's
backward-Euler diffusion and positive implicit-sink machinery.
"""
function Model(;
    grid = UniformGrid(N=64, H=50),
    parameters = KOmegaParameters(),
    background_diffusivities = BackgroundDiffusivities(),
    constants = Constants(),
    stepper = :BackwardEuler,
    bcs = ModelBoundaryConditions(
        eltype(grid);
        parameters=parameters,
    ),
    initial_tke = parameters.k_min,
    initial_omega = nothing,
    initial_epsilon =
        initial_omega === nothing ?
        parameters.epsilon_min :
        nothing,
    forcing = Forcing(),
)
    @assert grid.N ≥ 2 "KOmega.Model requires at least two vertical cells."
    validate_parameters(parameters)
    validate_background_diffusivities(background_diffusivities)
    @assert initial_omega === nothing || initial_epsilon === nothing "Specify only one of initial_omega and initial_epsilon."

    if initial_epsilon !== nothing
        initial_omega =
            max(initial_epsilon, parameters.epsilon_min) /
            (
                parameters.c_mu0^4 *
                max(initial_tke, parameters.k_min)
            )
    end

    solution = Solution((CellField(grid) for _ = 1:nsol)...)
    initialize_turbulence!(
        solution,
        grid,
        parameters,
        initial_tke,
        initial_omega,
    )

    Kphi = (
        U     = KU,
        V     = KV,
        T     = KT,
        S     = KS,
        k     = Kk,
        omega = Komega,
    )

    Rphi = (
        U     = RU,
        V     = RV,
        T     = RT,
        S     = RS,
        k     = Rk,
        omega = Romega,
    )

    Lphi = (
        U     = minuszero,
        V     = minuszero,
        T     = minuszero,
        S     = minuszero,
        k     = Lk,
        omega = Lomega,
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
        background_diffusivities,
        constants,
        state,
        forcing,
    )

    # Populate diagnostics and solution ghost cells consistently with the
    # configured boundary conditions before the first time step.
    OceanTurb.update!(model)
    return model
end

export Model,
       ModelBoundaryConditions,
       Forcing,
       KOmegaParameters,
       BackgroundDiffusivities,
       canuto_a_coefficients,
       canuto_a_reference_values,
       canuto_a_weak_equilibrium,
       canuto_a_quasi_equilibrium,
       canuto_a_compute_c3_stable,
       canuto_a_komega_coefficients,
       LogarithmicNeumannBoundaryConditions,
       logarithmic_omega_flux,
       BottomLogarithmicWallConditions,
       TopLogarithmicWallConditions,
       friction_velocity,
       logarithmic_tke_value,
       logarithmic_omega_value,
       update_state!,
       shear_production,
       buoyancy_production,
       dissipation,
       turbulent_length_scale

end # module KOmega
