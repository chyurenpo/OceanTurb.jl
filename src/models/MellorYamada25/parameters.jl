"""
Kantha and Clayson (1994) pressure-strain and scalar-flux coefficients
in the notation used by GOTM.

GOTM declares these constants as default-real Fortran literals before
assigning them to double-precision variables. The constructor preserves
that promotion so derived constants agree with the native executable.
"""
struct KanthaClaysonCoefficients{T} <: AbstractParameters
    cc1 :: T
    cc2 :: T
    cc3 :: T
    cc4 :: T
    cc5 :: T
    cc6 :: T

    ct1 :: T
    ct2 :: T
    ct3 :: T
    ct4 :: T
    ct5 :: T
    ctt :: T
end

function KanthaClaysonCoefficients(
    ::Type{T}=Float64,
) where T
    reference_real(x) = T(Float32(x))

    return KanthaClaysonCoefficients{T}(
        reference_real(6.0000),
        reference_real(0.3200),
        reference_real(0.0000),
        reference_real(0.0000),
        reference_real(0.0000),
        reference_real(0.0000),
        reference_real(3.7280),
        reference_real(0.7000),
        reference_real(0.7000),
        reference_real(0.0000),
        reference_real(0.2000),
        reference_real(0.6102),
    )
end

"""
Parameters for OceanTurb's GOTM-compatible Mellor-Yamada level-2.5
closure.

The canonical defaults reproduce the GOTM configuration used by the
free-convection reference cases:

- weak-equilibrium second-moment closure,
- Kantha-Clayson coefficients,
- `Ri_st = 0.23`,
- a Galperin coefficient of `0.53`, and
- dynamically derived `B1`, `E3`, and von Karman constant.
"""
struct MY25Parameters{T} <: AbstractParameters
    Ri_st       :: T
    c_mu0       :: T
    c_mu_shear_free :: T

    B1          :: T
    Sq          :: T
    Sl          :: T
    E1          :: T
    E2          :: T
    E3          :: T

    κ           :: T
    galperin    :: T
    length_limit :: Bool

    k_min       :: T
    eps_min     :: T
    q2l_min     :: T

    z0_bottom   :: T
    z0_top      :: T

    # :parabolic, :triangular, or :infinite_depth
    wall_length :: Symbol
end

"""
    MY25Parameters(; kwargs...)

Construct the canonical GOTM/Kantha-Clayson MY2.5 parameters. `B1`,
`E3`, and `kappa` are derived exactly as in GOTM unless explicitly
overridden. The Unicode keyword `κ` is retained as an alias for
`kappa`.
"""
function MY25Parameters(;
    Ri_st = 0.23,
    Sq = 0.2,
    Sl = 0.2,
    E1 = 1.8,
    E2 = 1.33,
    B1 = nothing,
    E3 = nothing,
    kappa = nothing,
    κ = nothing,
    galperin = 0.53,
    length_limit = true,
    k_min = 1e-10,
    eps_min = 1e-12,
    q2l_min = nothing,
    z0_bottom = 1.5e-3,
    z0_top = 2.0e-2,
    wall_length = :parabolic,
)
    kappa !== nothing && κ !== nothing &&
        throw(ArgumentError("Specify only one of kappa and κ."))
    κ !== nothing && (kappa = κ)

    coefficients = kantha_clayson_my25_coefficients(
        Ri_st=Ri_st,
        Sq=Sq,
        Sl=Sl,
        E1=E1,
        E2=E2,
    )

    B1 === nothing && (B1 = coefficients.B1)
    E3 === nothing && (E3 = coefficients.E3)
    kappa === nothing &&
        (kappa = sqrt((E2 - E1 + 1) / (Sl * B1)))

    c_mu0 = coefficients.c_mu0
    c_mu_shear_free = coefficients.c_mu_shear_free
    cde_value = c_mu0^3
    length_min = cde_value * k_min^(3/2) / eps_min
    q2l_min === nothing &&
        (q2l_min = 2 * k_min * length_min)

    T = promote_type(
        typeof(Ri_st),
        typeof(c_mu0),
        typeof(c_mu_shear_free),
        typeof(B1),
        typeof(Sq),
        typeof(Sl),
        typeof(E1),
        typeof(E2),
        typeof(E3),
        typeof(kappa),
        typeof(galperin),
        typeof(k_min),
        typeof(eps_min),
        typeof(q2l_min),
        typeof(z0_bottom),
        typeof(z0_top),
    )

    parameters = MY25Parameters{T}(
        T(Ri_st),
        T(c_mu0),
        T(c_mu_shear_free),
        T(B1),
        T(Sq),
        T(Sl),
        T(E1),
        T(E2),
        T(E3),
        T(kappa),
        T(galperin),
        length_limit,
        T(k_min),
        T(eps_min),
        T(q2l_min),
        T(z0_bottom),
        T(z0_top),
        wall_length,
    )

    validate_parameters(parameters)
    return parameters
end

"""
Molecular/background diffusivities in square metres per second.

The mean-flow values match GOTM's defaults. GOTM adds no separate
background diffusivity to the `q²/2` and `q²l` equations.
"""
Base.@kwdef struct BackgroundDiffusivities{T} <: AbstractParameters
    momentum :: T = 1.3e-6
    tracer   :: T = 1.4e-7
    tke      :: T = 0.0
    q2l      :: T = 0.0
end

@inline cde(p) = p.c_mu0^3
@inline cm0(p) = p.c_mu0
@inline minimum_length(p) =
    cde(p) * p.k_min^(3/2) / p.eps_min

function validate_parameters(p)
    @assert p.Ri_st > 0 "Ri_st must be positive."
    @assert p.c_mu0 > 0 "c_mu0 must be positive."
    @assert p.c_mu_shear_free > 0 "c_mu_shear_free must be positive."
    @assert p.B1 > 0 "B1 must be positive."
    @assert p.Sq > 0 "Sq must be positive."
    @assert p.Sl > 0 "Sl must be positive."
    @assert p.E1 ≥ 0 "E1 must be non-negative."
    @assert p.E2 ≥ 0 "E2 must be non-negative."
    @assert isfinite(p.E3) "E3 must be finite."
    @assert p.κ > 0 "kappa must be positive."
    @assert p.galperin > 0 "galperin must be positive."
    @assert p.k_min > 0 "k_min must be positive."
    @assert p.eps_min > 0 "eps_min must be positive."
    @assert p.q2l_min > 0 "q2l_min must be positive."
    @assert p.z0_bottom ≥ 0 "z0_bottom must be non-negative."
    @assert p.z0_top ≥ 0 "z0_top must be non-negative."
    p.wall_length in (:parabolic, :triangular, :infinite_depth) ||
        throw(ArgumentError(
            "wall_length must be :parabolic, :triangular, or :infinite_depth",
        ))
    return nothing
end

function validate_background_diffusivities(K0)
    @assert K0.momentum ≥ 0 "background momentum diffusivity must be non-negative."
    @assert K0.tracer ≥ 0 "background tracer diffusivity must be non-negative."
    @assert K0.tke ≥ 0 "background TKE diffusivity must be non-negative."
    @assert K0.q2l ≥ 0 "background q2l diffusivity must be non-negative."
    return nothing
end
