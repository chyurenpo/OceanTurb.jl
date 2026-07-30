"""
Return the dependent coefficients for the GOTM-compatible,
weak-equilibrium Canuto-A `k`-`omega` configuration.
"""
function canuto_a_komega_coefficients(;
    Ri_st = 0.25,
    c_omega1 = 0.555,
    c_omega2 = 0.833,
    sigma_omega = 2.0,
)
    T = promote_type(
        typeof(Ri_st),
        typeof(c_omega1),
        typeof(c_omega2),
        typeof(sigma_omega),
    )
    c_mu0, c_mu_shear_free =
        canuto_a_reference_values(T)
    kappa = c_mu0 * sqrt(
        T(sigma_omega) *
        (T(c_omega2) - T(c_omega1)),
    )
    c_omega3_stable = canuto_a_compute_c3_stable(
        T(c_omega1),
        T(c_omega2),
        T(Ri_st),
        c_mu0,
    )

    return (
        c_mu0=c_mu0,
        c_mu_shear_free=c_mu_shear_free,
        kappa=kappa,
        c_omega3_stable=c_omega3_stable,
    )
end

"""
Parameters for OceanTurb's GOTM-compatible oceanic `k`-`omega`
closure.

The default dependent coefficients are derived from GOTM's Canuto-A
second-moment closure and the target steady-state Richardson number.
"""
struct KOmegaParameters{T} <: AbstractParameters
    Ri_st               :: T
    c_mu0               :: T

    sigma_k             :: T
    sigma_omega         :: T

    c_omega1            :: T
    c_omega2            :: T
    c_omega3_stable     :: T
    c_omega3_unstable   :: T

    kappa               :: T
    galperin            :: T
    length_limit        :: Bool

    k_min               :: T
    omega_min           :: T
    epsilon_min         :: T
    length_min          :: T

    z0_bottom           :: T
    z0_top              :: T
end

"""
    KOmegaParameters(; kwargs...)

Construct the canonical GOTM/Canuto-A `k`-`omega` parameters.
`c_mu0`, `kappa`, and `c_omega3_stable` are derived from the primary
model constants unless explicitly overridden.
"""
function KOmegaParameters(;
    Ri_st = 0.25,
    c_omega1 = 0.555,
    c_omega2 = 0.833,
    c_omega3_unstable = 0.5,
    sigma_k = 2.0,
    sigma_omega = 2.0,
    c_mu0 = nothing,
    c_omega3_stable = nothing,
    kappa = nothing,
    galperin = 0.27,
    length_limit = true,
    k_min = 1e-10,
    omega_min = 1e-12,
    epsilon_min = 1e-12,
    length_min = 1e-12,
    z0_bottom = 1.5e-3,
    z0_top = 2.0e-2,
)
    coefficients = canuto_a_komega_coefficients(
        Ri_st=Ri_st,
        c_omega1=c_omega1,
        c_omega2=c_omega2,
        sigma_omega=sigma_omega,
    )

    c_mu0 === nothing &&
        (c_mu0 = coefficients.c_mu0)
    c_omega3_stable === nothing &&
        (c_omega3_stable =
            coefficients.c_omega3_stable)
    kappa === nothing &&
        (kappa = coefficients.kappa)

    T = promote_type(
        typeof(Ri_st),
        typeof(c_mu0),
        typeof(sigma_k),
        typeof(sigma_omega),
        typeof(c_omega1),
        typeof(c_omega2),
        typeof(c_omega3_stable),
        typeof(c_omega3_unstable),
        typeof(kappa),
        typeof(galperin),
        typeof(k_min),
        typeof(omega_min),
        typeof(epsilon_min),
        typeof(length_min),
        typeof(z0_bottom),
        typeof(z0_top),
    )

    return KOmegaParameters{T}(
        T(Ri_st),
        T(c_mu0),
        T(sigma_k),
        T(sigma_omega),
        T(c_omega1),
        T(c_omega2),
        T(c_omega3_stable),
        T(c_omega3_unstable),
        T(kappa),
        T(galperin),
        length_limit,
        T(k_min),
        T(omega_min),
        T(epsilon_min),
        T(length_min),
        T(z0_bottom),
        T(z0_top),
    )
end

"""
Molecular/background diffusivities in square metres per second.

Defaults match the GOTM free-convection reference configuration. GOTM
does not add separate background diffusivities to the `k` and `omega`
equations.
"""
Base.@kwdef struct BackgroundDiffusivities{T} <: AbstractParameters
    momentum :: T = 1.3e-6
    tracer   :: T = 1.4e-7
    tke      :: T = 0.0
    omega    :: T = 0.0
end

function validate_parameters(p)
    @assert p.Ri_st > 0 "Ri_st must be positive."
    @assert p.c_mu0 > 0 "c_mu0 must be positive."
    @assert p.sigma_k > 0 "sigma_k must be positive."
    @assert p.sigma_omega > 0 "sigma_omega must be positive."
    @assert p.c_omega1 ≥ 0 "c_omega1 must be non-negative."
    @assert p.c_omega2 > 0 "c_omega2 must be positive."
    @assert isfinite(p.c_omega3_stable) "c_omega3_stable must be finite."
    @assert isfinite(p.c_omega3_unstable) "c_omega3_unstable must be finite."
    @assert p.kappa > 0 "kappa must be positive."
    @assert p.galperin > 0 "galperin must be positive."
    @assert p.k_min > 0 "k_min must be positive."
    @assert p.omega_min > 0 "omega_min must be positive."
    @assert p.epsilon_min > 0 "epsilon_min must be positive."
    @assert p.length_min > 0 "length_min must be positive."
    @assert p.z0_bottom ≥ 0 "z0_bottom must be non-negative."
    @assert p.z0_top ≥ 0 "z0_top must be non-negative."
    return nothing
end

function validate_background_diffusivities(K0)
    @assert K0.momentum ≥ 0 "background momentum diffusivity must be non-negative."
    @assert K0.tracer ≥ 0 "background tracer diffusivity must be non-negative."
    @assert K0.tke ≥ 0 "background TKE diffusivity must be non-negative."
    @assert K0.omega ≥ 0 "background omega diffusivity must be non-negative."
    return nothing
end
