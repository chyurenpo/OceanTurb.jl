"""
Classical Mellor–Yamada level-2.5 constants.

The defaults follow the values exposed by GOTM for the MY q²/2 and q²ℓ
equations: B1=16.6, Sq=Sl=0.2, E1=1.8, E2=1.33, E3=1.8.
"""
Base.@kwdef struct MY25Parameters{T} <: AbstractParameters
    B1          :: T = 16.6
    Sq          :: T = 0.2
    Sl          :: T = 0.2
    E1          :: T = 1.8
    E2          :: T = 1.33
    E3          :: T = 1.8

    κ           :: T = 0.4
    galperin    :: T = 0.27
    length_limit :: Bool = true

    k_min       :: T = 1e-10
    eps_min     :: T = 1e-12
    q2l_min     :: T = 1e-14

    z0_bottom   :: T = 1e-4
    z0_top      :: T = 1e-4

    # :parabolic, :triangular, or :infinite_depth
    wall_length :: Symbol = :parabolic

    # Numerical guards for the algebraic stability functions.
    stability_floor :: T = 0.0
    denominator_floor :: T = 1e-14
end

"""
MY82 pressure-strain and scalar-flux coefficients used by GOTM's
second-order closure option.
"""
Base.@kwdef struct MY82Coefficients{T} <: AbstractParameters
    cc1 :: T = 6.0
    cc2 :: T = 0.32
    cc3 :: T = 0.0
    cc4 :: T = 0.0
    cc5 :: T = 0.0
    cc6 :: T = 0.0

    ct1 :: T = 3.728
    ct2 :: T = 0.0
    ct3 :: T = 0.0
    ct4 :: T = 0.0
    ct5 :: T = 0.0
    ctt :: T = 0.6102
end

"""
Background molecular/numerical diffusivities [m² s⁻¹].
"""
Base.@kwdef struct BackgroundDiffusivities{T} <: AbstractParameters
    momentum :: T = 1e-6
    tracer   :: T = 1e-7
    tke      :: T = 1e-7
    q2l      :: T = 1e-7
end

@inline cde(p) = 2 * sqrt(2) / p.B1
@inline cm0(p) = cbrt(cde(p))
@inline minimum_length(p) = cde(p) * p.k_min^(3/2) / p.eps_min
