@inline function kantha_clayson_derived_coefficients(
    c::KanthaClaysonCoefficients{T},
) where T
    # GOTM evaluates 2./3. in default-real precision before promoting it.
    two_thirds = T(Float32(2) / Float32(3))

    a1  = two_thirds - c.cc2 / T(2)
    a2  = one(T) - c.cc3 / T(2)
    a3  = one(T) - c.cc4 / T(2)
    a4  = c.cc5 / T(2)
    a5  = one(T) / T(2) - c.cc6 / T(2)

    at1 = one(T) - c.ct2
    at2 = one(T) - c.ct3
    at3 = T(2) * (one(T) - c.ct4)
    at4 = T(2) * (one(T) - c.ct5)
    at5 = T(2) * c.ctt * (one(T) - c.ct5)

    return (
        a1=a1, a2=a2, a3=a3, a4=a4, a5=a5,
        at1=at1, at2=at2, at3=at3, at4=at4, at5=at5,
    )
end

@inline function kantha_clayson_polynomials(
    ::Type{T}=Float64,
) where T
    c = KanthaClaysonCoefficients(T)
    a = kantha_clayson_derived_coefficients(c)

    N  = c.cc1 / T(2)
    Nt = c.ct1

    d0 = T(36) * N^3 * Nt^2
    d1 = T(84) * a.a5 * a.at3 * N^2 * Nt +
         T(36) * a.at5 * N^3 * Nt
    d2 = T(9) * (a.at2^2 - a.at1^2) * N^3 -
         T(12) * (a.a2^2 - T(3) * a.a3^2) * N * Nt^2
    d3 = T(12) * a.a5 * a.at3 *
         (a.a2 * a.at1 - T(3) * a.a3 * a.at2) * N +
         T(12) * a.a5 * a.at3 *
         (a.a3^2 - a.a2^2) * Nt +
         T(12) * a.at5 *
         (T(3) * a.a3^2 - a.a2^2) * N * Nt
    d4 = T(48) * a.a5^2 * a.at3^2 * N +
         T(36) * a.a5 * a.at3 * a.at5 * N^2
    d5 = T(3) *
         (a.a2^2 - T(3) * a.a3^2) *
         (a.at1^2 - a.at2^2) * N

    n0 = T(36) * a.a1 * N^2 * Nt^2
    n1 = -T(12) * a.a5 * a.at3 *
         (a.at1 + a.at2) * N^2 +
         T(8) * a.a5 * a.at3 *
         (T(6) * a.a1 - a.a2 - T(3) * a.a3) * N * Nt +
         T(36) * a.a1 * a.at5 * N^2 * Nt
    n2 = T(9) * a.a1 *
         (a.at2^2 - a.at1^2) * N^2

    nt0 = T(12) * a.at3 * N^3 * Nt
    nt1 = T(12) * a.a5 * a.at3^2 * N^2
    nt2 = T(9) * a.a1 * a.at3 *
          (a.at1 - a.at2) * N^2 +
          (
              T(6) * a.a1 *
              (a.a2 - T(3) * a.a3) -
              T(4) *
              (a.a2^2 - T(3) * a.a3^2)
          ) * a.at3 * N * Nt

    return (
        coefficients=c,
        derived=a,
        N=N,
        Nt=Nt,
        d0=d0, d1=d1, d2=d2, d3=d3, d4=d4, d5=d5,
        n0=n0, n1=n1, n2=n2,
        nt0=nt0, nt1=nt1, nt2=nt2,
    )
end

"""
Return GOTM's log-layer and shear-free stability-function values for
the Kantha-Clayson weak-equilibrium closure.
"""
@inline function kantha_clayson_reference_values(
    ::Type{T}=Float64,
) where T
    p = kantha_clayson_polynomials(T)
    a = p.derived

    c_mu0 = (
        (
            a.a2^2 - T(3) * a.a3^2 +
            T(3) * a.a1 * p.N
        ) / (T(3) * p.N^2)
    )^(T(1) / T(4))

    c_mu_shear_free =
        a.a1 / (p.N * c_mu0^3)
    return c_mu0, c_mu_shear_free
end

"""
Evaluate GOTM's quasi-equilibrium Kantha-Clayson stability functions.

GOTM uses this branch when deriving `E3` from the target steady-state
Richardson number, even though the prognostic MY2.5 calculation uses
the weak-equilibrium branch.
"""
@inline function kantha_clayson_quasi_equilibrium(
    alpha_N,
    c_mu0,
)
    T = promote_type(typeof(alpha_N), typeof(c_mu0))
    p = kantha_clayson_polynomials(T)

    an_min_numerator =
        -(p.d1 + p.nt0) +
        sqrt(
            (p.d1 + p.nt0)^2 -
            T(4) * p.d0 * (p.d4 + p.nt1),
        )
    an_min_denominator =
        T(2) * (p.d4 + p.nt1)
    an = max(
        T(alpha_N),
        T(0.5) *
        an_min_numerator / an_min_denominator,
    )

    tmp0 =
        -p.d0 - (p.d1 + p.nt0) * an -
        (p.d4 + p.nt1) * an^2
    tmp1 =
        -p.d2 + p.n0 +
        (p.n1 - p.d3 - p.nt2) * an
    tmp2 = p.n2 - p.d5

    alpha_M = if abs(tmp2) < T(1e-10)
        -tmp0 / tmp1
    else
        (
            -tmp1 +
            sqrt(tmp1^2 - T(4) * tmp0 * tmp2)
        ) / (T(2) * tmp2)
    end

    denominator =
        p.d0 + p.d1 * an + p.d2 * alpha_M +
        p.d3 * an * alpha_M + p.d4 * an^2 +
        p.d5 * alpha_M^2
    momentum_numerator =
        p.n0 + p.n1 * an + p.n2 * alpha_M
    tracer_numerator =
        p.nt0 + p.nt1 * an + p.nt2 * alpha_M
    inverse_c_mu0_cubed = inv(T(c_mu0)^3)

    c_mu =
        inverse_c_mu0_cubed *
        momentum_numerator / denominator
    c_mu_prime =
        inverse_c_mu0_cubed *
        tracer_numerator / denominator

    return c_mu, c_mu_prime
end

"""
Evaluate GOTM's local weak-equilibrium Kantha-Clayson stability
functions.

`alpha_M = (k / epsilon)^2 M²` and
`alpha_N = (k / epsilon)^2 N²`. The returned `alpha_M` and `alpha_N`
include GOTM's realizability limits.
"""
@inline function kantha_clayson_weak_equilibrium(
    alpha_M,
    alpha_N,
    c_mu0,
)
    T = promote_type(
        typeof(alpha_M),
        typeof(alpha_N),
        typeof(c_mu0),
    )
    p = kantha_clayson_polynomials(T)

    # alpha_mnb.F90 applies this floor before cmue_c.F90 is called.
    alpha_M_floor = T(Float32(1e-10))
    as = max(T(alpha_M), alpha_M_floor)

    an_min_numerator =
        -(p.d1 + p.nt0) +
        sqrt(
            (p.d1 + p.nt0)^2 -
            T(4) * p.d0 * (p.d4 + p.nt1),
        )
    an_min_denominator =
        T(2) * (p.d4 + p.nt1)
    an_min = an_min_numerator / an_min_denominator
    an = max(T(alpha_N), T(0.5) * an_min)

    as_max_numerator =
        p.d0 * p.n0 +
        (p.d0 * p.n1 + p.d1 * p.n0) * an +
        (p.d1 * p.n1 + p.d4 * p.n0) * an^2 +
        p.d4 * p.n1 * an^3
    as_max_denominator =
        p.d2 * p.n0 +
        (p.d2 * p.n1 + p.d3 * p.n0) * an +
        p.d3 * p.n1 * an^2
    as = min(as, as_max_numerator / as_max_denominator)

    denominator =
        p.d0 + p.d1 * an + p.d2 * as +
        p.d3 * an * as + p.d4 * an^2 +
        p.d5 * as^2
    momentum_numerator =
        p.n0 + p.n1 * an + p.n2 * as
    tracer_numerator =
        p.nt0 + p.nt1 * an + p.nt2 * as
    inverse_c_mu0_cubed = inv(T(c_mu0)^3)

    c_mu =
        inverse_c_mu0_cubed *
        momentum_numerator / denominator
    c_mu_prime =
        inverse_c_mu0_cubed *
        tracer_numerator / denominator

    return as, an, c_mu, c_mu_prime
end

"""
Derive the MY2.5 buoyancy coefficient `E3` using GOTM's
`compute_cpsi3` iteration.
"""
function kantha_clayson_compute_my25_E3(
    E1,
    Ri_st,
    c_mu0,
)
    T = promote_type(
        typeof(E1),
        typeof(Ri_st),
        typeof(c_mu0),
    )
    E2_equilibrium = one(T)
    finite_difference = T(1e-8)
    tolerance = T(1e-10)
    alpha_N = T(5)
    converged = false

    residual(an) = begin
        c_mu, c_mu_prime =
            kantha_clayson_quasi_equilibrium(
                an,
                T(c_mu0),
            )
        c_mu * an / T(Ri_st) -
        c_mu_prime * an -
        inv(T(c_mu0)^3)
    end

    for _ in 0:100
        f = residual(alpha_N)
        fp = residual(alpha_N + finite_difference)
        step =
            -f / ((fp - f) / finite_difference)
        alpha_N += T(0.5) * step

        abs(step) > T(100) &&
            error("GOTM MY2.5 E3 iteration did not converge.")

        if abs(step) < tolerance
            converged = true
            break
        end
    end

    converged ||
        error("GOTM MY2.5 E3 iteration exceeded 101 iterations.")

    c_mu, c_mu_prime =
        kantha_clayson_quasi_equilibrium(
            alpha_N,
            T(c_mu0),
        )

    return E2_equilibrium +
           (T(E1) - E2_equilibrium) / T(Ri_st) *
           c_mu / c_mu_prime
end

"""
Return the dependent coefficients for GOTM's canonical
Kantha-Clayson MY2.5 configuration.
"""
function kantha_clayson_my25_coefficients(;
    Ri_st = 0.23,
    Sq = 0.2,
    Sl = 0.2,
    E1 = 1.8,
    E2 = 1.33,
)
    T = promote_type(
        typeof(Ri_st),
        typeof(Sq),
        typeof(Sl),
        typeof(E1),
        typeof(E2),
    )
    c_mu0, c_mu_shear_free =
        kantha_clayson_reference_values(T)
    cde_value = c_mu0^3

    # GOTM evaluates 2.**1.5 in default-real precision.
    two_to_three_halves =
        T(Float32(2)^Float32(1.5))
    B1 = two_to_three_halves / cde_value
    kappa = sqrt(
        (T(E2) - T(E1) + one(T)) /
        (T(Sl) * B1),
    )
    E3 = kantha_clayson_compute_my25_E3(
        T(E1),
        T(Ri_st),
        c_mu0,
    )

    return (
        c_mu0=c_mu0,
        c_mu_shear_free=c_mu_shear_free,
        cde=cde_value,
        B1=B1,
        kappa=kappa,
        E3=E3,
    )
end
