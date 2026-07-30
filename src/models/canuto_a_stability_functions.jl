"""
Coefficients for the Canuto et al. (2001) version-A second-moment
closure in the notation used by GOTM.

The decimal values intentionally reproduce GOTM's `CHCD01A` parameter
set, including the promotion of default-real Fortran literals to
GOTM's double-precision `REALTYPE`.
"""
@inline function canuto_a_coefficients(::Type{T}=Float64) where T
    reference_real(x) = T(Float32(x))

    cc1 = reference_real(5.0000)
    cc2 = reference_real(0.8000)
    cc3 = reference_real(1.9680)
    cc4 = reference_real(1.1360)
    cc5 = reference_real(0.0000)
    cc6 = reference_real(0.4000)

    ct1 = reference_real(5.9500)
    ct2 = reference_real(0.6000)
    ct3 = reference_real(1.0000)
    ct4 = reference_real(0.0000)
    ct5 = reference_real(0.3333)
    ctt = reference_real(0.7200)

    # In GOTM, `2./3.` is evaluated in default-real precision before
    # assignment to the double-precision variable `a1`.
    two_thirds = T(Float32(2) / Float32(3))
    a1  = two_thirds - cc2 / T(2)
    a2  = one(T) - cc3 / T(2)
    a3  = one(T) - cc4 / T(2)
    a4  = cc5 / T(2)
    a5  = one(T) / T(2) - cc6 / T(2)

    at1 = one(T) - ct2
    at2 = one(T) - ct3
    at3 = T(2) * (one(T) - ct4)
    at4 = T(2) * (one(T) - ct5)
    at5 = T(2) * ctt * (one(T) - ct5)

    return (
        cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, cc5=cc5, cc6=cc6,
        ct1=ct1, ct2=ct2, ct3=ct3, ct4=ct4, ct5=ct5, ctt=ctt,
        a1=a1, a2=a2, a3=a3, a4=a4, a5=a5,
        at1=at1, at2=at2, at3=at3, at4=at4, at5=at5,
    )
end

"""
Evaluate GOTM's quasi-equilibrium Canuto-A stability functions for
`alpha_N`. GOTM uses this branch internally when deriving the stable
`c_omega3`, even when the prognostic run uses the weak-equilibrium
second-moment closure.
"""
@inline function canuto_a_quasi_equilibrium(
    alpha_N,
    c_mu0,
)
    T = promote_type(typeof(alpha_N), typeof(c_mu0))
    c = canuto_a_coefficients(T)

    N  = c.cc1 / T(2)
    Nt = c.ct1

    d0 = T(36) * N^3 * Nt^2
    d1 = T(84) * c.a5 * c.at3 * N^2 * Nt +
         T(36) * c.at5 * N^3 * Nt
    d2 = T(9) * (c.at2^2 - c.at1^2) * N^3 -
         T(12) * (c.a2^2 - T(3) * c.a3^2) * N * Nt^2
    d3 = T(12) * c.a5 * c.at3 *
         (c.a2 * c.at1 - T(3) * c.a3 * c.at2) * N +
         T(12) * c.a5 * c.at3 *
         (c.a3^2 - c.a2^2) * Nt +
         T(12) * c.at5 *
         (T(3) * c.a3^2 - c.a2^2) * N * Nt
    d4 = T(48) * c.a5^2 * c.at3^2 * N +
         T(36) * c.a5 * c.at3 * c.at5 * N^2
    d5 = T(3) *
         (c.a2^2 - T(3) * c.a3^2) *
         (c.at1^2 - c.at2^2) * N

    n0 = T(36) * c.a1 * N^2 * Nt^2
    n1 = -T(12) * c.a5 * c.at3 *
         (c.at1 + c.at2) * N^2 +
         T(8) * c.a5 * c.at3 *
         (T(6) * c.a1 - c.a2 - T(3) * c.a3) * N * Nt +
         T(36) * c.a1 * c.at5 * N^2 * Nt
    n2 = T(9) * c.a1 *
         (c.at2^2 - c.at1^2) * N^2

    nt0 = T(12) * c.at3 * N^3 * Nt
    nt1 = T(12) * c.a5 * c.at3^2 * N^2
    nt2 = T(9) * c.a1 * c.at3 *
          (c.at1 - c.at2) * N^2 +
          (
              T(6) * c.a1 * (c.a2 - T(3) * c.a3) -
              T(4) * (c.a2^2 - T(3) * c.a3^2)
          ) * c.at3 * N * Nt

    an_min_numerator =
        -(d1 + nt0) +
        sqrt((d1 + nt0)^2 - T(4) * d0 * (d4 + nt1))
    an_min_denominator = T(2) * (d4 + nt1)
    an = max(
        T(alpha_N),
        T(0.5) * an_min_numerator / an_min_denominator,
    )

    tmp0 =
        -d0 - (d1 + nt0) * an -
        (d4 + nt1) * an^2
    tmp1 =
        -d2 + n0 +
        (n1 - d3 - nt2) * an
    tmp2 = n2 - d5

    alpha_M = if abs(tmp2) < T(1e-10)
        -tmp0 / tmp1
    else
        (
            -tmp1 +
            sqrt(tmp1^2 - T(4) * tmp0 * tmp2)
        ) / (T(2) * tmp2)
    end

    denominator =
        d0 + d1 * an + d2 * alpha_M +
        d3 * an * alpha_M + d4 * an^2 +
        d5 * alpha_M^2
    momentum_numerator =
        n0 + n1 * an + n2 * alpha_M
    tracer_numerator =
        nt0 + nt1 * an + nt2 * alpha_M
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
Return GOTM's log-layer and shear-free stability-function values for
the weak-equilibrium Canuto-A closure.
"""
@inline function canuto_a_reference_values(::Type{T}=Float64) where T
    c = canuto_a_coefficients(T)
    N = c.cc1 / T(2)

    c_mu0 = (
        (
            c.a2^2 - T(3) * c.a3^2 +
            T(3) * c.a1 * N
        ) / (T(3) * N^2)
    )^(T(1) / T(4))

    c_mu_shear_free = c.a1 / (N * c_mu0^3)
    return c_mu0, c_mu_shear_free
end

"""
Evaluate GOTM's local weak-equilibrium Canuto-A stability functions.

`alpha_M = (k / epsilon)^2 M^2` and
`alpha_N = (k / epsilon)^2 N^2`.
"""
@inline function canuto_a_weak_equilibrium(
    alpha_M,
    alpha_N,
    c_mu0,
)
    T = promote_type(
        typeof(alpha_M),
        typeof(alpha_N),
        typeof(c_mu0),
    )
    c = canuto_a_coefficients(T)

    N  = c.cc1 / T(2)
    Nt = c.ct1

    d0 = T(36) * N^3 * Nt^2
    d1 = T(84) * c.a5 * c.at3 * N^2 * Nt +
         T(36) * c.at5 * N^3 * Nt
    d2 = T(9) * (c.at2^2 - c.at1^2) * N^3 -
         T(12) * (c.a2^2 - T(3) * c.a3^2) * N * Nt^2
    d3 = T(12) * c.a5 * c.at3 *
         (c.a2 * c.at1 - T(3) * c.a3 * c.at2) * N +
         T(12) * c.a5 * c.at3 *
         (c.a3^2 - c.a2^2) * Nt +
         T(12) * c.at5 *
         (T(3) * c.a3^2 - c.a2^2) * N * Nt
    d4 = T(48) * c.a5^2 * c.at3^2 * N +
         T(36) * c.a5 * c.at3 * c.at5 * N^2
    d5 = T(3) *
         (c.a2^2 - T(3) * c.a3^2) *
         (c.at1^2 - c.at2^2) * N

    n0 = T(36) * c.a1 * N^2 * Nt^2
    n1 = -T(12) * c.a5 * c.at3 *
         (c.at1 + c.at2) * N^2 +
         T(8) * c.a5 * c.at3 *
         (T(6) * c.a1 - c.a2 - T(3) * c.a3) * N * Nt +
         T(36) * c.a1 * c.at5 * N^2 * Nt
    n2 = T(9) * c.a1 * (c.at2^2 - c.at1^2) * N^2

    nt0 = T(12) * c.at3 * N^3 * Nt
    nt1 = T(12) * c.a5 * c.at3^2 * N^2
    nt2 = T(9) * c.a1 * c.at3 *
          (c.at1 - c.at2) * N^2 +
          (
              T(6) * c.a1 * (c.a2 - T(3) * c.a3) -
              T(4) * (c.a2^2 - T(3) * c.a3^2)
          ) * c.at3 * N * Nt

    # GOTM first clips alpha_M in alpha_mnb.F90, then applies the
    # weak-equilibrium realizability limits in cmue_c.F90.
    as = max(T(alpha_M), T(1e-10))

    an_min_numerator =
        -(d1 + nt0) +
        sqrt((d1 + nt0)^2 - T(4) * d0 * (d4 + nt1))
    an_min_denominator = T(2) * (d4 + nt1)
    an_min = an_min_numerator / an_min_denominator
    an = max(T(alpha_N), T(0.5) * an_min)

    as_max_numerator =
        d0 * n0 +
        (d0 * n1 + d1 * n0) * an +
        (d1 * n1 + d4 * n0) * an^2 +
        d4 * n1 * an^3
    as_max_denominator =
        d2 * n0 +
        (d2 * n1 + d3 * n0) * an +
        d3 * n1 * an^2
    as = min(as, as_max_numerator / as_max_denominator)

    denominator =
        d0 + d1 * an + d2 * as +
        d3 * an * as + d4 * an^2 + d5 * as^2
    momentum_numerator = n0 + n1 * an + n2 * as
    tracer_numerator = nt0 + nt1 * an + nt2 * as
    inverse_c_mu0_cubed = inv(T(c_mu0)^3)

    c_mu = inverse_c_mu0_cubed *
           momentum_numerator / denominator
    c_mu_prime = inverse_c_mu0_cubed *
                 tracer_numerator / denominator

    return c_mu, c_mu_prime
end

"""
Compute the stable second-equation buoyancy coefficient used by GOTM
for a target steady-state gradient Richardson number.
"""
function canuto_a_compute_c3_stable(
    c1,
    c2,
    Ri_st,
    c_mu0,
)
    T = promote_type(
        typeof(c1),
        typeof(c2),
        typeof(Ri_st),
        typeof(c_mu0),
    )
    finite_difference = T(1e-8)
    tolerance = T(1e-10)
    alpha_N = T(5)
    converged = false

    residual(an) = begin
        c_mu, c_mu_prime =
            canuto_a_quasi_equilibrium(
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
        step = -f / ((fp - f) / finite_difference)
        alpha_N += T(0.5) * step

        abs(step) > T(100) &&
            error("GOTM c_omega3 iteration did not converge.")

        if abs(step) < tolerance
            converged = true
            break
        end
    end

    converged ||
        error("GOTM c_omega3 iteration exceeded 101 iterations.")

    c_mu, c_mu_prime =
        canuto_a_quasi_equilibrium(
            alpha_N,
            T(c_mu0),
        )

    return T(c2) +
           (T(c1) - T(c2)) / T(Ri_st) *
           c_mu / c_mu_prime
end
