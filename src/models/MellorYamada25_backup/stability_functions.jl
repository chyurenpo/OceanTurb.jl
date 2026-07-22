"""
    my82_quasi_equilibrium(alpha_N, parameters, coefficients)

Return `(alpha_M, alpha_N_limited, c_mu, c_mu_prime)` for the
Mellor–Yamada (1982) quasi-equilibrium second-order closure.

This is an independent Julia transcription of the algebra in GOTM's
`cmue_d.F90`; it is intentionally isolated from the prognostic
equations so it can be tested separately.
"""
@inline function my82_quasi_equilibrium(
    alpha_N,
    p::MY25Parameters,
    c::MY82Coefficients = MY82Coefficients(),
)
    # Derived second-order coefficients used in GOTM.
    a1  = 2/3 - c.cc2/2
    a2  = 1   - c.cc3/2
    a3  = 1   - c.cc4/2
    a4  = c.cc5/2
    a5  = 1/2 - c.cc6/2

    at1 = 1 - c.ct2
    at2 = 1 - c.ct3
    at3 = 2 * (1 - c.ct4)
    at4 = 2 * (1 - c.ct5)
    at5 = 2 * c.ctt * (1 - c.ct5)

    N  = 0.5 * c.cc1
    Nt = c.ct1

    d0 = 36 * N^3 * Nt^2
    d1 = 84 * a5 * at3 * N^2 * Nt + 36 * at5 * N^3 * Nt
    d2 = 9 * (at2^2 - at1^2) * N^3 -
         12 * (a2^2 - 3*a3^2) * N * Nt^2
    d3 = 12 * a5 * at3 * (a2*at1 - 3*a3*at2) * N +
         12 * a5 * at3 * (a3^2 - a2^2) * Nt +
         12 * at5 * (3*a3^2 - a2^2) * N * Nt
    d4 = 48 * a5^2 * at3^2 * N +
         36 * a5 * at3 * at5 * N^2
    d5 = 3 * (a2^2 - 3*a3^2) * (at1^2 - at2^2) * N

    n0 = 36 * a1 * N^2 * Nt^2
    n1 = -12 * a5 * at3 * (at1 + at2) * N^2 +
          8 * a5 * at3 * (6*a1 - a2 - 3*a3) * N * Nt +
         36 * a1 * at5 * N^2 * Nt
    n2 = 9 * a1 * (at2^2 - at1^2) * N^2

    nt0 = 12 * at3 * N^3 * Nt
    nt1 = 12 * a5 * at3^2 * N^2
    nt2 = 9 * a1 * at3 * (at1 - at2) * N^2 +
          (6*a1*(a2 - 3*a3) - 4*(a2^2 - 3*a3^2)) * at3 * N * Nt

    # GOTM's lower bound prevents negative equilibrium shear number
    # under sufficiently strong unstable stratification.
    discriminant_min = max(
        zero(alpha_N),
        (d1 + nt0)^2 - 4*d0*(d4 + nt1),
    )
    denominator_min = 2 * (d4 + nt1)
    alpha_N_min = if abs(denominator_min) > p.denominator_floor
        (-(d1 + nt0) + sqrt(discriminant_min)) / denominator_min
    else
        -typemax(typeof(alpha_N))
    end

    an = max(alpha_N, 0.5 * alpha_N_min)

    tmp0 = -d0 - (d1 + nt0)*an - (d4 + nt1)*an^2
    tmp1 = -d2 + n0 + (n1 - d3 - nt2)*an
    tmp2 = n2 - d5

    am = if abs(tmp2) < 1e-10
        -tmp0 / copysign(max(abs(tmp1), p.denominator_floor), tmp1)
    else
        discriminant = max(zero(an), tmp1^2 - 4*tmp0*tmp2)
        (-tmp1 + sqrt(discriminant)) / (2*tmp2)
    end

    am = max(am, zero(am))

    dcm  = d0 + d1*an + d2*am + d3*an*am + d4*an^2 + d5*am^2
    ncm  = n0 + n1*an + n2*am
    ncmp = nt0 + nt1*an + nt2*am

    safe_dcm = copysign(max(abs(dcm), p.denominator_floor), dcm)
    cm3_inv = inv(cm0(p)^3)

    c_mu       = cm3_inv * ncm  / safe_dcm
    c_mu_prime = cm3_inv * ncmp / safe_dcm

    c_mu = isfinite(c_mu) ? max(c_mu, p.stability_floor) : p.stability_floor
    c_mu_prime = isfinite(c_mu_prime) ?
                 max(c_mu_prime, p.stability_floor) :
                 p.stability_floor

    return am, an, c_mu, c_mu_prime
end
