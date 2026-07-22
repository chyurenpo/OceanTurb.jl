# Mean-flow equations.
@inline RU(m, i) = @inbounds  m.constants.f * m.solution.V[i]
@inline RV(m, i) = @inbounds -m.constants.f * m.solution.U[i]
@inline RT(m, i) = 0
@inline RS(m, i) = 0

@inline shear_production(m, i) = @inbounds m.state.P[i]
@inline buoyancy_production(m, i) = @inbounds m.state.B[i]
@inline dissipation(m, i) = @inbounds m.state.ε[i]

"""
Explicit source for the TKE equation.

GOTM splitting:
- if P+B > 0, put P+B in the explicit source;
- otherwise put P in the explicit source and treat -B implicitly.
"""
@inline function Re(m, i)
    P = shear_production(m, i)
    B = buoyancy_production(m, i)
    return P + B > 0 ? P + B : P
end

"""
Positive implicit sink coefficient in OceanTurb's convention
`∂t e = ... + Re - Le*e`.
"""
@inline function Le(m, i)
    P = shear_production(m, i)
    B = buoyancy_production(m, i)
    ε = dissipation(m, i)
    @inbounds e = m.solution.e[i]
    e = max(e, m.parameters.k_min)

    return P + B > 0 ? ε/e : (ε - B)/e
end

@inline function q2l_terms(m, i)
    p = m.parameters

    @inbounds begin
        ℓ = m.state.ℓ[i]
        Lz = m.state.Lz[i]
        e = m.solution.e[i]
        P = m.state.P[i]
        B = m.state.B[i]
    end

    production = ℓ * p.E1 * P
    buoyancy = p.E3 * ℓ * B

    q3 = (2 * max(e, p.k_min))^(3/2)
    destruction = q3 / p.B1 * (1 + p.E2 * (ℓ/Lz)^2)

    return production, buoyancy, destruction
end

"""
Explicit source for the prognostic `q2l = q²ℓ` equation.
"""
@inline function Rq2l(m, i)
    production, buoyancy, _ = q2l_terms(m, i)
    return production + buoyancy > 0 ?
           production + buoyancy :
           production
end

"""
Positive implicit sink coefficient in
`∂t q2l = ... + Rq2l - Lq2l*q2l`.
"""
@inline function Lq2l(m, i)
    production, buoyancy, destruction = q2l_terms(m, i)
    @inbounds q2l = m.solution.q2l[i]
    q2l = max(q2l, m.parameters.q2l_min)

    return production + buoyancy > 0 ?
           destruction/q2l :
           (destruction - buoyancy)/q2l
end
