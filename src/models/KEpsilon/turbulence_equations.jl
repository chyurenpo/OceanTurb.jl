# Mean-flow equations.
@inline RU(m, i) = @inbounds  m.constants.f * m.solution.V[i] + m.forcing.U(m, i)
@inline RV(m, i) = @inbounds -m.constants.f * m.solution.U[i] + m.forcing.V(m, i)
@inline RT(m, i) = m.forcing.T(m, i)
@inline RS(m, i) = m.forcing.S(m, i)

@inline shear_production(m, i) = @inbounds m.state.P[i]
@inline buoyancy_production(m, i) = @inbounds m.state.B[i]
@inline dissipation(m, i) = @inbounds m.solution.epsilon[i]
@inline turbulent_length_scale(m, i) = @inbounds m.state.ell[i]

# Explicit TKE source. If P+B is negative, stable buoyancy destruction
# is moved to the positive implicit sink coefficient Lk.
@inline function Rk(m, i)
    P = shear_production(m, i)
    B = buoyancy_production(m, i)
    return P + B > 0 ? P + B : P
end

# Positive implicit sink coefficient in
# partial_t k = ... + Rk - Lk*k.
@inline function Lk(m, i)
    P = shear_production(m, i)
    B = buoyancy_production(m, i)
    epsilon = dissipation(m, i)
    @inbounds k = max(m.solution.k[i], m.parameters.k_min)

    return P + B > 0 ? epsilon/k : (epsilon - B)/k
end

@inline function epsilon_terms(m, i)
    p = m.parameters

    @inbounds begin
        k = max(m.solution.k[i], p.k_min)
        epsilon = max(m.solution.epsilon[i], p.epsilon_min)
        P = m.state.P[i]
        B = m.state.B[i]
    end

    c_epsilon3 = B > 0 ?
                 p.c_epsilon3_unstable :
                 p.c_epsilon3_stable

    epsilon_over_k = epsilon / k
    production = epsilon_over_k * p.c_epsilon1 * P
    buoyancy = epsilon_over_k * c_epsilon3 * B
    destruction = epsilon_over_k * p.c_epsilon2 * epsilon

    return production, buoyancy, destruction
end

# Explicit source for the dissipation-rate equation.
@inline function Repsilon(m, i)
    production, buoyancy, _ = epsilon_terms(m, i)
    return production + buoyancy > 0 ?
           production + buoyancy :
           production
end

# Positive implicit sink coefficient in
# partial_t epsilon = ... + Repsilon - Lepsilon*epsilon.
@inline function Lepsilon(m, i)
    production, buoyancy, destruction = epsilon_terms(m, i)
    @inbounds epsilon = max(
        m.solution.epsilon[i],
        m.parameters.epsilon_min,
    )

    return production + buoyancy > 0 ?
           destruction/epsilon :
           (destruction - buoyancy)/epsilon
end
