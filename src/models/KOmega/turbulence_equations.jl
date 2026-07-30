# Mean-flow equations.
@inline RU(m, i) = @inbounds  m.constants.f * m.solution.V[i] + m.forcing.U(m, i)
@inline RV(m, i) = @inbounds -m.constants.f * m.solution.U[i] + m.forcing.V(m, i)
@inline RT(m, i) = m.forcing.T(m, i)
@inline RS(m, i) = m.forcing.S(m, i)

@inline shear_production(m, i) = @inbounds m.state.P[i]
@inline buoyancy_production(m, i) = @inbounds m.state.B[i]
@inline dissipation(m, i) = @inbounds m.state.epsilon[i]
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

@inline function omega_terms(m, i)
    p = m.parameters

    @inbounds begin
        k = max(m.solution.k[i], p.k_min)
        omega = m.solution.omega[i]
        P = m.state.P[i]
        B = m.state.B[i]
        epsilon = m.state.epsilon[i]
    end

    c_omega3 = B > 0 ?
               p.c_omega3_unstable :
               p.c_omega3_stable

    omega_over_k = omega / k
    production = omega_over_k * p.c_omega1 * P
    buoyancy = omega_over_k * c_omega3 * B
    destruction = omega_over_k * p.c_omega2 * epsilon

    return production, buoyancy, destruction
end

# Explicit source for the inverse-time-scale equation.
@inline function Romega(m, i)
    production, buoyancy, _ = omega_terms(m, i)
    return production + buoyancy > 0 ?
           production + buoyancy :
           production
end

# Positive implicit sink coefficient in
# partial_t omega = ... + Romega - Lomega*omega.
@inline function Lomega(m, i)
    production, buoyancy, destruction = omega_terms(m, i)
    @inbounds omega = max(
        m.solution.omega[i],
        m.parameters.omega_min,
    )

    return production + buoyancy > 0 ?
           destruction/omega :
           (destruction - buoyancy)/omega
end
