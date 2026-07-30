# Mean-flow equations.
@inline RU(m, i) = @inbounds  m.constants.f * m.solution.V[i] + m.forcing.U(m, i)
@inline RV(m, i) = @inbounds -m.constants.f * m.solution.U[i] + m.forcing.V(m, i)
@inline RT(m, i) = m.forcing.T(m, i)
@inline RS(m, i) = m.forcing.S(m, i)

@inline shear_production(m, i) = @inbounds m.state.P[i]
@inline buoyancy_production(m, i) = @inbounds m.state.B[i]
@inline dissipation(m, i) = @inbounds m.state.epsilon[i]

# Positive production and destabilizing buoyancy flux are explicit.
@inline function Re(m, i)
    P = shear_production(m, i)
    B = buoyancy_production(m, i)
    return P + max(B, zero(B)) + m.forcing.e(m, i)
end

# Stable buoyancy destruction, physical dissipation, and the maintained
# CATKE implementation's near-bottom TKE flux are positive implicit sinks.
@inline function Le(m, i)
    @inbounds e = max(m.solution.e[i], m.parameters.minimum_tke)
    B = buoyancy_production(m, i)
    physical_dissipation = sqrt(e) / max(
        m.state.ell_D[i],
        m.parameters.minimum_length,
    )
    buoyancy_destruction =
        e > m.parameters.minimum_tke ?
        max(-B, zero(B)) / e :
        zero(e)

    bottom_dissipation = i == 1 ?
        m.tke_equation.C_bottom_dissipation *
        sqrt(e) / Δf(m.grid, 1) :
        zero(e)

    return physical_dissipation +
           buoyancy_destruction +
           bottom_dissipation
end
