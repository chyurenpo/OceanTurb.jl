"""
    friction_velocity(momentum_flux_u, momentum_flux_v)

Return `u_tau = (Qu^2 + Qv^2)^(1/4)` when `Qu` and `Qv` are
kinematic momentum fluxes.
"""
@inline friction_velocity(momentum_flux_u, momentum_flux_v) =
    (momentum_flux_u^2 + momentum_flux_v^2)^(1/4)

"""
GOTM logarithmic-wall Dirichlet value for turbulent kinetic energy.
"""
@inline logarithmic_tke_value(u_tau, p::KEpsilonParameters) =
    u_tau^2 / p.c_mu0^2

"""
GOTM logarithmic-wall Dirichlet value for turbulent kinetic energy
dissipation. `distance` is the distance represented by the boundary
stencil.
"""
@inline logarithmic_epsilon_value(
    k,
    distance,
    z0,
    p::KEpsilonParameters,
) =
    p.c_mu0^3 * max(k, p.k_min)^(3/2) /
    (
        p.kappa *
        max(distance + z0, p.length_min)
    )

"""
    logarithmic_epsilon_flux(k, distance, z0, parameters)

Return the inward logarithmic-wall flux used by GOTM for a Neumann
boundary condition on `epsilon`.
"""
@inline logarithmic_epsilon_flux(
    k,
    distance,
    z0,
    p::KEpsilonParameters,
) =
    p.c_mu0^4 * max(k, p.k_min)^2 /
    (
        p.sigma_epsilon *
        max(distance + z0, p.length_min)
    )

@inline evaluate_wall_input(value::Number, model) = value
@inline evaluate_wall_input(func, model) = func(model)

"""
    BottomLogarithmicWallConditions(u_tau; parameters)

Return `(k_bc, epsilon_bc)` value boundary conditions for a logarithmic
bottom wall. `u_tau` may be a number or a callable `u_tau(model)`.
"""
function BottomLogarithmicWallConditions(
    u_tau;
    parameters = KEpsilonParameters(),
)
    k_bc = ValueBoundaryCondition(
        m -> logarithmic_tke_value(
            evaluate_wall_input(u_tau, m),
            parameters,
        ),
    )

    epsilon_bc = ValueBoundaryCondition(
        m -> logarithmic_epsilon_value(
            m.solution.k[1],
            Δf(m.grid, 1),
            parameters.z0_bottom,
            parameters,
        ),
    )

    return k_bc, epsilon_bc
end

"""
    TopLogarithmicWallConditions(u_tau; parameters)

Return `(k_bc, epsilon_bc)` value boundary conditions for a logarithmic
upper wall.
"""
function TopLogarithmicWallConditions(
    u_tau;
    parameters = KEpsilonParameters(),
)
    k_bc = ValueBoundaryCondition(
        m -> logarithmic_tke_value(
            evaluate_wall_input(u_tau, m),
            parameters,
        ),
    )

    epsilon_bc = ValueBoundaryCondition(
        m -> logarithmic_epsilon_value(
            m.solution.k[m.grid.N],
            Δf(m.grid, m.grid.N),
            parameters.z0_top,
            parameters,
        ),
    )

    return k_bc, epsilon_bc
end

"""
    LogarithmicNeumannBoundaryConditions([FT=Float64]; parameters)

Return `(k_bcs, epsilon_bcs)` matching GOTM's logarithmic Neumann
boundary conditions. GOTM prescribes zero boundary flux for `k` and a
positive inward flux
`c_mu0^4 * k^2 / (sigma_epsilon * (Delta z / 2 + z0))` for `epsilon`
at both walls. OceanTurb's top-flux sign convention is
outward-positive, so the top condition carries a minus sign.
"""
function LogarithmicNeumannBoundaryConditions(
    FT=Float64;
    parameters=KEpsilonParameters(),
)
    k_bcs = ZeroFluxBoundaryConditions(FT)
    epsilon_bcs = OceanTurb.BoundaryConditions(
        FT;
        bottom=FluxBoundaryCondition(
            m -> logarithmic_epsilon_flux(
                m.solution.k[1],
                Δf(m.grid, 1) / 2,
                parameters.z0_bottom,
                parameters,
            ),
        ),
        top=FluxBoundaryCondition(
            m -> -logarithmic_epsilon_flux(
                m.solution.k[m.grid.N],
                Δf(m.grid, m.grid.N) / 2,
                parameters.z0_top,
                parameters,
            ),
        ),
    )

    return k_bcs, epsilon_bcs
end
