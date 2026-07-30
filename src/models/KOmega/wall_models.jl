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
@inline logarithmic_tke_value(u_tau, p::KOmegaParameters) =
    u_tau^2 / p.c_mu0^2

"""
GOTM logarithmic-wall Dirichlet value for the inverse turbulent time
scale. `distance` is the distance represented by the boundary stencil.
"""
@inline logarithmic_omega_value(k, distance, z0, p::KOmegaParameters) =
    sqrt(max(k, p.k_min)) /
    (
        p.c_mu0 *
        p.kappa *
        max(distance + z0, p.length_min)
    )

"""
    logarithmic_omega_flux(k, distance, z0, parameters)

Return the inward logarithmic-wall flux used by GOTM for a Neumann
boundary condition on `omega`.
"""
@inline logarithmic_omega_flux(
    k,
    distance,
    z0,
    p::KOmegaParameters,
) =
    max(k, p.k_min) /
    (
        p.sigma_omega *
        max(distance + z0, p.length_min)
    )

@inline evaluate_wall_input(value::Number, model) = value
@inline evaluate_wall_input(func, model) = func(model)

"""
    BottomLogarithmicWallConditions(u_tau; parameters)

Return `(k_bc, omega_bc)` value boundary conditions for a logarithmic
bottom wall. `u_tau` may be a number or a callable `u_tau(model)`.
"""
function BottomLogarithmicWallConditions(
    u_tau;
    parameters = KOmegaParameters(),
)
    k_bc = ValueBoundaryCondition(
        m -> logarithmic_tke_value(
            evaluate_wall_input(u_tau, m),
            parameters,
        ),
    )

    omega_bc = ValueBoundaryCondition(
        m -> logarithmic_omega_value(
            m.solution.k[1],
            Δf(m.grid, 1),
            parameters.z0_bottom,
            parameters,
        ),
    )

    return k_bc, omega_bc
end

"""
    TopLogarithmicWallConditions(u_tau; parameters)

Return `(k_bc, omega_bc)` value boundary conditions for a logarithmic
upper wall.
"""
function TopLogarithmicWallConditions(
    u_tau;
    parameters = KOmegaParameters(),
)
    k_bc = ValueBoundaryCondition(
        m -> logarithmic_tke_value(
            evaluate_wall_input(u_tau, m),
            parameters,
        ),
    )

    omega_bc = ValueBoundaryCondition(
        m -> logarithmic_omega_value(
            m.solution.k[m.grid.N],
            Δf(m.grid, m.grid.N),
            parameters.z0_top,
            parameters,
        ),
    )

    return k_bc, omega_bc
end

"""
    LogarithmicNeumannBoundaryConditions([FT=Float64]; parameters)

Return `(k_bcs, omega_bcs)` matching GOTM's logarithmic Neumann
boundary conditions. GOTM prescribes zero boundary flux for `k` and a
positive inward flux `k / (sigma_omega * (Delta z / 2 + z0))` for
`omega` at both walls. OceanTurb's top-flux sign convention is
outward-positive, so the top condition carries a minus sign.
"""
function LogarithmicNeumannBoundaryConditions(
    FT=Float64;
    parameters=KOmegaParameters(),
)
    k_bcs = ZeroFluxBoundaryConditions(FT)
    omega_bcs = OceanTurb.BoundaryConditions(
        FT;
        bottom=FluxBoundaryCondition(
            m -> logarithmic_omega_flux(
                m.solution.k[1],
                Δf(m.grid, 1) / 2,
                parameters.z0_bottom,
                parameters,
            ),
        ),
        top=FluxBoundaryCondition(
            m -> -logarithmic_omega_flux(
                m.solution.k[m.grid.N],
                Δf(m.grid, m.grid.N) / 2,
                parameters.z0_top,
                parameters,
            ),
        ),
    )

    return k_bcs, omega_bcs
end
