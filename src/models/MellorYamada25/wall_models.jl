"""
    friction_velocity(momentum_flux_u, momentum_flux_v)

Return `uτ = (Qu² + Qv²)^(1/4)` when `Qu` and `Qv` are kinematic
momentum fluxes [m² s⁻²].
"""
@inline friction_velocity(Qu, Qv) = (Qu^2 + Qv^2)^(1/4)

"""
GOTM logarithmic-wall Dirichlet value for `e=q²/2`.
"""
@inline function logarithmic_tke_value(
    uτ,
    p::MY25Parameters,
)
    T = promote_type(typeof(uτ), typeof(p.B1))
    exponent = T(Float32(2) / Float32(3))
    return T(0.5) * p.B1^exponent * uτ^2
end

"""
GOTM logarithmic-wall Dirichlet value for `q²ℓ`.
`distance` is the distance represented by the boundary stencil.
"""
@inline logarithmic_q2l_value(e, distance, z0, p::MY25Parameters) =
    2 * p.κ * e * (distance + z0)

"""
    logarithmic_q2l_flux(k, distance, z0, parameters)

Return GOTM's logarithmic Neumann boundary contribution for `q²l`.
GOTM defines this contribution as negative at both walls.
"""
@inline function logarithmic_q2l_flux(
    k,
    distance,
    z0,
    p::MY25Parameters,
)
    T = promote_type(
        typeof(k),
        typeof(distance),
        typeof(z0),
        typeof(p.Sl),
    )
    two_sqrt_two =
        T(Float32(2) * sqrt(Float32(2)))
    return -two_sqrt_two *
           p.Sl * p.κ^2 *
           max(k, p.k_min)^(T(3) / T(2)) *
           (distance + z0)
end

@inline evaluate_wall_input(value::Number, model) = value
@inline evaluate_wall_input(func, model) = func(model)

"""
    LogarithmicNeumannBoundaryConditions([FT=Float64]; parameters)

Return `(e_bcs, q2l_bcs)` matching GOTM's logarithmic Neumann
conditions. TKE has zero boundary flux. The `q²l` contribution is a
sink at both walls; OceanTurb's outward-positive top convention
therefore requires the opposite sign at the top.
"""
function LogarithmicNeumannBoundaryConditions(
    FT=Float64;
    parameters=MY25Parameters(),
)
    e_bcs = ZeroFluxBoundaryConditions(FT)
    q2l_bcs = OceanTurb.BoundaryConditions(
        FT;
        bottom=FluxBoundaryCondition(
            m -> logarithmic_q2l_flux(
                m.solution.e[1],
                Δf(m.grid, 1) / 2,
                parameters.z0_bottom,
                parameters,
            ),
        ),
        top=FluxBoundaryCondition(
            m -> -logarithmic_q2l_flux(
                m.solution.e[m.grid.N],
                Δf(m.grid, m.grid.N) / 2,
                parameters.z0_top,
                parameters,
            ),
        ),
    )

    return e_bcs, q2l_bcs
end

"""
    BottomLogarithmicWallConditions(uτ; parameters)

Return `(e_bc, q2l_bc)` ValueBoundaryConditions for a logarithmic
bottom wall. `uτ` may be a number or a callable `uτ(model)`.

The q²ℓ value uses the first-cell TKE and one full cell thickness,
matching GOTM's Dirichlet-wall stencil convention.
"""
function BottomLogarithmicWallConditions(
    uτ;
    parameters = MY25Parameters(),
)
    e_bc = ValueBoundaryCondition(
        m -> logarithmic_tke_value(
            evaluate_wall_input(uτ, m),
            parameters,
        ),
    )

    q2l_bc = ValueBoundaryCondition(
        m -> logarithmic_q2l_value(
            m.solution.e[1],
            Δf(m.grid, 1),
            parameters.z0_bottom,
            parameters,
        ),
    )

    return e_bc, q2l_bc
end

"""
    TopLogarithmicWallConditions(uτ; parameters)

Return `(e_bc, q2l_bc)` ValueBoundaryConditions for a logarithmic
upper wall. For an ocean free surface, a TKE-flux condition may be
more appropriate; this helper is provided for controlled GOTM
wall-condition comparisons.
"""
function TopLogarithmicWallConditions(
    uτ;
    parameters = MY25Parameters(),
)
    e_bc = ValueBoundaryCondition(
        m -> logarithmic_tke_value(
            evaluate_wall_input(uτ, m),
            parameters,
        ),
    )

    q2l_bc = ValueBoundaryCondition(
        m -> logarithmic_q2l_value(
            m.solution.e[m.grid.N],
            Δf(m.grid, m.grid.N),
            parameters.z0_top,
            parameters,
        ),
    )

    return e_bc, q2l_bc
end
