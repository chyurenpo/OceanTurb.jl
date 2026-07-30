"""
    friction_velocity(momentum_flux_u, momentum_flux_v)

Return the ocean-side friction velocity associated with two kinematic
surface momentum fluxes.
"""
@inline friction_velocity(momentum_flux_u, momentum_flux_v) =
    (momentum_flux_u^2 + momentum_flux_v^2)^(1/4)

@inline friction_velocity(m::Model) =
    friction_velocity(
        m.state.surface_momentum_flux_u,
        m.state.surface_momentum_flux_v,
    )

"""
    surface_buoyancy_flux(model)

Return the surface buoyancy loss. Positive values destabilize the
water column and drive convection.
"""
@inline surface_buoyancy_flux(m::Model) =
    m.state.surface_buoyancy_flux

"""
    surface_tke_flux(model)

Return CATKE's calibrated upward TKE flux. It is negative when wind
stress or surface buoyancy loss injects TKE into the ocean.
"""
@inline function surface_tke_flux(m::Model)
    u_tau = friction_velocity(m)
    convective_velocity_cubed =
        Δf(m.grid, m.grid.N) *
        max(surface_buoyancy_flux(m), zero(eltype(m.grid)))

    return -m.tke_equation.C_surface_shear * u_tau^3 -
           m.tke_equation.C_surface_convection *
           convective_velocity_cubed
end

"""
    CATKETKEBoundaryConditions([FT=Float64])

Return CATKE's default no-flux bottom and calibrated surface-flux TKE
boundary conditions.
"""
function CATKETKEBoundaryConditions(FT=Float64)
    return OceanTurb.BoundaryConditions(
        FT;
        bottom = FluxBoundaryCondition(-zero(FT)),
        top = FluxBoundaryCondition(surface_tke_flux),
    )
end
