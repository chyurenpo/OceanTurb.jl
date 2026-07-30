"""
    stability_scale(Ri, unstable, low, high, transition, width)

CATKE's piecewise-linear Richardson-number stability function.
"""
@inline function stability_scale(
    Ri,
    unstable,
    low,
    high,
    transition,
    width,
)
    ramp = clamp((Ri - transition) / width, zero(Ri), one(Ri))
    stable_scale = low + (high - low) * ramp
    return Ri < 0 ? unstable : stable_scale
end

@inline gradient_richardson_number(N2, M2) =
    N2 == 0 ? zero(N2) : N2 / max(M2, eps(one(M2)))

@inline function cell_turbulent_velocity(m, i)
    @inbounds e = max(m.solution.e[i], m.parameters.minimum_tke)
    return sqrt(e)
end

@inline function face_turbulent_velocity(m, i)
    return (
        cell_turbulent_velocity(m, i-1) +
        cell_turbulent_velocity(m, i)
    ) / 2
end

@inline function face_three_halves_tke(m, i)
    return (
        cell_turbulent_velocity(m, i-1)^3 +
        cell_turbulent_velocity(m, i)^3
    ) / 2
end

@inline function face_stable_base_length(m, i)
    p = m.parameters
    ml = m.mixing_length
    w = face_turbulent_velocity(m, i)
    @inbounds N2 = m.state.N2_face[i]

    surface_distance = ml.C_surface * (-m.grid.zf[i])
    bottom_distance = ml.C_bottom * (m.grid.zf[i] + m.grid.H)
    stratification_length = N2 > 0 ? w / sqrt(N2) : Inf

    return max(
        p.minimum_length,
        min(surface_distance, bottom_distance, stratification_length),
    )
end

@inline function cell_stable_base_length(m, i)
    p = m.parameters
    ml = m.mixing_length
    w = cell_turbulent_velocity(m, i)
    @inbounds N2 = m.state.N2[i]

    surface_distance = ml.C_surface * (-m.grid.zc[i])
    bottom_distance = ml.C_bottom * (m.grid.zc[i] + m.grid.H)
    stratification_length = N2 > 0 ? w / sqrt(N2) : Inf

    return max(
        p.minimum_length,
        min(surface_distance, bottom_distance, stratification_length),
    )
end

@inline function face_convective_length(
    m,
    i,
    C_convective,
    C_entrainment,
)
    p = m.parameters
    ml = m.mixing_length
    Jb = m.state.filtered_surface_buoyancy_flux
    Jb <= p.minimum_convective_buoyancy_flux &&
        return zero(Jb)

    denominator = Jb + p.minimum_convective_buoyancy_flux

    w = face_turbulent_velocity(m, i)
    w3 = face_three_halves_tke(m, i)
    @inbounds begin
        M2 = m.state.M2_face[i]
        N2 = m.state.N2_face[i]
        N2_above = i < m.grid.N ?
                   m.state.N2_face[i+1] :
                   N2
    end

    convective_length = C_convective * w3 / denominator
    depth = -m.grid.zf[i]
    flux_Ri = depth * w * M2 / denominator
    shear_factor = max(zero(flux_Ri), one(flux_Ri) - ml.C_sheared_plume * flux_Ri)
    convective_length *= shear_factor

    entrainment_length =
        C_entrainment * Jb /
        (w * max(N2, zero(N2)) + p.minimum_convective_buoyancy_flux)

    convecting = N2 < 0
    entraining =
        N2 > 0 &&
        N2_above < 0

    return convecting ? convective_length :
           entraining ? entrainment_length :
           zero(w)
end

@inline function cell_convective_length(
    m,
    i,
    C_convective,
    C_entrainment,
)
    p = m.parameters
    ml = m.mixing_length
    Jb = m.state.filtered_surface_buoyancy_flux
    Jb <= p.minimum_convective_buoyancy_flux &&
        return zero(Jb)

    denominator = Jb + p.minimum_convective_buoyancy_flux

    w = cell_turbulent_velocity(m, i)
    @inbounds begin
        M2 = m.state.M2[i]
        N2 = m.state.N2[i]
        N2_above = i < m.grid.N ? m.state.N2[i+1] : N2
    end

    convective_length = C_convective * w^3 / denominator
    depth = -m.grid.zc[i]
    flux_Ri = depth * w * M2 / denominator
    shear_factor = max(zero(flux_Ri), one(flux_Ri) - ml.C_sheared_plume * flux_Ri)
    convective_length *= shear_factor

    entrainment_length =
        C_entrainment * Jb /
        (w * max(N2, zero(N2)) + p.minimum_convective_buoyancy_flux)

    convecting = N2 < 0
    entraining =
        N2 > 0 &&
        N2_above < 0

    return convecting ? convective_length :
           entraining ? entrainment_length :
           zero(w)
end

@inline function face_total_length(
    m,
    i,
    C_unstable,
    C_low,
    C_high,
    C_convective,
    C_entrainment,
)
    ml = m.mixing_length
    @inbounds Ri = m.state.Ri_face[i]
    scale = stability_scale(
        Ri,
        C_unstable,
        C_low,
        C_high,
        ml.Ri_transition,
        ml.Ri_width,
    )

    stable_length = scale * face_stable_base_length(m, i)
    convective_length = face_convective_length(
        m,
        i,
        C_convective,
        C_entrainment,
    )

    return min(m.grid.H, max(stable_length, convective_length))
end

@inline function momentum_mixing_length(m, i)
    ml = m.mixing_length
    return face_total_length(
        m,
        i,
        ml.C_momentum_unstable,
        ml.C_momentum_low,
        ml.C_momentum_high,
        ml.C_momentum_convective,
        ml.C_momentum_entrainment,
    )
end

@inline function tracer_mixing_length(m, i)
    ml = m.mixing_length
    return face_total_length(
        m,
        i,
        ml.C_tracer_unstable,
        ml.C_tracer_low,
        ml.C_tracer_high,
        ml.C_tracer_convective,
        ml.C_tracer_entrainment,
    )
end

@inline function tke_mixing_length(m, i)
    ml = m.mixing_length
    return face_total_length(
        m,
        i,
        ml.C_tke_unstable,
        ml.C_tke_low,
        ml.C_tke_high,
        ml.C_tke_convective,
        ml.C_tke_entrainment,
    )
end

@inline function dissipation_length(m, i)
    ml = m.mixing_length
    te = m.tke_equation
    @inbounds Ri = m.state.Ri[i]
    scale = stability_scale(
        Ri,
        te.C_dissipation_unstable,
        te.C_dissipation_low,
        te.C_dissipation_high,
        ml.Ri_transition,
        ml.Ri_width,
    )

    stable_length = cell_stable_base_length(m, i) / scale
    convective_length = cell_convective_length(
        m,
        i,
        te.C_dissipation_convective,
        te.C_dissipation_entrainment,
    )

    return min(m.grid.H, max(stable_length, convective_length))
end
