mutable struct State{F, C, T}
    surface_momentum_flux_u       :: T
    surface_momentum_flux_v       :: T
    surface_temperature_flux      :: T
    surface_salinity_flux         :: T
    surface_buoyancy_flux         :: T
    filtered_surface_buoyancy_flux :: T
    last_filter_time              :: T

    M2_face :: F
    N2_face :: F
    Ri_face :: F

    M2      :: C
    N2      :: C
    Ri      :: C
    P       :: C
    B       :: C
    ell_D   :: C
    epsilon :: C

    ell_u   :: F
    ell_c   :: F
    ell_e   :: F
    nu_t    :: F
    kappa_t :: F
    K_e     :: F
end

function State(grid)
    face_fields = ntuple(_ -> FaceField(grid), 9)
    cell_fields = ntuple(_ -> CellField(grid), 7)
    T = eltype(grid)

    return State(
        (zero(T) for _ = 1:7)...,
        face_fields[1],
        face_fields[2],
        face_fields[3],
        cell_fields...,
        face_fields[4],
        face_fields[5],
        face_fields[6],
        face_fields[7],
        face_fields[8],
        face_fields[9],
    )
end

@inline function copy_cell_ghosts!(field, N)
    @inbounds begin
        field[0] = field[1]
        field[N+1] = field[N]
    end
    return nothing
end

function initialize_tke!(solution, grid, p, initial_tke)
    e0 = max(initial_tke, p.minimum_tke)
    set!(solution.e, e0)
    copy_cell_ghosts!(solution.e, grid.N)
    return nothing
end

function clip_tke!(m)
    e_min = m.parameters.minimum_tke
    for i in eachindex(m.solution.e)
        @inbounds m.solution.e[i] = max(m.solution.e[i], e_min)
    end
    return nothing
end

@inline shear_squared_at_face(m, i) =
    ∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2

@inline buoyancy_frequency_squared_at_face(m, i) =
    m.constants.g * (
        m.constants.α * ∂z(m.solution.T, i) -
        m.constants.β * ∂z(m.solution.S, i)
    )

function update_surface_fluxes!(m)
    state = m.state
    state.surface_momentum_flux_u = getbc(m, m.bcs.U.top)
    state.surface_momentum_flux_v = getbc(m, m.bcs.V.top)
    state.surface_temperature_flux = getbc(m, m.bcs.T.top)
    state.surface_salinity_flux = getbc(m, m.bcs.S.top)
    state.surface_buoyancy_flux =
        m.constants.g * (
            m.constants.α * state.surface_temperature_flux -
            m.constants.β * state.surface_salinity_flux
        )
    return nothing
end

function update_gradient_diagnostics!(m)
    N = m.grid.N

    for i = 2:N
        @inbounds begin
            M2 = shear_squared_at_face(m, i)
            N2 = buoyancy_frequency_squared_at_face(m, i)
            m.state.M2_face[i] = M2
            m.state.N2_face[i] = N2
            m.state.Ri_face[i] = gradient_richardson_number(N2, M2)
        end
    end

    for field in (
        m.state.M2_face,
        m.state.N2_face,
        m.state.Ri_face,
    )
        @inbounds begin
            field[1] = field[2]
            field[N+1] = field[N]
        end
    end

    for i = 1:N
        @inbounds begin
            M2 = oncell(m.state.M2_face, i)
            N2 = oncell(m.state.N2_face, i)
            m.state.M2[i] = M2
            m.state.N2[i] = N2
            m.state.Ri[i] = oncell(m.state.Ri_face, i)
        end
    end

    return nothing
end

function update_filtered_surface_buoyancy_flux!(m)
    state = m.state
    current_time = m.clock.time
    elapsed_time = current_time - state.last_filter_time

    if elapsed_time > 0
        J_min = m.parameters.minimum_convective_buoyancy_flux
        J_old = state.filtered_surface_buoyancy_flux
        J_raw = state.surface_buoyancy_flux
        J_fast = max(J_min, J_old, J_raw)

        ell_D_top = max(
            m.parameters.minimum_length,
            dissipation_length(m, m.grid.N),
        )
        mixing_time = cbrt(ell_D_top^2 / J_fast)
        weight = elapsed_time / mixing_time

        state.filtered_surface_buoyancy_flux =
            (J_old + weight * J_raw) / (1 + weight)
        state.last_filter_time = current_time
    elseif elapsed_time < 0
        state.last_filter_time = current_time
    end

    return nothing
end

function update_length_scales_and_diffusivities!(m)
    p = m.parameters
    N = m.grid.N

    for i = 2:N
        @inbounds begin
            ell_u = momentum_mixing_length(m, i)
            ell_c = tracer_mixing_length(m, i)
            ell_e = tke_mixing_length(m, i)
            w = face_turbulent_velocity(m, i)

            m.state.ell_u[i] = ell_u
            m.state.ell_c[i] = ell_c
            m.state.ell_e[i] = ell_e

            m.state.nu_t[i] =
                min(p.maximum_viscosity, ell_u * w)
            m.state.kappa_t[i] =
                min(p.maximum_tracer_diffusivity, ell_c * w)
            m.state.K_e[i] =
                min(p.maximum_tke_diffusivity, ell_e * w)
        end
    end

    for field in (
        m.state.ell_u,
        m.state.ell_c,
        m.state.ell_e,
        m.state.nu_t,
        m.state.kappa_t,
        m.state.K_e,
    )
        @inbounds begin
            field[1] = field[2]
            field[N+1] = field[N]
        end
    end

    for i = 1:N
        @inbounds begin
            ell_D = dissipation_length(m, i)
            e = m.solution.e[i]
            m.state.ell_D[i] = ell_D
            m.state.epsilon[i] = e^(3/2) / ell_D

            lower_face_P = m.state.nu_t[i] * m.state.M2_face[i]
            upper_face_P = m.state.nu_t[i+1] * m.state.M2_face[i+1]
            lower_face_B = -m.state.kappa_t[i] * m.state.N2_face[i]
            upper_face_B = -m.state.kappa_t[i+1] * m.state.N2_face[i+1]

            m.state.P[i] = (lower_face_P + upper_face_P) / 2
            m.state.B[i] = (lower_face_B + upper_face_B) / 2
        end
    end

    for field in (
        m.state.M2,
        m.state.N2,
        m.state.Ri,
        m.state.P,
        m.state.B,
        m.state.ell_D,
        m.state.epsilon,
    )
        copy_cell_ghosts!(field, N)
    end

    return nothing
end

function update_state!(m)
    clip_tke!(m)
    update_surface_fluxes!(m)
    update_gradient_diagnostics!(m)
    update_filtered_surface_buoyancy_flux!(m)
    update_length_scales_and_diffusivities!(m)
    return nothing
end
