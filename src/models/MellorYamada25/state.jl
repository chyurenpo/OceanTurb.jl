mutable struct State{F}
    ℓ       :: F
    ε       :: F
    Lz      :: F
    alpha_M :: F
    alpha_N :: F
    c_mu    :: F
    c_mu_p  :: F

    M2      :: F
    N2      :: F
    P       :: F
    B       :: F

    νt      :: F
    κt      :: F
    Kk      :: F
    Kq2l    :: F
end

function State(grid)
    fields = ntuple(_ -> CellField(grid), 15)
    return State(fields...)
end

@inline function copy_neumann_ghosts!(field, N)
    @inbounds begin
        field[0] = field[1]
        field[N+1] = field[N]
    end
    return nothing
end

function initialize_turbulence!(
    solution,
    grid,
    p,
    initial_tke,
    initial_epsilon,
)
    e0 = max(initial_tke, p.k_min)
    ε0 = max(initial_epsilon, p.eps_min)
    ℓ0 = cde(p) * e0^(3/2) / ε0
    q2l0 = max(2 * e0 * ℓ0, p.q2l_min)

    for i in eachindex(solution.e)
        @inbounds begin
            solution.e[i] = e0
            solution.q2l[i] = q2l0
        end
    end

    copy_neumann_ghosts!(solution.e, grid.N)
    copy_neumann_ghosts!(solution.q2l, grid.N)
    return nothing
end

@inline shear_squared_at_face(m, i) =
    ∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2

@inline buoyancy_frequency_squared_at_face(m, i) =
    m.constants.g * (
        m.constants.α * ∂z(m.solution.T, i) -
        m.constants.β * ∂z(m.solution.S, i)
    )

# Avoid depending on not-yet-filled solution ghost cells while updating
# the diagnostic closure state.
@inline function cell_value_from_interior_faces(f, m, i)
    N = m.grid.N
    if i == 1
        return f(m, 2)
    elseif i == N
        return f(m, N)
    else
        return (f(m, i) + f(m, i+1)) / 2
    end
end

@inline cell_shear_squared(m, i) =
    cell_value_from_interior_faces(shear_squared_at_face, m, i)

@inline cell_buoyancy_frequency_squared(m, i) =
    cell_value_from_interior_faces(buoyancy_frequency_squared_at_face, m, i)

@inline function diagnostic_wall_length(m, i)
    p = m.parameters
    z = m.grid.zc[i]

    db = z + m.grid.H
    ds = -z

    db0 = db + p.z0_bottom
    ds0 = ds + p.z0_top

    if p.wall_length === :parabolic
        return p.κ * ds0 * db0 / (ds0 + db0)
    elseif p.wall_length === :triangular
        return p.κ * min(ds0, db0)
    elseif p.wall_length === :infinite_depth
        return p.κ * ds0
    else
        throw(ArgumentError("wall_length must be :parabolic, :triangular, or :infinite_depth"))
    end
end

@inline function diagnose_length_and_epsilon(e, q2l, N2, p)
    ℓmin = minimum_length(p)
    e = max(e, p.k_min)
    ℓ = q2l / (2 * e)

    if p.length_limit && N2 > 0
        ℓcrit = p.galperin * sqrt(2 * e / N2)
        ℓ = min(ℓ, ℓcrit)
    end

    ℓ = max(ℓ, ℓmin)
    ε = cde(p) * e^(3/2) / ℓ

    # GOTM resets L consistently if the dissipation floor is active.
    if ε < p.eps_min
        ε = p.eps_min
        ℓ = cde(p) * e^(3/2) / ε
    end

    return ℓ, ε
end

function clip_prognostic_turbulence!(m)
    p = m.parameters
    ℓmin = minimum_length(p)

    for i in eachindex(m.solution.e)
        @inbounds begin
            m.solution.e[i] = max(m.solution.e[i], p.k_min)
            local_q2l_min = max(p.q2l_min, 2 * m.solution.e[i] * ℓmin)
            m.solution.q2l[i] = max(m.solution.q2l[i], local_q2l_min)
        end
    end

    return nothing
end

function update_state!(m)
    p = m.parameters

    clip_prognostic_turbulence!(m)

    for i in eachindex(m.solution.e)
        @inbounds begin
            e = m.solution.e[i]
            q2l = m.solution.q2l[i]

            M2 = cell_shear_squared(m, i)
            N2 = cell_buoyancy_frequency_squared(m, i)
            Lz = max(diagnostic_wall_length(m, i), minimum_length(p))
            ℓ, ε =
                diagnose_length_and_epsilon(e, q2l, N2, p)

            alpha_M_raw = (e / ε)^2 * M2
            alpha_N_raw = (e / ε)^2 * N2
            alpha_M, alpha_N, c_mu, c_mu_p =
                kantha_clayson_weak_equilibrium(
                    alpha_M_raw,
                    alpha_N_raw,
                    p.c_mu0,
                )

            sqrt_e = sqrt(e)
            q = sqrt(2 * e)

            νt = c_mu   * sqrt_e * ℓ
            κt = c_mu_p * sqrt_e * ℓ
            Kk = p.Sq * q * ℓ
            Kq2l = p.Sl * q * ℓ

            P = νt * M2
            B = -κt * N2

            m.state.ℓ[i]       = ℓ
            m.state.ε[i]       = ε
            m.state.Lz[i]      = Lz
            m.state.alpha_M[i] = alpha_M
            m.state.alpha_N[i] = alpha_N
            m.state.c_mu[i]    = c_mu
            m.state.c_mu_p[i]  = c_mu_p

            m.state.M2[i]      = M2
            m.state.N2[i]      = N2
            m.state.P[i]       = P
            m.state.B[i]       = B

            m.state.νt[i]      = νt
            m.state.κt[i]      = κt
            m.state.Kk[i]      = Kk
            m.state.Kq2l[i]    = Kq2l
        end
    end

    N = m.grid.N
    for field in (
        m.state.ℓ,
        m.state.ε,
        m.state.Lz,
        m.state.alpha_M,
        m.state.alpha_N,
        m.state.c_mu,
        m.state.c_mu_p,
        m.state.M2,
        m.state.N2,
        m.state.P,
        m.state.B,
        m.state.νt,
        m.state.κt,
        m.state.Kk,
        m.state.Kq2l,
    )
        copy_neumann_ghosts!(field, N)
    end

    return nothing
end
