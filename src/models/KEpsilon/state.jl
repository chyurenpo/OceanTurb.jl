mutable struct State{F}
    ell      :: F

    M2       :: F
    N2       :: F
    P        :: F
    B        :: F

    nu_t     :: F
    kappa_t  :: F
    c_mu     :: F
    c_mu_prime :: F
    Kk       :: F
    Kepsilon :: F
end

function State(grid)
    fields = ntuple(_ -> CellField(grid), 11)
    return State(fields...)
end

@inline function turbulent_stability_functions(
    p,
    M2,
    N2,
    k,
    epsilon,
)
    turbulent_timescale_squared = (k / epsilon)^2
    alpha_M = turbulent_timescale_squared * M2
    alpha_N = turbulent_timescale_squared * N2
    return canuto_a_weak_equilibrium(
        alpha_M,
        alpha_N,
        p.c_mu0,
    )
end

@inline function copy_neumann_ghosts!(field, N)
    @inbounds begin
        field[0] = field[1]
        field[N+1] = field[N]
    end
    return nothing
end

@inline function initial_length_scale(grid, i, p)
    z = grid.zc[i]
    distance_bottom = z + grid.H + p.z0_bottom
    distance_top = -z + p.z0_top
    nearest_boundary = min(distance_bottom, distance_top)
    return max(p.length_min, p.kappa * nearest_boundary)
end

@inline epsilon_from_length(k, ell, p) =
    p.c_mu0^3 * k^(3/2) / max(ell, p.length_min)

function initialize_turbulence!(
    solution,
    grid,
    p,
    initial_tke,
    initial_epsilon,
)
    k0 = max(initial_tke, p.k_min)

    for i in eachindex(solution.k)
        epsilon0 = if initial_epsilon === nothing
            epsilon_from_length(
                k0,
                initial_length_scale(grid, i, p),
                p,
            )
        else
            initial_epsilon
        end

        @inbounds begin
            solution.k[i] = k0
            solution.epsilon[i] = max(epsilon0, p.epsilon_min)
        end
    end

    copy_neumann_ghosts!(solution.k, grid.N)
    copy_neumann_ghosts!(solution.epsilon, grid.N)
    return nothing
end

@inline shear_squared_at_face(m, i) =
    ∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2

@inline buoyancy_frequency_squared_at_face(m, i) =
    m.constants.g * (
        m.constants.α * ∂z(m.solution.T, i) -
        m.constants.β * ∂z(m.solution.S, i)
    )

# Avoid using not-yet-filled solution ghost cells while diagnosing the
# closure state.
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

@inline function length_limited_epsilon(k, epsilon, N2, p)
    limited_epsilon = max(epsilon, p.epsilon_min)

    if p.length_limit && N2 > 0
        epsilon_for_maximum_length =
            p.c_mu0^3 * k * sqrt(N2) /
            (p.galperin * sqrt(2))
        limited_epsilon = max(
            limited_epsilon,
            epsilon_for_maximum_length,
        )
    end

    return limited_epsilon
end

function clip_prognostic_turbulence!(m)
    p = m.parameters

    for i in eachindex(m.solution.k)
        @inbounds begin
            m.solution.k[i] = max(m.solution.k[i], p.k_min)
            m.solution.epsilon[i] = max(
                m.solution.epsilon[i],
                p.epsilon_min,
            )
        end
    end

    return nothing
end

function update_state!(m)
    p = m.parameters
    clip_prognostic_turbulence!(m)

    for i in eachindex(m.solution.k)
        @inbounds begin
            k = m.solution.k[i]
            M2 = cell_shear_squared(m, i)
            N2 = cell_buoyancy_frequency_squared(m, i)

            epsilon = length_limited_epsilon(
                k,
                m.solution.epsilon[i],
                N2,
                p,
            )
            m.solution.epsilon[i] = epsilon

            sqrt_k = sqrt(k)
            ell = p.c_mu0^3 * k * sqrt_k / epsilon

            c_mu, c_mu_prime = turbulent_stability_functions(
                p,
                M2,
                N2,
                k,
                epsilon,
            )

            nu_t = c_mu * sqrt_k * ell
            kappa_t = c_mu_prime * sqrt_k * ell
            Kk = nu_t / p.sigma_k
            Kepsilon = nu_t / p.sigma_epsilon

            P = nu_t * M2
            B = -kappa_t * N2

            m.state.ell[i] = ell

            m.state.M2[i] = M2
            m.state.N2[i] = N2
            m.state.P[i]  = P
            m.state.B[i]  = B

            m.state.nu_t[i]     = nu_t
            m.state.kappa_t[i]  = kappa_t
            m.state.c_mu[i]     = c_mu
            m.state.c_mu_prime[i] = c_mu_prime
            m.state.Kk[i]       = Kk
            m.state.Kepsilon[i] = Kepsilon
        end
    end

    N = m.grid.N
    for field in (
        m.state.ell,
        m.state.M2,
        m.state.N2,
        m.state.P,
        m.state.B,
        m.state.nu_t,
        m.state.kappa_t,
        m.state.c_mu,
        m.state.c_mu_prime,
        m.state.Kk,
        m.state.Kepsilon,
    )
        copy_neumann_ghosts!(field, N)
    end

    return nothing
end
