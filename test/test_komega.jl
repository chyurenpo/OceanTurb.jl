#
# Tests for the oceanic k-omega turbulence closure.
#

const KW = OceanTurb.KOmega

@testset "KOmega time steppers" begin
    for stepper in (:ForwardEuler, :BackwardEuler)
        model = KW.Model(
            grid = UniformGrid(N=8, H=8.0),
            stepper = stepper,
        )
        time_step!(model, 1e-6, 2)
        KW.update_state!(model)

        @test iter(model) == 2
        @test all(isfinite, data(model.solution.k))
        @test all(isfinite, data(model.solution.omega))
    end
end

@testset "KOmega construction and diagnostic identities" begin
    grid = UniformGrid(N=16, H=20.0)
    parameters = KW.KOmegaParameters(
        length_limit = false,
        epsilon_min = 1e-16,
    )
    model = KW.Model(
        grid = grid,
        parameters = parameters,
        stepper = :BackwardEuler,
    )

    @test length(model.solution) == 6
    @test model.grid === grid
    @test iter(model) == 0
    @test time(model) == 0

    set!(
        model.solution,
        U = z -> 0.01 * (z + grid.H),
        V = 0.0,
        T = z -> 10.0 + 0.01*z,
        S = 35.0,
        k = 1e-5,
        omega = 0.2,
    )
    KW.update_state!(model)

    for i in eachindex(model.solution.k)
        k = model.solution.k[i]
        omega = model.solution.omega[i]

        @test model.state.epsilon[i] ≈
              parameters.c_mu0^4 * k * omega
        @test model.state.ell[i] ≈
              sqrt(k) / (parameters.c_mu0 * omega)
        @test model.state.nu_t[i] ≈
              model.state.c_mu[i] *
              sqrt(k) *
              model.state.ell[i]
        @test model.state.kappa_t[i] ≈
              model.state.c_mu_prime[i] *
              sqrt(k) *
              model.state.ell[i]
        @test model.state.Kk[i] ≈
              model.state.nu_t[i] / parameters.sigma_k
        @test model.state.Komega[i] ≈
              model.state.nu_t[i] / parameters.sigma_omega
    end
end

@testset "KOmega positivity and implicit time stepping" begin
    grid = UniformGrid(N=16, H=20.0)
    model = KW.Model(grid=grid, stepper=:BackwardEuler)

    set!(
        model.solution,
        U = z -> 0.01 * (z + grid.H),
        V = 0.0,
        T = z -> 10.0 + 0.02*z,
        S = 35.0,
        k = -1.0,
        omega = -1.0,
    )

    KW.update_state!(model)
    @test minimum(model.solution.k) ≥ model.parameters.k_min
    @test minimum(model.solution.omega) ≥ model.parameters.omega_min
    @test minimum(model.state.epsilon) ≥ model.parameters.epsilon_min

    time_step!(model, 0.1, 5)
    KW.update_state!(model)

    @test iter(model) == 5
    @test time(model) ≈ 0.5
    @test all(isfinite, data(model.solution.k))
    @test all(isfinite, data(model.solution.omega))
    @test all(isfinite, data(model.state.nu_t))
    @test all(isfinite, data(model.state.kappa_t))
    @test minimum(model.solution.k) ≥ model.parameters.k_min
    @test minimum(model.solution.omega) ≥ model.parameters.omega_min
    @test minimum(model.state.nu_t) ≥ 0
    @test minimum(model.state.kappa_t) ≥ 0
end

@testset "KOmega homogeneous decay" begin
    grid = UniformGrid(N=4, H=4.0)
    parameters = KW.KOmegaParameters(
        length_limit = false,
        epsilon_min = 1e-20,
    )
    zero_turbulence_flux_bcs =
        KW.ModelBoundaryConditions(
            k=ZeroFluxBoundaryConditions(),
            omega=ZeroFluxBoundaryConditions(),
            parameters=parameters,
        )
    model = KW.Model(
        grid = grid,
        parameters = parameters,
        bcs = zero_turbulence_flux_bcs,
        stepper = :BackwardEuler,
        initial_tke = 1e-3,
        initial_omega = 0.1,
    )

    k0 = model.solution.k[1]
    omega0 = model.solution.omega[1]
    dt = 1e-3
    number_of_steps = 1000
    final_time = dt * number_of_steps

    time_step!(model, dt, number_of_steps)

    decay_factor = 1 + parameters.c_omega2 *
                       parameters.c_mu0^4 *
                       omega0 *
                       final_time
    exact_omega = omega0 / decay_factor
    exact_k = k0 * decay_factor^(-1 / parameters.c_omega2)

    @test all(isapprox.(
        data(model.solution.omega),
        exact_omega;
        rtol = 5e-12,
    ))
    @test all(isapprox.(
        data(model.solution.k),
        exact_k;
        rtol = 5e-6,
    ))
end

@testset "KOmega stable-stratification length limit" begin
    grid = UniformGrid(N=12, H=12.0)
    parameters = KW.KOmegaParameters(
        length_limit = true,
        epsilon_min = 1e-20,
        galperin = 0.27,
    )
    model = KW.Model(grid=grid, parameters=parameters)

    set!(
        model.solution,
        T = z -> 10.0 + 0.05*z,
        S = 35.0,
        k = 1e-4,
        omega = 1e-6,
    )
    KW.update_state!(model)

    for i in eachindex(model.solution.k)
        @test model.state.N2[i] > 0
        critical_length = parameters.galperin *
                          sqrt(2 * model.solution.k[i] / model.state.N2[i])
        @test model.state.ell[i] ≤ critical_length * (1 + 10eps())
    end
end

@testset "KOmega mean-field forcing" begin
    grid = UniformGrid(N=8, H=8.0)
    constants = Constants(f=0.0)

    default_model = KW.Model(grid=grid, constants=constants)
    @test KW.RU(default_model, 1) == 0.0
    @test KW.RV(default_model, 1) == 0.0
    @test KW.RT(default_model, 1) == 0.0
    @test KW.RS(default_model, 1) == 0.0

    forcing = KW.Forcing(
        U = (m, i) -> 1.0,
        V = (m, i) -> 2.0,
        T = (m, i) -> 3.0,
        S = (m, i) -> 4.0,
    )
    model = KW.Model(
        grid = grid,
        constants = constants,
        forcing = forcing,
        stepper = :BackwardEuler,
    )

    @test KW.RU(model, 1) == 1.0
    @test KW.RV(model, 1) == 2.0
    @test KW.RT(model, 1) == 3.0
    @test KW.RS(model, 1) == 4.0

    set!(model.solution, U=0.0, V=0.0, T=0.0, S=0.0)
    dt = 0.01
    time_step!(model, dt, 1)

    @test all(isapprox.(data(model.solution.U), dt; atol=1e-12))
    @test all(isapprox.(data(model.solution.V), 2dt; atol=1e-12))
    @test all(isapprox.(data(model.solution.T), 3dt; atol=1e-12))
    @test all(isapprox.(data(model.solution.S), 4dt; atol=1e-12))
end

@testset "KOmega logarithmic wall conditions" begin
    grid = UniformGrid(N=8, H=8.0)
    parameters = KW.KOmegaParameters(
        z0_bottom = 1e-3,
        z0_top = 2e-3,
    )
    model = KW.Model(grid=grid, parameters=parameters)
    set!(model.solution, k=2e-4)

    u_tau = 0.01
    k_bottom, omega_bottom =
        KW.BottomLogarithmicWallConditions(
            u_tau;
            parameters = parameters,
        )

    @test getbc(model, k_bottom) ≈
          KW.logarithmic_tke_value(u_tau, parameters)
    @test getbc(model, omega_bottom) ≈
          KW.logarithmic_omega_value(
              model.solution.k[1],
              Δf(grid, 1),
              parameters.z0_bottom,
              parameters,
          )

    k_top, omega_top =
        KW.TopLogarithmicWallConditions(
            m -> 2u_tau;
            parameters = parameters,
        )

    @test getbc(model, k_top) ≈
          KW.logarithmic_tke_value(2u_tau, parameters)
    @test getbc(model, omega_top) ≈
          KW.logarithmic_omega_value(
              model.solution.k[grid.N],
              Δf(grid, grid.N),
              parameters.z0_top,
              parameters,
          )

    wall_bcs = KW.ModelBoundaryConditions(
        k = BoundaryConditions(
            bottom = k_bottom,
            top = FluxBoundaryCondition(0.0),
        ),
        omega = BoundaryConditions(
            bottom = omega_bottom,
            top = FluxBoundaryCondition(0.0),
        ),
    )
    wall_model = KW.Model(
        grid = grid,
        parameters = parameters,
        bcs = wall_bcs,
        stepper = :BackwardEuler,
    )
    time_step!(wall_model, 0.01, 1)

    @test all(isfinite, data(wall_model.solution.k))
    @test all(isfinite, data(wall_model.solution.omega))
end

@testset "KOmega source-sink splitting" begin
    parameters = KW.KOmegaParameters(
        c_omega3_stable = 1.0,
        length_limit = false,
        epsilon_min = 1e-20,
    )
    model = KW.Model(
        grid = UniformGrid(N=4, H=4.0),
        parameters = parameters,
    )

    i = 2
    model.solution.k[i] = 1.0
    model.solution.omega[i] = 3.0
    model.state.P[i] = 0.1
    model.state.B[i] = -1.0
    model.state.epsilon[i] = 0.25

    production = 3.0 * parameters.c_omega1 * 0.1
    buoyancy = 3.0 * parameters.c_omega3_stable * -1.0
    destruction = 3.0 * parameters.c_omega2 * 0.25

    @test production + buoyancy < 0
    @test KW.Romega(model, i) ≈ production
    @test KW.Lomega(model, i) ≈
          (destruction - buoyancy) / model.solution.omega[i]
    @test KW.Rk(model, i) ≈ model.state.P[i]
    @test KW.Lk(model, i) ≈
          (model.state.epsilon[i] - model.state.B[i]) /
          model.solution.k[i]
end

@testset "KOmega parameter validation" begin
    bad_parameters = KW.KOmegaParameters(c_mu0=-0.5477)
    @test_throws AssertionError KW.Model(
        grid = UniformGrid(N=4, H=4.0),
        parameters = bad_parameters,
    )
end

@testset "KOmega canonical GOTM configuration" begin
    coefficients = KW.canuto_a_komega_coefficients()

    # Reference values printed by the native GOTM k-omega runs. These
    # also test the default-real literal promotion in the Fortran port.
    @test coefficients.c_mu0 ≈
          0.52646471083532442 atol=2e-16
    @test coefficients.c_mu_shear_free ≈
          0.73100601349000971 atol=2e-16
    @test coefficients.kappa ≈
          0.39256055713193655 atol=2e-16
    @test coefficients.c_omega3_stable ≈
          -0.63861160645254555 atol=3e-16

    parameters = KW.KOmegaParameters()
    @test parameters.c_mu0 == coefficients.c_mu0
    @test parameters.c_omega3_stable ==
          coefficients.c_omega3_stable
    @test parameters.z0_bottom == 1.5e-3
    @test parameters.z0_top == 2e-2

    background = KW.BackgroundDiffusivities()
    @test background.momentum == 1.3e-6
    @test background.tracer == 1.4e-7
    @test background.tke == 0
    @test background.omega == 0

    grid = UniformGrid(N=8, H=8.0)
    model = KW.Model(grid=grid)

    expected_omega =
        parameters.epsilon_min /
        (parameters.c_mu0^4 * parameters.k_min)
    @test all(data(model.solution.k) .==
              parameters.k_min)
    @test all(data(model.solution.omega) .≈
              expected_omega)

    # With zero shear and stratification, alpha_M is clipped to
    # GOTM's small positive floor and c_mu approaches cmsf.
    @test model.state.c_mu[1] ≈
          coefficients.c_mu_shear_free rtol=1e-10
    @test model.state.c_mu_prime[1] > 0

    expected_bottom_flux =
        KW.logarithmic_omega_flux(
            model.solution.k[1],
            Δf(grid, 1) / 2,
            parameters.z0_bottom,
            parameters,
        )
    expected_top_flux =
        KW.logarithmic_omega_flux(
            model.solution.k[grid.N],
            Δf(grid, grid.N) / 2,
            parameters.z0_top,
            parameters,
        )

    @test getbc(model, model.bcs.k.bottom) == 0
    @test getbc(model, model.bcs.k.top) == 0
    @test getbc(model, model.bcs.omega.bottom) ≈
          expected_bottom_flux
    @test getbc(model, model.bcs.omega.top) ≈
          -expected_top_flux

    time_step!(model, 0.01, 2)
    @test all(isfinite, data(model.solution.k))
    @test all(isfinite, data(model.solution.omega))
    @test all(isfinite, data(model.state.nu_t))
    @test all(isfinite, data(model.state.kappa_t))

    @test_throws AssertionError KW.Model(
        grid=grid,
        initial_omega=1e-2,
        initial_epsilon=1e-12,
    )

    @test_throws MethodError KW.KOmegaParameters(
        stability_function=:constant,
    )
end
