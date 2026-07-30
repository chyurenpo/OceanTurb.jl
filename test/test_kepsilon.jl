#
# Tests for the oceanic k-epsilon turbulence closure.
#

const KE = OceanTurb.KEpsilon

@testset "KEpsilon time steppers" begin
    for stepper in (:ForwardEuler, :BackwardEuler)
        model = KE.Model(
            grid = UniformGrid(N=8, H=8.0),
            stepper = stepper,
        )
        time_step!(model, 1e-6, 2)
        KE.update_state!(model)

        @test iter(model) == 2
        @test all(isfinite, data(model.solution.k))
        @test all(isfinite, data(model.solution.epsilon))
    end
end

@testset "KEpsilon construction and diagnostic identities" begin
    grid = UniformGrid(N=16, H=20.0)
    parameters = KE.KEpsilonParameters(
        length_limit = false,
        epsilon_min = 1e-16,
    )
    model = KE.Model(
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
        epsilon = 1e-6,
    )
    KE.update_state!(model)

    for i in eachindex(model.solution.k)
        k = model.solution.k[i]
        epsilon = model.solution.epsilon[i]

        @test model.state.ell[i] ≈
              parameters.c_mu0^3 * k^(3/2) / epsilon
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
        @test model.state.Kepsilon[i] ≈
              model.state.nu_t[i] / parameters.sigma_epsilon
    end
end

@testset "KEpsilon positivity and implicit time stepping" begin
    grid = UniformGrid(N=16, H=20.0)
    model = KE.Model(grid=grid, stepper=:BackwardEuler)

    set!(
        model.solution,
        U = z -> 0.01 * (z + grid.H),
        V = 0.0,
        T = z -> 10.0 + 0.02*z,
        S = 35.0,
        k = -1.0,
        epsilon = -1.0,
    )

    KE.update_state!(model)
    @test minimum(model.solution.k) ≥ model.parameters.k_min
    @test minimum(model.solution.epsilon) ≥
          model.parameters.epsilon_min

    time_step!(model, 0.1, 5)
    KE.update_state!(model)

    @test iter(model) == 5
    @test time(model) ≈ 0.5
    @test all(isfinite, data(model.solution.k))
    @test all(isfinite, data(model.solution.epsilon))
    @test all(isfinite, data(model.state.nu_t))
    @test all(isfinite, data(model.state.kappa_t))
    @test minimum(model.solution.k) ≥ model.parameters.k_min
    @test minimum(model.solution.epsilon) ≥
          model.parameters.epsilon_min
    @test minimum(model.state.nu_t) ≥ 0
    @test minimum(model.state.kappa_t) ≥ 0
end

@testset "KEpsilon homogeneous decay" begin
    grid = UniformGrid(N=4, H=4.0)
    parameters = KE.KEpsilonParameters(
        length_limit = false,
        epsilon_min = 1e-20,
    )
    background = KE.BackgroundDiffusivities(
        momentum = 0.0,
        tracer = 0.0,
        tke = 0.0,
        epsilon = 0.0,
    )
    zero_turbulence_flux_bcs =
        KE.ModelBoundaryConditions(
            k=ZeroFluxBoundaryConditions(),
            epsilon=ZeroFluxBoundaryConditions(),
            parameters=parameters,
        )
    model = KE.Model(
        grid = grid,
        parameters = parameters,
        background_diffusivities = background,
        bcs = zero_turbulence_flux_bcs,
        stepper = :BackwardEuler,
        initial_tke = 1e-3,
        initial_epsilon = 1e-4,
    )

    k0 = model.solution.k[1]
    epsilon0 = model.solution.epsilon[1]
    dt = 1e-3
    number_of_steps = 1000
    final_time = dt * number_of_steps

    time_step!(model, dt, number_of_steps)

    decay_factor = 1 +
                   (parameters.c_epsilon2 - 1) *
                   (epsilon0 / k0) *
                   final_time
    exact_k = k0 *
              decay_factor^(-1 / (parameters.c_epsilon2 - 1))
    exact_epsilon = epsilon0 *
                    decay_factor^(
                        -parameters.c_epsilon2 /
                        (parameters.c_epsilon2 - 1)
                    )

    @test all(isapprox.(
        data(model.solution.k),
        exact_k;
        rtol = 2e-4,
    ))
    @test all(isapprox.(
        data(model.solution.epsilon),
        exact_epsilon;
        rtol = 2e-4,
    ))
end

@testset "KEpsilon stable-stratification length limit" begin
    grid = UniformGrid(N=12, H=12.0)
    parameters = KE.KEpsilonParameters(
        length_limit = true,
        epsilon_min = 1e-20,
        galperin = 0.27,
    )
    model = KE.Model(grid=grid, parameters=parameters)

    set!(
        model.solution,
        T = z -> 10.0 + 0.05*z,
        S = 35.0,
        k = 1e-4,
        epsilon = 1e-12,
    )
    KE.update_state!(model)

    for i in eachindex(model.solution.k)
        @test model.state.N2[i] > 0
        critical_length = parameters.galperin *
                          sqrt(2 * model.solution.k[i] / model.state.N2[i])
        @test model.state.ell[i] ≤ critical_length * (1 + 10eps())
    end
end

@testset "KEpsilon mean-field forcing" begin
    grid = UniformGrid(N=8, H=8.0)
    constants = Constants(f=0.0)

    default_model = KE.Model(grid=grid, constants=constants)
    @test KE.RU(default_model, 1) == 0.0
    @test KE.RV(default_model, 1) == 0.0
    @test KE.RT(default_model, 1) == 0.0
    @test KE.RS(default_model, 1) == 0.0

    forcing = KE.Forcing(
        U = (m, i) -> 1.0,
        V = (m, i) -> 2.0,
        T = (m, i) -> 3.0,
        S = (m, i) -> 4.0,
    )
    model = KE.Model(
        grid = grid,
        constants = constants,
        forcing = forcing,
        stepper = :BackwardEuler,
    )

    @test KE.RU(model, 1) == 1.0
    @test KE.RV(model, 1) == 2.0
    @test KE.RT(model, 1) == 3.0
    @test KE.RS(model, 1) == 4.0

    set!(model.solution, U=0.0, V=0.0, T=0.0, S=0.0)
    dt = 0.01
    time_step!(model, dt, 1)

    @test all(isapprox.(data(model.solution.U), dt; atol=1e-12))
    @test all(isapprox.(data(model.solution.V), 2dt; atol=1e-12))
    @test all(isapprox.(data(model.solution.T), 3dt; atol=1e-12))
    @test all(isapprox.(data(model.solution.S), 4dt; atol=1e-12))
end

@testset "KEpsilon logarithmic wall conditions" begin
    grid = UniformGrid(N=8, H=8.0)
    parameters = KE.KEpsilonParameters(
        z0_bottom = 1e-3,
        z0_top = 2e-3,
    )
    model = KE.Model(grid=grid, parameters=parameters)
    set!(model.solution, k=2e-4)

    u_tau = 0.01
    k_bottom, epsilon_bottom =
        KE.BottomLogarithmicWallConditions(
            u_tau;
            parameters = parameters,
        )

    @test getbc(model, k_bottom) ≈
          KE.logarithmic_tke_value(u_tau, parameters)
    @test getbc(model, epsilon_bottom) ≈
          KE.logarithmic_epsilon_value(
              model.solution.k[1],
              Δf(grid, 1),
              parameters.z0_bottom,
              parameters,
          )

    k_top, epsilon_top =
        KE.TopLogarithmicWallConditions(
            m -> 2u_tau;
            parameters = parameters,
        )

    @test getbc(model, k_top) ≈
          KE.logarithmic_tke_value(2u_tau, parameters)
    @test getbc(model, epsilon_top) ≈
          KE.logarithmic_epsilon_value(
              model.solution.k[grid.N],
              Δf(grid, grid.N),
              parameters.z0_top,
              parameters,
          )

    wall_bcs = KE.ModelBoundaryConditions(
        k = BoundaryConditions(
            bottom = k_bottom,
            top = FluxBoundaryCondition(0.0),
        ),
        epsilon = BoundaryConditions(
            bottom = epsilon_bottom,
            top = FluxBoundaryCondition(0.0),
        ),
    )
    wall_model = KE.Model(
        grid = grid,
        parameters = parameters,
        bcs = wall_bcs,
        stepper = :BackwardEuler,
    )
    time_step!(wall_model, 0.01, 1)

    @test all(isfinite, data(wall_model.solution.k))
    @test all(isfinite, data(wall_model.solution.epsilon))
end

@testset "KEpsilon source-sink splitting" begin
    parameters = KE.KEpsilonParameters(
        c_epsilon3_stable = 1.0,
        length_limit = false,
        epsilon_min = 1e-20,
    )
    model = KE.Model(
        grid = UniformGrid(N=4, H=4.0),
        parameters = parameters,
    )

    i = 2
    model.solution.k[i] = 1.0
    model.solution.epsilon[i] = 3.0
    model.state.P[i] = 0.1
    model.state.B[i] = -1.0

    production = 3.0 * parameters.c_epsilon1 * 0.1
    buoyancy = 3.0 * parameters.c_epsilon3_stable * -1.0
    destruction = 3.0 * parameters.c_epsilon2 * 3.0

    @test production + buoyancy < 0
    @test KE.Repsilon(model, i) ≈ production
    @test KE.Lepsilon(model, i) ≈
          (destruction - buoyancy) / model.solution.epsilon[i]
    @test KE.Rk(model, i) ≈ model.state.P[i]
    @test KE.Lk(model, i) ≈
          (model.solution.epsilon[i] - model.state.B[i]) /
          model.solution.k[i]
end

@testset "KEpsilon parameter validation" begin
    bad_parameters = KE.KEpsilonParameters(c_mu0=-0.5477)
    @test_throws AssertionError KE.Model(
        grid = UniformGrid(N=4, H=4.0),
        parameters = bad_parameters,
    )

    bad_background = KE.BackgroundDiffusivities(epsilon=-1e-7)
    @test_throws AssertionError KE.Model(
        grid = UniformGrid(N=4, H=4.0),
        background_diffusivities = bad_background,
    )
end

@testset "KEpsilon canonical GOTM configuration" begin
    coefficients =
        KE.canuto_a_kepsilon_coefficients()

    # Reference values printed by the native GOTM k-epsilon runs.
    # These also test the default-real literal promotion in the
    # Fortran port.
    @test coefficients.c_mu0 ≈
          0.52646471083532442 atol=2e-16
    @test coefficients.c_mu_shear_free ≈
          0.73100601349000971 atol=2e-16
    @test coefficients.kappa ≈
          0.41587379967373655 atol=2e-16
    @test coefficients.c_epsilon3_stable ≈
          -0.62091212624900072 atol=3e-16

    parameters = KE.KEpsilonParameters()
    @test parameters.c_mu0 == coefficients.c_mu0
    @test parameters.c_epsilon3_stable ==
          coefficients.c_epsilon3_stable
    @test parameters.sigma_k == 1.0
    @test parameters.sigma_epsilon == 1.3
    @test parameters.c_epsilon1 == 1.44
    @test parameters.c_epsilon2 == 1.92
    @test parameters.c_epsilon3_unstable == 1.5
    @test parameters.z0_bottom == 1.5e-3
    @test parameters.z0_top == 2e-2

    background = KE.BackgroundDiffusivities()
    @test background.momentum == 1.3e-6
    @test background.tracer == 1.4e-7
    @test background.tke == 0
    @test background.epsilon == 0

    grid = UniformGrid(N=8, H=8.0)
    model = KE.Model(grid=grid)

    @test all(data(model.solution.k) .==
              parameters.k_min)
    @test all(data(model.solution.epsilon) .==
              parameters.epsilon_min)

    # With zero shear and stratification, alpha_M is clipped to
    # GOTM's small positive floor and c_mu approaches cmsf.
    @test model.state.c_mu[1] ≈
          coefficients.c_mu_shear_free rtol=1e-10
    @test model.state.c_mu_prime[1] > 0

    expected_bottom_flux =
        KE.logarithmic_epsilon_flux(
            model.solution.k[1],
            Δf(grid, 1) / 2,
            parameters.z0_bottom,
            parameters,
        )
    expected_top_flux =
        KE.logarithmic_epsilon_flux(
            model.solution.k[grid.N],
            Δf(grid, grid.N) / 2,
            parameters.z0_top,
            parameters,
        )

    @test getbc(model, model.bcs.k.bottom) == 0
    @test getbc(model, model.bcs.k.top) == 0
    @test getbc(model, model.bcs.epsilon.bottom) ≈
          expected_bottom_flux
    @test getbc(model, model.bcs.epsilon.top) ≈
          -expected_top_flux

    time_step!(model, 0.01, 2)
    @test all(isfinite, data(model.solution.k))
    @test all(isfinite, data(model.solution.epsilon))
    @test all(isfinite, data(model.state.nu_t))
    @test all(isfinite, data(model.state.kappa_t))

    @test_throws MethodError KE.KEpsilonParameters(
        turbulent_Prandtl=0.74,
    )
end
