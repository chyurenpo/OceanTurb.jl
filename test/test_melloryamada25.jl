#
# Regression tests for the GOTM-compatible Mellor-Yamada 2.5 closure.
#

const MY25 = OceanTurb.MellorYamada25

@testset "MellorYamada25 smoke test" begin
    grid = UniformGrid(N=16, H=20.0)
    model = MY25.Model(grid=grid, stepper=:BackwardEuler)

    set!(
        model.solution,
        U = z -> 0.01 * (z + grid.H),
        V = 0.0,
        T = z -> 10.0 + 0.02*z,
        S = 35.0,
        e = 1e-6,
    )

    for i in eachindex(model.solution.e)
        model.solution.q2l[i] = 2 * model.solution.e[i] * 0.1
    end

    time_step!(model, 0.1, 2)
    MY25.update_state!(model)

    @test all(isfinite, data(model.solution.e))
    @test all(isfinite, data(model.solution.q2l))
    @test minimum(model.solution.e) ≥ model.parameters.k_min
    @test minimum(model.solution.q2l) ≥ model.parameters.q2l_min
    @test minimum(model.state.νt) ≥ 0
    @test minimum(model.state.κt) ≥ 0
end

@testset "MellorYamada25 canonical GOTM configuration" begin
    coefficients = MY25.kantha_clayson_my25_coefficients()

    # Values printed by the native GOTM MY2.5 free-convection runs.
    @test coefficients.c_mu0 ≈
          0.55490553813915444 atol=2e-16
    @test coefficients.c_mu_shear_free ≈
          0.98842545407876936 atol=2e-16
    @test coefficients.B1 ≈
          16.553422823077717 atol=2e-15
    @test coefficients.kappa ≈
          0.40010965626841194 atol=2e-16
    @test coefficients.E3 ≈
          4.8258559267504708 atol=3e-15

    parameters = MY25.MY25Parameters()
    @test parameters.Ri_st == 0.23
    @test parameters.c_mu0 == coefficients.c_mu0
    @test parameters.c_mu_shear_free ==
          coefficients.c_mu_shear_free
    @test parameters.B1 == coefficients.B1
    @test parameters.E1 == 1.8
    @test parameters.E2 == 1.33
    @test parameters.E3 == coefficients.E3
    @test parameters.κ == coefficients.kappa
    @test parameters.galperin == 0.53
    @test parameters.z0_bottom == 1.5e-3
    @test parameters.z0_top == 2e-2

    second_moment =
        MY25.KanthaClaysonCoefficients()
    @test second_moment.ct2 == Float64(Float32(0.7))
    @test second_moment.ct3 == Float64(Float32(0.7))
    @test second_moment.ct5 == Float64(Float32(0.2))

    background = MY25.BackgroundDiffusivities()
    @test background.momentum == 1.3e-6
    @test background.tracer == 1.4e-7
    @test background.tke == 0
    @test background.q2l == 0

    grid = UniformGrid(N=8, H=8.0)
    model = MY25.Model(grid=grid)
    initial_length =
        parameters.c_mu0^3 *
        parameters.k_min^(3/2) /
        parameters.eps_min
    initial_q2l =
        2 * parameters.k_min * initial_length

    @test all(data(model.solution.e) .==
              parameters.k_min)
    @test all(data(model.solution.q2l) .≈
              initial_q2l)
    @test model.state.ℓ[1] ≈ initial_length
    @test model.state.ε[1] == parameters.eps_min

    # With zero shear and stratification, GOTM clips alpha_M to a
    # small positive value and approaches the shear-free coefficient.
    @test model.state.c_mu[1] ≈
          coefficients.c_mu_shear_free rtol=1e-10
    @test model.state.c_mu_p[1] > 0
end

@testset "MellorYamada25 GOTM stability functions" begin
    parameters = MY25.MY25Parameters()

    # Independent regression values from GOTM's cmue_c and cmue_d
    # algebra with the KC94 coefficient set.
    am, an, c_mu, c_mu_prime =
        MY25.kantha_clayson_weak_equilibrium(
            0.2,
            -0.1,
            parameters.c_mu0,
        )
    @test am == 0.2
    @test an == -0.1
    @test c_mu ≈ 1.0013647321828738 atol=2e-15
    @test c_mu_prime ≈
          1.075643181464473 atol=2e-15

    am, an, c_mu, c_mu_prime =
        MY25.kantha_clayson_weak_equilibrium(
            1.0,
            0.2,
            parameters.c_mu0,
        )
    @test am == 1.0
    @test an == 0.2
    @test c_mu ≈ 0.8755764829624187 atol=2e-15
    @test c_mu_prime ≈
          0.923969466601721 atol=2e-15

    c_mu_quasi, c_mu_prime_quasi =
        MY25.kantha_clayson_quasi_equilibrium(
            0.2,
            parameters.c_mu0,
        )
    @test c_mu_quasi ≈
          0.51377972232852 atol=2e-15
    @test c_mu_prime_quasi ≈
          0.6412072317143511 atol=2e-15

    # The canonical closure uses actual local shear rather than deriving
    # alpha_M diagnostically from alpha_N.
    _, _, c_mu_no_shear, _ =
        MY25.kantha_clayson_weak_equilibrium(
            0.0,
            0.2,
            parameters.c_mu0,
        )
    @test c_mu_no_shear != c_mu
end

@testset "MellorYamada25 logarithmic Neumann boundaries" begin
    parameters = MY25.MY25Parameters()
    grid = UniformGrid(N=8, H=8.0)
    model = MY25.Model(
        grid=grid,
        parameters=parameters,
    )

    expected_bottom_flux =
        MY25.logarithmic_q2l_flux(
            model.solution.e[1],
            Δf(grid, 1) / 2,
            parameters.z0_bottom,
            parameters,
        )
    expected_top_flux =
        MY25.logarithmic_q2l_flux(
            model.solution.e[grid.N],
            Δf(grid, grid.N) / 2,
            parameters.z0_top,
            parameters,
        )

    @test getbc(model, model.bcs.e.bottom) == 0
    @test getbc(model, model.bcs.e.top) == 0
    @test getbc(model, model.bcs.q2l.bottom) ≈
          expected_bottom_flux
    @test getbc(model, model.bcs.q2l.top) ≈
          -expected_top_flux
    @test expected_bottom_flux < 0
    @test expected_top_flux < 0
end

@testset "MellorYamada25 mean-field forcing" begin
    grid = UniformGrid(N=8, H=8.0)
    constants = Constants(f=0.0)
    forcing = MY25.Forcing(
        U = (m, i) -> 1.0,
        V = (m, i) -> 2.0,
        T = (m, i) -> 3.0,
        S = (m, i) -> 4.0,
    )
    model = MY25.Model(
        grid = grid,
        constants = constants,
        forcing = forcing,
        stepper = :BackwardEuler,
    )

    @test MY25.RU(model, 1) == 1.0
    @test MY25.RV(model, 1) == 2.0
    @test MY25.RT(model, 1) == 3.0
    @test MY25.RS(model, 1) == 4.0
end
