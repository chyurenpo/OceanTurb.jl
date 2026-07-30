#
# Tests for Convective Adjustment and Turbulent Kinetic Energy (CATKE).
#

const CA = OceanTurb.CATKE

@testset "CATKE time steppers" begin
    for stepper in (:ForwardEuler, :BackwardEuler)
        model = CA.Model(
            grid = UniformGrid(N=8, H=8.0),
            stepper = stepper,
        )
        time_step!(model, 1e-6, 2)
        CA.update_state!(model)

        @test iter(model) == 2
        @test all(isfinite, data(model.solution.e))
        @test all(isfinite, data(model.state.nu_t))
        @test all(isfinite, data(model.state.kappa_t))
    end
end

@testset "CATKE construction and diagnostic identities" begin
    grid = UniformGrid(N=16, H=20.0)
    model = CA.Model(
        grid = grid,
        stepper = :BackwardEuler,
        initial_tke = 4e-4,
    )

    @test length(model.solution) == 5
    @test model.grid === grid
    @test iter(model) == 0
    @test time(model) == 0

    set!(
        model.solution,
        U = z -> 0.01 * (z + grid.H),
        V = 0.0,
        T = z -> 10.0 + 0.01*z,
        S = 35.0,
        e = 4e-4,
    )
    CA.update_state!(model)

    for i = 2:grid.N
        w = CA.face_turbulent_velocity(model, i)
        @test model.state.nu_t[i] ≈ model.state.ell_u[i] * w
        @test model.state.kappa_t[i] ≈ model.state.ell_c[i] * w
        @test model.state.K_e[i] ≈ model.state.ell_e[i] * w
        @test model.state.Ri_face[i] ≈
              model.state.N2_face[i] / model.state.M2_face[i]
        @test CA.face_three_halves_tke(model, i) ≈
              (
                  model.solution.e[i-1]^(3/2) +
                  model.solution.e[i]^(3/2)
              ) / 2
    end

    for i in eachindex(model.solution.e)
        @test model.state.Ri[i] ≈
              oncell(model.state.Ri_face, i)
        @test model.state.epsilon[i] ≈
              model.solution.e[i]^(3/2) / model.state.ell_D[i]
        @test model.state.ell_D[i] > 0
    end
end

@testset "CATKE calibrated stability function" begin
    ml = CA.CATKEMixingLength()
    unstable = ml.C_tracer_unstable
    low = ml.C_tracer_low
    high = ml.C_tracer_high
    transition = ml.Ri_transition
    width = ml.Ri_width

    @test CA.stability_scale(
        -1.0, unstable, low, high, transition, width,
    ) == unstable
    @test CA.stability_scale(
        0.0, unstable, low, high, transition, width,
    ) == low
    @test CA.stability_scale(
        transition, unstable, low, high, transition, width,
    ) == low
    @test CA.stability_scale(
        transition + width/2,
        unstable,
        low,
        high,
        transition,
        width,
    ) ≈ (low + high) / 2
    @test CA.stability_scale(
        transition + width,
        unstable,
        low,
        high,
        transition,
        width,
    ) ≈ high

    # Regularizing the zero-shear limit avoids +Inf and -Inf producing
    # NaN when Richardson numbers are interpolated across an entrainment
    # interface.
    @test isfinite(CA.gradient_richardson_number(1e-5, 0.0))
    @test isfinite(CA.gradient_richardson_number(-1e-5, 0.0))
end

@testset "CATKE convective and penetrative mixing lengths" begin
    grid = UniformGrid(N=12, H=12.0)
    parameters = CA.CATKEParameters(
        minimum_tke = 1e-12,
        minimum_convective_buoyancy_flux = 1e-11,
    )
    model = CA.Model(
        grid = grid,
        parameters = parameters,
        initial_tke = 1e-5,
    )

    set!(
        model.solution,
        U = 0.0,
        V = 0.0,
        T = z -> 10.0 - 0.01*z,
        S = 35.0,
        e = 1e-5,
    )
    model.state.filtered_surface_buoyancy_flux = 1e-7
    CA.update_state!(model)

    i = 6
    w = CA.face_turbulent_velocity(model, i)
    expected_convective_length =
        model.mixing_length.C_tracer_convective * w^3 /
        (
            model.state.filtered_surface_buoyancy_flux +
            parameters.minimum_convective_buoyancy_flux
        )

    @test model.state.N2_face[i] < 0
    @test model.state.ell_c[i] ≈ expected_convective_length
    @test model.state.kappa_t[i] ≈
          expected_convective_length * w

    # Directly exercise the one-cell penetrative condition.
    model.state.N2_face[i] = 1e-4
    model.state.N2_face[i+1] = -1e-4
    model.state.M2_face[i] = 0.0
    expected_entrainment_length =
        model.mixing_length.C_tracer_entrainment *
        model.state.filtered_surface_buoyancy_flux /
        (
            w * model.state.N2_face[i] +
            parameters.minimum_convective_buoyancy_flux
        )

    @test CA.face_convective_length(
        model,
        i,
        model.mixing_length.C_tracer_convective,
        model.mixing_length.C_tracer_entrainment,
    ) ≈ expected_entrainment_length
end

@testset "CATKE surface forcing and filtered buoyancy flux" begin
    grid = UniformGrid(N=8, H=8.0)
    Jb = 1e-7
    momentum_flux_u = -1e-4
    constants = Constants()

    bcs = CA.ModelBoundaryConditions(
        U = BoundaryConditions(
            top = FluxBoundaryCondition(momentum_flux_u),
        ),
        T = BoundaryConditions(
            top = FluxBoundaryCondition(
                Jb / (constants.g * constants.α),
            ),
        ),
    )
    model = CA.Model(
        grid = grid,
        constants = constants,
        bcs = bcs,
    )

    u_tau = CA.friction_velocity(momentum_flux_u, 0.0)
    expected_tke_flux =
        -model.tke_equation.C_surface_shear * u_tau^3 -
        model.tke_equation.C_surface_convection *
        Δf(grid, grid.N) * Jb

    @test CA.surface_buoyancy_flux(model) ≈ Jb
    @test CA.surface_tke_flux(model) ≈ expected_tke_flux
    @test getbc(model, model.bcs.e.top) ≈ expected_tke_flux
    @test model.state.filtered_surface_buoyancy_flux == 0

    ell_D_top = CA.dissipation_length(model, grid.N)
    elapsed_time = 10.0
    mixing_time = cbrt(
        ell_D_top^2 /
        max(
            model.parameters.minimum_convective_buoyancy_flux,
            Jb,
        ),
    )
    weight = elapsed_time / mixing_time
    expected_filtered_flux = weight * Jb / (1 + weight)

    model.clock.time = elapsed_time
    CA.update_state!(model)

    @test model.state.filtered_surface_buoyancy_flux ≈
          expected_filtered_flux
    @test 0 < model.state.filtered_surface_buoyancy_flux < Jb
end

@testset "CATKE positivity and source-sink splitting" begin
    parameters = CA.CATKEParameters(minimum_tke=1e-12)
    tke_equation = CA.CATKEEquation(C_bottom_dissipation=0.0)
    model = CA.Model(
        grid = UniformGrid(N=4, H=4.0),
        parameters = parameters,
        tke_equation = tke_equation,
        stepper = :BackwardEuler,
    )

    set!(model.solution, e=-1.0)
    CA.update_state!(model)
    @test minimum(model.solution.e) ≥ parameters.minimum_tke

    i = 2
    model.solution.e[i] = 1.0
    model.state.P[i] = 0.1
    model.state.B[i] = -0.2
    model.state.ell_D[i] = 2.0

    @test CA.Re(model, i) ≈ 0.1
    @test CA.Le(model, i) ≈ 0.5 + 0.2

    model.state.B[i] = 0.2
    @test CA.Re(model, i) ≈ 0.3
    @test CA.Le(model, i) ≈ 0.5

    model.solution.e[i] = parameters.minimum_tke
    model.state.B[i] = -0.2
    model.state.ell_D[i] = 2.0
    @test CA.Le(model, i) ≈
          sqrt(parameters.minimum_tke) / 2

    time_step!(model, 1e-3, 2)
    CA.update_state!(model)
    @test all(isfinite, data(model.solution.e))
    @test minimum(model.solution.e) ≥ parameters.minimum_tke
end

@testset "CATKE backward-Euler homogeneous decay" begin
    grid = UniformGrid(N=6, H=6.0)
    background = CA.BackgroundDiffusivities(
        momentum = 0.0,
        tracer = 0.0,
        tke = 0.0,
    )
    tke_equation = CA.CATKEEquation(C_bottom_dissipation=0.0)
    parameters = CA.CATKEParameters(maximum_tke_diffusivity=0.0)
    model = CA.Model(
        grid = grid,
        parameters = parameters,
        background_diffusivities = background,
        tke_equation = tke_equation,
        stepper = :BackwardEuler,
        initial_tke = 1e-4,
    )

    CA.update_state!(model)
    e0 = copy(data(model.solution.e))
    decay_rate = [
        sqrt(e0[i]) / model.state.ell_D[i]
        for i = eachindex(e0)
    ]
    dt = 1e-3

    time_step!(model, dt, 1)
    expected = e0 ./ (1 .+ dt .* decay_rate)

    @test all(isapprox.(
        data(model.solution.e),
        expected;
        rtol = 5e-13,
    ))
end

@testset "CATKE mean-field and TKE forcing" begin
    forcing = CA.Forcing(
        U = (m, i) -> 1.0,
        V = (m, i) -> 2.0,
        T = (m, i) -> 3.0,
        S = (m, i) -> 4.0,
        e = (m, i) -> 5.0,
    )
    tke_equation = CA.CATKEEquation(C_bottom_dissipation=0.0)
    model = CA.Model(
        grid = UniformGrid(N=8, H=8.0),
        constants = Constants(f=0.0),
        forcing = forcing,
        tke_equation = tke_equation,
        stepper = :BackwardEuler,
    )

    i = 1
    @test CA.RU(model, i) == 1.0
    @test CA.RV(model, i) == 2.0
    @test CA.RT(model, i) == 3.0
    @test CA.RS(model, i) == 4.0

    # The TKE source also includes diagnosed production terms.
    @test CA.Re(model, i) ≥ 5.0
end

@testset "CATKE diffusivity caps and parameter validation" begin
    parameters = CA.CATKEParameters(
        maximum_viscosity = 1e-5,
        maximum_tracer_diffusivity = 2e-5,
        maximum_tke_diffusivity = 3e-5,
    )
    model = CA.Model(
        grid = UniformGrid(N=8, H=8.0),
        parameters = parameters,
        initial_tke = 1.0,
    )

    @test maximum(model.state.nu_t) ≤
          parameters.maximum_viscosity
    @test maximum(model.state.kappa_t) ≤
          parameters.maximum_tracer_diffusivity
    @test maximum(model.state.K_e) ≤
          parameters.maximum_tke_diffusivity

    ml = CA.CATKEMixingLength()
    te = CA.CATKEEquation()
    @test ml.C_surface == 1.131
    @test ml.C_tracer_convective == 4.793
    @test ml.C_tracer_entrainment == 0.112
    @test te.C_dissipation_convective == 3.254
    @test te.C_surface_shear == 3.179
    @test te.C_surface_convection == 0.383

    @test_throws AssertionError CA.Model(
        grid = UniformGrid(N=4, H=4.0),
        parameters = CA.CATKEParameters(minimum_tke=-1e-9),
    )
    @test_throws AssertionError CA.Model(
        grid = UniformGrid(N=4, H=4.0),
        mixing_length = CA.CATKEMixingLength(Ri_width=0.0),
    )
end
