using DelimitedFiles
using OceanTurb
using Printf

const CA = OceanTurb.CATKE
const KP = OceanTurb.KPP
const MY = OceanTurb.MellorYamada25
const KW = OceanTurb.KOmega
const KE = OceanTurb.KEpsilon

const experiment_dir = @__DIR__
const data_dir = joinpath(experiment_dir, "data")

const H = 256.0
const Nz = 128
const Δt = parse(
    Float64,
    get(ENV, "OCEANTURB_FREE_CONVECTION_DT", "10.0"),
)
const start_time = 10minute
const initial_tke = 1.0e-9

Δt > 0 ||
    error("OCEANTURB_FREE_CONVECTION_DT must be positive.")

const g = 9.81
const α = 2.0e-4
const Tref = 10.0
const N²_deep = 2.0e-6

const tracer_center = -96.0
const tracer_width = 8.0

const all_suites = (
    (hours=6,  Jb=9.6e-7, tracer_time=15minute, label="extreme forcing"),
    (hours=12, Jb=4.8e-7, tracer_time=30minute, label="strong forcing"),
    (hours=24, Jb=2.4e-7, tracer_time=1hour,    label="medium forcing"),
    (hours=48, Jb=1.2e-7, tracer_time=2hour,   label="weak forcing"),
    (hours=72, Jb=8.7e-8, tracer_time=4hour,   label="very weak forcing"),
)

function selected_suites()
    requested = strip(get(ENV, "OCEANTURB_SUITES", ""))
    isempty(requested) && return all_suites
    hours = Set(parse.(Int, split(requested, ",")))
    suites = filter(case -> case.hours in hours, all_suites)
    length(suites) == length(hours) ||
        error("OCEANTURB_SUITES contains an unsupported suite.")
    return suites
end

function timestep_token(dt)
    token = replace(@sprintf("%g", dt), "." => "p")
    return "dt$(token)s"
end

function read_les_case(hours)
    path = joinpath(
        data_dir,
        @sprintf("les_free_convection_%02dh.csv", hours),
    )
    isfile(path) ||
        error("Missing prepared LES profile $path. Run prepare_les_profiles.jl.")

    values, header = readdlm(path, ',', Float64, '\n'; header=true)
    size(values) == (Nz, 5) ||
        error("Unexpected dimensions $(size(values)) in $path.")

    return (
        path=path,
        z=collect(values[:, 1]),
        initial_b=collect(values[:, 2]),
        initial_c=collect(values[:, 3]),
        final_b=collect(values[:, 4]),
        final_c=collect(values[:, 5]),
    )
end

function make_boundary_conditions(namespace, temperature_flux)
    return namespace.ModelBoundaryConditions(
        U=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        V=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        T=BoundaryConditions(
            bottom=GradientBoundaryCondition(N²_deep / (g * α)),
            top=FluxBoundaryCondition(temperature_flux),
        ),
        S=BoundaryConditions(),
    )
end

function make_komega_boundary_conditions(parameters, temperature_flux)
    return KW.ModelBoundaryConditions(
        U=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        V=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        T=BoundaryConditions(
            bottom=GradientBoundaryCondition(N²_deep / (g * α)),
            top=FluxBoundaryCondition(temperature_flux),
        ),
        S=BoundaryConditions(),
        parameters=parameters,
    )
end

function make_kepsilon_boundary_conditions(parameters, temperature_flux)
    return KE.ModelBoundaryConditions(
        U=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        V=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        T=BoundaryConditions(
            bottom=GradientBoundaryCondition(N²_deep / (g * α)),
            top=FluxBoundaryCondition(temperature_flux),
        ),
        S=BoundaryConditions(),
        parameters=parameters,
    )
end

function build_models(case, les)
    grid = UniformGrid(N=Nz, H=H)
    maximum(abs.(collect(grid.zc) .- les.z)) < 1e-10 ||
        error("Prepared LES coordinates do not match the OceanTurb grid.")

    constants = Constants(α=α, β=0.0, g=g, f=1.0e-4)
    temperature_flux = case.Jb / (g * α)

    μplus = 1 / case.tracer_time
    μminus = μplus * sqrt(2π) * tracer_width / H
    passive_source = (model, i) ->
        μplus * exp(
            -(model.grid.zc[i] - tracer_center)^2 /
            (2 * tracer_width^2),
        ) - μminus

    catke_parameters = CA.CATKEParameters(minimum_tke=initial_tke)
    my25_parameters = MY.MY25Parameters()
    komega_parameters = KW.KOmegaParameters()
    kepsilon_parameters = KE.KEpsilonParameters()

    catke = CA.Model(
        grid=grid,
        constants=constants,
        parameters=catke_parameters,
        bcs=make_boundary_conditions(CA, temperature_flux),
        forcing=CA.Forcing(S=passive_source),
        stepper=:BackwardEuler,
        mixing_length=CA.CATKEMixingLength(C_bottom=1.0e12),
        tke_equation=CA.CATKEEquation(C_bottom_dissipation=0.0),
        initial_tke=initial_tke,
    )

    kpp = KP.Model(
        grid=grid,
        constants=constants,
        bcs=make_boundary_conditions(KP, temperature_flux),
        forcing=KP.Forcing(S=passive_source),
        stepper=:BackwardEuler,
    )

    # The LES archive has no TKE, q²ℓ, omega, or epsilon observations.
    # All prognostic closures therefore receive the paper's CATKE TKE floor.
    # The two-equation closures retain their canonical native second-variable
    # floors, matching a missing-field initialization rather than inventing
    # unobserved LES length scales.
    my25 = MY.Model(
        grid=grid,
        constants=constants,
        parameters=my25_parameters,
        bcs=make_boundary_conditions(MY, temperature_flux),
        forcing=MY.Forcing(S=passive_source),
        stepper=:BackwardEuler,
        initial_tke=initial_tke,
    )

    komega = KW.Model(
        grid=grid,
        constants=constants,
        parameters=komega_parameters,
        bcs=make_komega_boundary_conditions(
            komega_parameters,
            temperature_flux,
        ),
        forcing=KW.Forcing(S=passive_source),
        stepper=:BackwardEuler,
        initial_tke=initial_tke,
    )

    kepsilon = KE.Model(
        grid=grid,
        constants=constants,
        parameters=kepsilon_parameters,
        bcs=make_kepsilon_boundary_conditions(
            kepsilon_parameters,
            temperature_flux,
        ),
        forcing=KE.Forcing(S=passive_source),
        stepper=:BackwardEuler,
        initial_tke=initial_tke,
    )

    initial_temperature = Tref .+ les.initial_b ./ (g * α)
    models = (
        catke=catke,
        kpp=kpp,
        my25=my25,
        komega=komega,
        kepsilon=kepsilon,
    )

    for model in models
        OceanTurb.set!(
            model.solution,
            U=0.0,
            V=0.0,
            T=initial_temperature,
            S=les.initial_c,
        )
        model.clock.time = start_time
        model.clock.iter = 0
        OceanTurb.update!(model)

        initialized_b =
            g * α .* (collect(data(model.solution.T)) .- Tref)
        initialized_c = collect(data(model.solution.S))
        maximum(abs.(initialized_b .- les.initial_b)) < 1e-12 ||
            error("A closure did not receive the common LES buoyancy profile.")
        maximum(abs.(initialized_c .- les.initial_c)) < 1e-12 ||
            error("A closure did not receive the common LES tracer profile.")
    end

    return models
end

function validate_model(model, name, case_hours)
    for field in model.solution
        all(isfinite, data(field)) ||
            error("$name produced non-finite output in the $case_hours-hour case.")
    end
    isapprox(model.clock.time, case_hours * hour; atol=1e-8) ||
        error("$name did not reach the requested absolute final time.")
    return nothing
end

function initial_turbulence_diagnostics(models)
    return (
        catke_k=extrema(collect(data(models.catke.solution.e))),
        my25_k=extrema(collect(data(models.my25.solution.e))),
        my25_length=extrema(collect(data(models.my25.state.ℓ))),
        komega_k=extrema(collect(data(models.komega.solution.k))),
        komega_length=extrema(collect(data(models.komega.state.ell))),
        kepsilon_k=extrema(collect(data(models.kepsilon.solution.k))),
        kepsilon_length=extrema(collect(data(models.kepsilon.state.ell))),
    )
end

function closure_profiles(models)
    return (
        catke=(
            b=g * α .* (collect(data(models.catke.solution.T)) .- Tref),
            c=collect(data(models.catke.solution.S)),
            k=collect(data(models.catke.solution.e)),
        ),
        kpp=(
            b=g * α .* (collect(data(models.kpp.solution.T)) .- Tref),
            c=collect(data(models.kpp.solution.S)),
        ),
        my25=(
            b=g * α .* (collect(data(models.my25.solution.T)) .- Tref),
            c=collect(data(models.my25.solution.S)),
            k=collect(data(models.my25.solution.e)),
            second=collect(data(models.my25.solution.q2l)),
            length=collect(data(models.my25.state.ℓ)),
        ),
        komega=(
            b=g * α .* (collect(data(models.komega.solution.T)) .- Tref),
            c=collect(data(models.komega.solution.S)),
            k=collect(data(models.komega.solution.k)),
            second=collect(data(models.komega.solution.omega)),
            length=collect(data(models.komega.state.ell)),
        ),
        kepsilon=(
            b=g * α .* (collect(data(models.kepsilon.solution.T)) .- Tref),
            c=collect(data(models.kepsilon.solution.S)),
            k=collect(data(models.kepsilon.solution.k)),
            second=collect(data(models.kepsilon.solution.epsilon)),
            length=collect(data(models.kepsilon.state.ell)),
        ),
    )
end

function mixing_depths(models, z)
    return (
        catke=z[argmax(collect(data(models.catke.state.N2)))],
        kpp=-models.kpp.state.h,
        my25=z[argmax(collect(data(models.my25.state.N2)))],
        komega=z[argmax(collect(data(models.komega.state.N2)))],
        kepsilon=z[argmax(collect(data(models.kepsilon.state.N2)))],
    )
end

function write_profiles(path, z, les, profiles)
    open(path, "w") do io
        println(
            io,
            "z_m,initial_les_buoyancy_m_s-2,initial_les_passive_tracer," *
            "final_les_buoyancy_m_s-2,final_les_passive_tracer," *
            "catke_buoyancy_m_s-2,kpp_buoyancy_m_s-2," *
            "my25_buoyancy_m_s-2,komega_buoyancy_m_s-2," *
            "kepsilon_buoyancy_m_s-2,catke_passive_tracer," *
            "kpp_passive_tracer,my25_passive_tracer," *
            "komega_passive_tracer,kepsilon_passive_tracer," *
            "catke_tke_m2_s-2,my25_tke_m2_s-2,my25_q2l_m3_s-2," *
            "my25_length_m,komega_tke_m2_s-2,komega_omega_s-1," *
            "komega_length_m,kepsilon_tke_m2_s-2," *
            "kepsilon_epsilon_m2_s-3,kepsilon_length_m",
        )
        for i in eachindex(z)
            values = (
                z[i],
                les.initial_b[i],
                les.initial_c[i],
                les.final_b[i],
                les.final_c[i],
                profiles.catke.b[i],
                profiles.kpp.b[i],
                profiles.my25.b[i],
                profiles.komega.b[i],
                profiles.kepsilon.b[i],
                profiles.catke.c[i],
                profiles.kpp.c[i],
                profiles.my25.c[i],
                profiles.komega.c[i],
                profiles.kepsilon.c[i],
                profiles.catke.k[i],
                profiles.my25.k[i],
                profiles.my25.second[i],
                profiles.my25.length[i],
                profiles.komega.k[i],
                profiles.komega.second[i],
                profiles.komega.length[i],
                profiles.kepsilon.k[i],
                profiles.kepsilon.second[i],
                profiles.kepsilon.length[i],
            )
            println(io, join((@sprintf("%.10e", value) for value in values), ","))
        end
    end
end

const output_dir = joinpath(experiment_dir, "output", timestep_token(Δt))
mkpath(output_dir)
const summary_path = joinpath(output_dir, "run_summary.txt")

open(summary_path, "w") do summary
    println(summary, "OceanTurb free-convection LES example")
    println(summary, "Initialization: official LES b and c profiles at t=600 s")
    @printf(summary, "Grid: Nz=%d, H=%.1f m, dz=%.1f m\n", Nz, H, H / Nz)
    @printf(summary, "Common time step: %.6f s\n", Δt)
    @printf(summary, "Common prognostic-closure initial TKE: %.10e m2 s-2\n", initial_tke)
    println(summary, "CATKE bottom extensions: disabled")
    println(summary, "Two-equation second variables: canonical native floor initialization")
    println(summary)
    println(
        summary,
        "hours,Jb_m2_s-3,CATKE_base_z_m,KPP_base_z_m," *
        "MY25_base_z_m,komega_base_z_m,kepsilon_base_z_m",
    )

    @printf(
        "Running from the official LES t=600 s snapshot: Nz=%d, dz=%.1f m, dt=%g s\n",
        Nz,
        H / Nz,
        Δt,
    )

    for case in selected_suites()
        les = read_les_case(case.hours)
        models = build_models(case, les)
        initial_diagnostics = initial_turbulence_diagnostics(models)

        @printf(
            "\n%2d h case: Jb=%.2e, integrating %.3f h (600 s to %.0f s)\n",
            case.hours,
            case.Jb,
            (case.hours * hour - start_time) / hour,
            case.hours * hour,
        )
        @printf(
            "  initial length ranges [m]: MY2.5 %.3e–%.3e, k-omega %.3e–%.3e, k-epsilon %.3e–%.3e\n",
            initial_diagnostics.my25_length...,
            initial_diagnostics.komega_length...,
            initial_diagnostics.kepsilon_length...,
        )

        for (name, model) in pairs(models)
            elapsed = @elapsed run_until!(model, Δt, case.hours * hour)
            OceanTurb.update!(model)
            validate_model(model, string(name), case.hours)
            @printf("  %-8s completed in %.2f s\n", string(name), elapsed)
        end

        z = collect(nodes(models.catke.solution.T))
        profiles = closure_profiles(models)
        depths = mixing_depths(models, z)

        profile_path = joinpath(
            output_dir,
            @sprintf("free_convection_%02dh.csv", case.hours),
        )
        write_profiles(profile_path, z, les, profiles)

        @printf(
            summary,
            "%d,%.10e,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            case.hours,
            case.Jb,
            depths.catke,
            depths.kpp,
            depths.my25,
            depths.komega,
            depths.kepsilon,
        )
        flush(summary)

        @printf(
            "  mixed-layer base z [m]: CATKE %.1f, KPP %.1f, MY2.5 %.1f, k-omega %.1f, k-epsilon %.1f\n",
            depths.catke,
            depths.kpp,
            depths.my25,
            depths.komega,
            depths.kepsilon,
        )
        println("  wrote $profile_path")
    end
end

println("\nAll selected cases completed.")
println("Summary: $summary_path")
