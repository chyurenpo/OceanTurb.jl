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
    get(ENV, "OCEANTURB_STRONG_WIND_DT", "10.0"),
)
const start_time = 10minute
const initial_tke = 1.0e-9

Δt > 0 ||
    error("OCEANTURB_STRONG_WIND_DT must be positive.")

const g = 9.81
const α = 2.0e-4
const Tref = 10.0
const f = 1.0e-4
const N²_deep = 2.0e-6

const tracer_center = -96.0
const tracer_width = 8.0

const all_suites = (
    (hours=6,  τx=-1.4e-3, tracer_time=15minute, label="extreme forcing"),
    (hours=12, τx=-9.0e-4, tracer_time=30minute, label="strong forcing"),
    (hours=24, τx=-6.8e-4, tracer_time=1hour,    label="medium forcing"),
    (hours=48, τx=-4.5e-4, tracer_time=2hour,   label="weak forcing"),
    (hours=72, τx=-4.1e-4, tracer_time=4hour,   label="very weak forcing"),
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
        @sprintf("les_strong_wind_%02dh.csv", hours),
    )
    isfile(path) ||
        error("Missing prepared LES profile $path. Run prepare_les_profiles.jl.")

    values, _ = readdlm(path, ',', Float64, '\n'; header=true)
    size(values) == (Nz, 9) ||
        error("Unexpected dimensions $(size(values)) in $path.")

    return (
        path=path,
        z=collect(values[:, 1]),
        initial_u=collect(values[:, 2]),
        initial_v=collect(values[:, 3]),
        initial_b=collect(values[:, 4]),
        initial_c=collect(values[:, 5]),
        final_u=collect(values[:, 6]),
        final_v=collect(values[:, 7]),
        final_b=collect(values[:, 8]),
        final_c=collect(values[:, 9]),
    )
end

function make_boundary_conditions(namespace, momentum_flux)
    return namespace.ModelBoundaryConditions(
        U=BoundaryConditions(top=FluxBoundaryCondition(momentum_flux)),
        V=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        T=BoundaryConditions(
            bottom=GradientBoundaryCondition(N²_deep / (g * α)),
            top=FluxBoundaryCondition(0.0),
        ),
        S=BoundaryConditions(),
    )
end

function make_komega_boundary_conditions(parameters, momentum_flux)
    return KW.ModelBoundaryConditions(
        U=BoundaryConditions(top=FluxBoundaryCondition(momentum_flux)),
        V=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        T=BoundaryConditions(
            bottom=GradientBoundaryCondition(N²_deep / (g * α)),
            top=FluxBoundaryCondition(0.0),
        ),
        S=BoundaryConditions(),
        parameters=parameters,
    )
end

function make_kepsilon_boundary_conditions(parameters, momentum_flux)
    return KE.ModelBoundaryConditions(
        U=BoundaryConditions(top=FluxBoundaryCondition(momentum_flux)),
        V=BoundaryConditions(top=FluxBoundaryCondition(0.0)),
        T=BoundaryConditions(
            bottom=GradientBoundaryCondition(N²_deep / (g * α)),
            top=FluxBoundaryCondition(0.0),
        ),
        S=BoundaryConditions(),
        parameters=parameters,
    )
end

function build_models(case, les)
    grid = UniformGrid(N=Nz, H=H)
    maximum(abs.(collect(grid.zc) .- les.z)) < 1e-10 ||
        error("Prepared LES coordinates do not match the OceanTurb grid.")

    constants = Constants(α=α, β=0.0, g=g, f=f)

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
        bcs=make_boundary_conditions(CA, case.τx),
        forcing=CA.Forcing(S=passive_source),
        stepper=:BackwardEuler,
        mixing_length=CA.CATKEMixingLength(C_bottom=1.0e12),
        tke_equation=CA.CATKEEquation(C_bottom_dissipation=0.0),
        initial_tke=initial_tke,
    )

    kpp = KP.Model(
        grid=grid,
        constants=constants,
        bcs=make_boundary_conditions(KP, case.τx),
        forcing=KP.Forcing(S=passive_source),
        stepper=:BackwardEuler,
    )

    my25 = MY.Model(
        grid=grid,
        constants=constants,
        parameters=my25_parameters,
        bcs=make_boundary_conditions(MY, case.τx),
        forcing=MY.Forcing(S=passive_source),
        stepper=:BackwardEuler,
        initial_tke=initial_tke,
    )

    komega = KW.Model(
        grid=grid,
        constants=constants,
        parameters=komega_parameters,
        bcs=make_komega_boundary_conditions(komega_parameters, case.τx),
        forcing=KW.Forcing(S=passive_source),
        stepper=:BackwardEuler,
        initial_tke=initial_tke,
    )

    kepsilon = KE.Model(
        grid=grid,
        constants=constants,
        parameters=kepsilon_parameters,
        bcs=make_kepsilon_boundary_conditions(kepsilon_parameters, case.τx),
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
            U=les.initial_u,
            V=les.initial_v,
            T=initial_temperature,
            S=les.initial_c,
        )
        model.clock.time = start_time
        model.clock.iter = 0
        OceanTurb.update!(model)

        initialized = (
            u=collect(data(model.solution.U)),
            v=collect(data(model.solution.V)),
            b=g * α .* (collect(data(model.solution.T)) .- Tref),
            c=collect(data(model.solution.S)),
        )
        maximum(abs.(initialized.u .- les.initial_u)) < 1e-12 ||
            error("A closure did not receive the common LES u profile.")
        maximum(abs.(initialized.v .- les.initial_v)) < 1e-12 ||
            error("A closure did not receive the common LES v profile.")
        maximum(abs.(initialized.b .- les.initial_b)) < 1e-12 ||
            error("A closure did not receive the common LES buoyancy profile.")
        maximum(abs.(initialized.c .- les.initial_c)) < 1e-12 ||
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

function closure_profiles(models)
    basic(model) = (
        u=collect(data(model.solution.U)),
        v=collect(data(model.solution.V)),
        b=g * α .* (collect(data(model.solution.T)) .- Tref),
        c=collect(data(model.solution.S)),
    )

    return (
        catke=merge(
            basic(models.catke),
            (k=collect(data(models.catke.solution.e)),),
        ),
        kpp=basic(models.kpp),
        my25=merge(
            basic(models.my25),
            (
                k=collect(data(models.my25.solution.e)),
                second=collect(data(models.my25.solution.q2l)),
                length=collect(data(models.my25.state.ℓ)),
                kappa=
                    models.my25.background_diffusivities.tracer .+
                    collect(data(models.my25.state.κt)),
            ),
        ),
        komega=merge(
            basic(models.komega),
            (
                k=collect(data(models.komega.solution.k)),
                second=collect(data(models.komega.solution.omega)),
                length=collect(data(models.komega.state.ell)),
                kappa=
                    models.komega.background_diffusivities.tracer .+
                    collect(data(models.komega.state.kappa_t)),
            ),
        ),
        kepsilon=merge(
            basic(models.kepsilon),
            (
                k=collect(data(models.kepsilon.solution.k)),
                second=collect(data(models.kepsilon.solution.epsilon)),
                length=collect(data(models.kepsilon.state.ell)),
                kappa=
                    models.kepsilon.background_diffusivities.tracer .+
                    collect(data(models.kepsilon.state.kappa_t)),
            ),
        ),
    )
end

function mixing_depths(models, z)
    function buoyancy_gradient_base(model)
        temperature = collect(data(model.solution.T))
        buoyancy = g * α .* (temperature .- Tref)
        N² = similar(buoyancy)
        N²[1] = (buoyancy[2] - buoyancy[1]) / (z[2] - z[1])
        for i in 2:length(z)-1
            N²[i] = (buoyancy[i+1] - buoyancy[i-1]) / (z[i+1] - z[i-1])
        end
        N²[end] =
            (buoyancy[end] - buoyancy[end-1]) / (z[end] - z[end-1])
        return z[argmax(N²)]
    end

    diagnosed = (
        catke=buoyancy_gradient_base(models.catke),
        kpp=buoyancy_gradient_base(models.kpp),
        my25=buoyancy_gradient_base(models.my25),
        komega=buoyancy_gradient_base(models.komega),
        kepsilon=buoyancy_gradient_base(models.kepsilon),
    )
    return merge(diagnosed, (kpp_internal=-models.kpp.state.h,))
end

function write_profiles(path, z, les, p)
    open(path, "w") do io
        println(
            io,
            "z_m,initial_les_u_m_s-1,initial_les_v_m_s-1," *
            "initial_les_buoyancy_m_s-2,initial_les_passive_tracer," *
            "final_les_u_m_s-1,final_les_v_m_s-1," *
            "final_les_buoyancy_m_s-2,final_les_passive_tracer," *
            "catke_u_m_s-1,kpp_u_m_s-1,my25_u_m_s-1," *
            "komega_u_m_s-1,kepsilon_u_m_s-1," *
            "catke_v_m_s-1,kpp_v_m_s-1,my25_v_m_s-1," *
            "komega_v_m_s-1,kepsilon_v_m_s-1," *
            "catke_buoyancy_m_s-2,kpp_buoyancy_m_s-2," *
            "my25_buoyancy_m_s-2,komega_buoyancy_m_s-2," *
            "kepsilon_buoyancy_m_s-2," *
            "catke_passive_tracer,kpp_passive_tracer," *
            "my25_passive_tracer,komega_passive_tracer," *
            "kepsilon_passive_tracer," *
            "catke_tke_m2_s-2,my25_tke_m2_s-2,my25_q2l_m3_s-2," *
            "my25_length_m,my25_tracer_diffusivity_m2_s-1," *
            "komega_tke_m2_s-2,komega_omega_s-1,komega_length_m," *
            "komega_tracer_diffusivity_m2_s-1," *
            "kepsilon_tke_m2_s-2,kepsilon_epsilon_m2_s-3," *
            "kepsilon_length_m,kepsilon_tracer_diffusivity_m2_s-1",
        )

        for i in eachindex(z)
            values = (
                z[i],
                les.initial_u[i], les.initial_v[i],
                les.initial_b[i], les.initial_c[i],
                les.final_u[i], les.final_v[i],
                les.final_b[i], les.final_c[i],
                p.catke.u[i], p.kpp.u[i], p.my25.u[i],
                p.komega.u[i], p.kepsilon.u[i],
                p.catke.v[i], p.kpp.v[i], p.my25.v[i],
                p.komega.v[i], p.kepsilon.v[i],
                p.catke.b[i], p.kpp.b[i], p.my25.b[i],
                p.komega.b[i], p.kepsilon.b[i],
                p.catke.c[i], p.kpp.c[i], p.my25.c[i],
                p.komega.c[i], p.kepsilon.c[i],
                p.catke.k[i],
                p.my25.k[i], p.my25.second[i],
                p.my25.length[i], p.my25.kappa[i],
                p.komega.k[i], p.komega.second[i],
                p.komega.length[i], p.komega.kappa[i],
                p.kepsilon.k[i], p.kepsilon.second[i],
                p.kepsilon.length[i], p.kepsilon.kappa[i],
            )
            println(io, join((@sprintf("%.10e", value) for value in values), ","))
        end
    end
end

const output_dir = joinpath(experiment_dir, "output", timestep_token(Δt))
mkpath(output_dir)
const summary_path = joinpath(output_dir, "run_summary.txt")

open(summary_path, "w") do summary
    println(summary, "OceanTurb rotating strong-wind LES example")
    println(summary, "Initialization: official LES u, v, b, and c profiles at t=600 s")
    @printf(summary, "Grid: Nz=%d, H=%.1f m, dz=%.1f m\n", Nz, H, H / Nz)
    @printf(summary, "Common time step: %.6f s\n", Δt)
    @printf(summary, "Coriolis parameter: %.10e s-1\n", f)
    println(summary, "Surface buoyancy flux: zero")
    @printf(summary, "Common prognostic-closure initial TKE: %.10e m2 s-2\n", initial_tke)
    println(summary, "CATKE bottom extensions: disabled")
    println(summary, "Two-equation second variables: canonical native floor initialization")
    println(summary, "OceanTurb models do not include an explicit Langmuir/Stokes-drift enhancement.")
    println(summary)
    println(
        summary,
        "hours,tau_x_m2_s-2,CATKE_base_z_m,KPP_base_z_m," *
        "MY25_base_z_m,komega_base_z_m,kepsilon_base_z_m,KPP_internal_h_base_z_m",
    )

    @printf(
        "Running rotating strong-wind example from official LES t=600 s: Nz=%d, dz=%.1f m, dt=%g s\n",
        Nz,
        H / Nz,
        Δt,
    )

    for case in selected_suites()
        les = read_les_case(case.hours)
        models = build_models(case, les)

        @printf(
            "\n%2d h case: tau_x=% .2e, integrating %.3f h (600 s to %.0f s)\n",
            case.hours,
            case.τx,
            (case.hours * hour - start_time) / hour,
            case.hours * hour,
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
            @sprintf("strong_wind_%02dh.csv", case.hours),
        )
        write_profiles(profile_path, z, les, profiles)

        @printf(
            summary,
            "%d,%.10e,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            case.hours,
            case.τx,
            depths.catke,
            depths.kpp,
            depths.my25,
            depths.komega,
            depths.kepsilon,
            depths.kpp_internal,
        )
        flush(summary)

        @printf(
            "  common max-N2 base z [m]: CATKE %.1f, KPP %.1f, MY2.5 %.1f, k-omega %.1f, k-epsilon %.1f\n",
            depths.catke,
            depths.kpp,
            depths.my25,
            depths.komega,
            depths.kepsilon,
        )
        println("  wrote $profile_path")
    end
end

println("\nAll selected strong-wind cases completed.")
println("Summary: $summary_path")
