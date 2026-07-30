using JLD2
using Printf

const experiment_dir = @__DIR__
const les_root = joinpath(
    experiment_dir,
    "..",
    "free_convection_les",
    "input",
    "les_1m",
)
const output_dir = joinpath(experiment_dir, "data")

const suites = (6, 12, 24, 48, 72)
const start_time = 10 * 60.0
const target_depth = 256.0
const target_cells = 128
const target_dz = target_depth / target_cells
const target_z = collect(
    range(
        -target_depth + target_dz / 2,
        stop=-target_dz / 2,
        length=target_cells,
    ),
)

function interpolate_profile(source_z, source_values, destination_z)
    issorted(source_z) || error("Source grid must be ordered bottom to top.")
    length(source_z) == length(source_values) ||
        error("Source coordinates and profile have different lengths.")

    destination_values = similar(destination_z)
    for i in eachindex(destination_z)
        z = destination_z[i]
        left = clamp(searchsortedlast(source_z, z), 1, length(source_z) - 1)
        z₁, z₂ = source_z[left], source_z[left + 1]
        v₁, v₂ = source_values[left], source_values[left + 1]
        destination_values[i] = v₁ + (z - z₁) * (v₂ - v₁) / (z₂ - z₁)
    end
    return destination_values
end

function time_key(file, requested_time)
    time_group = file["timeseries/t"]
    keys_without_serialized =
        filter(!=("serialized"), collect(keys(time_group)))
    pairs = [
        (Float64(file["timeseries/t/$key"]), key)
        for key in keys_without_serialized
    ]
    matches = filter(pair -> isapprox(pair[1], requested_time; atol=1e-8), pairs)
    length(matches) == 1 ||
        error("Expected exactly one snapshot at $requested_time seconds.")
    return only(matches)[2]
end

function interior_profile(file, field, key)
    Nz = Int(file["grid/Nz"])
    halo = Int(file["grid/Hz"])
    array = file["timeseries/$field/$key"]
    return vec(array[1, 1, halo + 1:halo + Nz])
end

mkpath(output_dir)

manifest_path = joinpath(output_dir, "les_profile_manifest.csv")
open(manifest_path, "w") do manifest
    println(
        manifest,
        "suite_hours,source_file,start_time_s,final_time_s," *
        "buoyancy_flux_m2_s-3,momentum_flux_m2_s-2," *
        "coriolis_parameter_s-1,tracer_forcing_timescale_s," *
        "source_Nz,target_Nz,target_dz_m",
    )

    for hours in suites
        source_path = joinpath(
            les_root,
            "$(hours)_hour_suite",
            "strong_wind_with_tracer_instantaneous_statistics.jld2",
        )
        isfile(source_path) || error("Missing official LES file: $source_path")

        jldopen(source_path, "r") do file
            Nz = Int(file["grid/Nz"])
            halo = Int(file["grid/Hz"])
            source_z_all = file["grid/zᵃᵃᶜ"]
            source_z = source_z_all[halo + 1:halo + Nz]

            initial_key = time_key(file, start_time)
            final_time = hours * 3600.0
            final_key = time_key(file, final_time)

            initial = NamedTuple{
                (:u, :v, :b, :c),
            }(
                Tuple(
                    interpolate_profile(
                        source_z,
                        interior_profile(file, field, initial_key),
                        target_z,
                    )
                    for field in ("u", "v", "b", "c")
                ),
            )
            final = NamedTuple{
                (:u, :v, :b, :c),
            }(
                Tuple(
                    interpolate_profile(
                        source_z,
                        interior_profile(file, field, final_key),
                        target_z,
                    )
                    for field in ("u", "v", "b", "c")
                ),
            )

            output_path = joinpath(
                output_dir,
                @sprintf("les_strong_wind_%02dh.csv", hours),
            )
            open(output_path, "w") do io
                println(
                    io,
                    "z_m,initial_u_m_s-1,initial_v_m_s-1," *
                    "initial_buoyancy_m_s-2,initial_passive_tracer," *
                    "final_les_u_m_s-1,final_les_v_m_s-1," *
                    "final_les_buoyancy_m_s-2,final_les_passive_tracer",
                )
                for i in eachindex(target_z)
                    @printf(
                        io,
                        "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                        target_z[i],
                        initial.u[i],
                        initial.v[i],
                        initial.b[i],
                        initial.c[i],
                        final.u[i],
                        final.v[i],
                        final.b[i],
                        final.c[i],
                    )
                end
            end

            @printf(
                manifest,
                "%d,%s,%.1f,%.1f,%.10e,%.10e,%.10e,%.1f,%d,%d,%.6f\n",
                hours,
                relpath(source_path, experiment_dir),
                start_time,
                final_time,
                file["parameters/buoyancy_flux"],
                file["parameters/momentum_flux"],
                file["parameters/coriolis_parameter"],
                file["parameters/tracer_forcing_timescale"],
                Nz,
                target_cells,
                target_dz,
            )

            @printf(
                "Prepared %2d h: tau_x=% .2e m2 s-2, t=600 s key=%s, final key=%s\n",
                hours,
                file["parameters/momentum_flux"],
                initial_key,
                final_key,
            )
        end
    end
end

println("LES profiles written to $output_dir")
