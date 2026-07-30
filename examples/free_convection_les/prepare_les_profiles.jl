using DelimitedFiles
using JLD2
using Printf

const experiment_dir = @__DIR__
const les_root = joinpath(experiment_dir, "input", "les_1m")
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

"""
Linearly interpolate (or, for the first target point, extrapolate) a
horizontally averaged LES profile from its native stretched grid to the
uniform OceanTurb cell centers.

This is the one-dimensional equivalent of the interpolation performed by
Oceananigans/ParameterEstimocean's `regrid!` in the archived comparison.
"""
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
        "buoyancy_flux_m2_s-3,tracer_forcing_timescale_s," *
        "source_Nz,target_Nz,target_dz_m,max_initial_interpolation_residual",
    )

    for hours in suites
        source_path = joinpath(
            les_root,
            "$(hours)_hour_suite",
            "free_convection_with_tracer_instantaneous_statistics.jld2",
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

            initial_b_native =
                interior_profile(file, "b", initial_key)
            initial_c_native =
                interior_profile(file, "c", initial_key)
            final_b_native =
                interior_profile(file, "b", final_key)
            final_c_native =
                interior_profile(file, "c", final_key)

            initial_b =
                interpolate_profile(source_z, initial_b_native, target_z)
            initial_c =
                interpolate_profile(source_z, initial_c_native, target_z)
            final_b =
                interpolate_profile(source_z, final_b_native, target_z)
            final_c =
                interpolate_profile(source_z, final_c_native, target_z)

            # Re-interpolate the target profile to the native centers as a
            # compact diagnostic of interpolation loss over the common range.
            common = findall(
                z -> first(target_z) <= z <= last(target_z),
                source_z,
            )
            reconstructed_b = interpolate_profile(
                target_z,
                initial_b,
                source_z[common],
            )
            residual = maximum(
                abs.(
                    reconstructed_b .-
                    initial_b_native[common]
                ),
            )

            output_path = joinpath(
                output_dir,
                @sprintf("les_free_convection_%02dh.csv", hours),
            )
            open(output_path, "w") do io
                println(
                    io,
                    "z_m,initial_buoyancy_m_s-2,initial_passive_tracer," *
                    "final_les_buoyancy_m_s-2,final_les_passive_tracer",
                )
                for i in eachindex(target_z)
                    @printf(
                        io,
                        "%.10e,%.10e,%.10e,%.10e,%.10e\n",
                        target_z[i],
                        initial_b[i],
                        initial_c[i],
                        final_b[i],
                        final_c[i],
                    )
                end
            end

            @printf(
                manifest,
                "%d,%s,%.1f,%.1f,%.10e,%.1f,%d,%d,%.6f,%.10e\n",
                hours,
                relpath(source_path, experiment_dir),
                start_time,
                final_time,
                file["parameters/buoyancy_flux"],
                file["parameters/tracer_forcing_timescale"],
                Nz,
                target_cells,
                target_dz,
                residual,
            )

            @printf(
                "Prepared %2d h: t=600 s key=%s, final key=%s, max b interpolation residual=%.3e m s^-2\n",
                hours,
                initial_key,
                final_key,
                residual,
            )
        end
    end
end

println("LES profiles written to $output_dir")
