"""
General CATKE bounds and numerical parameters.
"""
Base.@kwdef struct CATKEParameters{T} <: AbstractParameters
    minimum_tke                      :: T = 1e-9
    minimum_convective_buoyancy_flux :: T = 1e-11
    minimum_length                   :: T = 1e-12

    maximum_viscosity          :: T = Inf
    maximum_tracer_diffusivity :: T = Inf
    maximum_tke_diffusivity    :: T = Inf
end

"""
Calibrated parameters for CATKE's diagnostic mixing lengths.

The defaults follow the maintained Oceananigans implementation of the
Wagner et al. (2025) closure. `high`, `low`, and `unstable` refer to
the corresponding gradient-Richardson-number regimes.
"""
Base.@kwdef struct CATKEMixingLength{T} <: AbstractParameters
    C_surface       :: T = 1.131
    C_bottom        :: T = 0.28
    C_sheared_plume :: T = 0.505

    Ri_width      :: T = 1.02
    Ri_transition :: T = 0.254

    C_momentum_high        :: T = 0.242
    C_momentum_low         :: T = 0.361
    C_momentum_unstable    :: T = 0.370
    C_momentum_convective  :: T = 3.705
    C_momentum_entrainment :: T = 0.0

    C_tracer_high        :: T = 0.098
    C_tracer_low         :: T = 0.369
    C_tracer_unstable    :: T = 0.572
    C_tracer_convective  :: T = 4.793
    C_tracer_entrainment :: T = 0.112

    C_tke_high        :: T = 0.548
    C_tke_low         :: T = 7.863
    C_tke_unstable    :: T = 1.447
    C_tke_convective  :: T = 3.642
    C_tke_entrainment :: T = 0.0
end

"""
Calibrated parameters for CATKE's TKE dissipation and surface flux.
"""
Base.@kwdef struct CATKEEquation{T} <: AbstractParameters
    C_dissipation_high        :: T = 0.579
    C_dissipation_low         :: T = 1.604
    C_dissipation_unstable    :: T = 0.923
    C_dissipation_convective  :: T = 3.254
    C_dissipation_entrainment :: T = 0.0

    C_surface_shear      :: T = 3.179
    C_surface_convection :: T = 0.383
    C_bottom_dissipation :: T = 1.0
end

"""
Optional additive background diffusivities in square metres per second.
CATKE's default minimum TKE already produces stratification-dependent
background mixing, so these default to zero.
"""
Base.@kwdef struct BackgroundDiffusivities{T} <: AbstractParameters
    momentum :: T = 0.0
    tracer   :: T = 0.0
    tke      :: T = 0.0
end

function validate_parameters(p)
    @assert p.minimum_tke > 0 "minimum_tke must be positive."
    @assert p.minimum_convective_buoyancy_flux > 0 "minimum_convective_buoyancy_flux must be positive."
    @assert p.minimum_length > 0 "minimum_length must be positive."
    @assert p.maximum_viscosity ≥ 0 "maximum_viscosity must be non-negative."
    @assert p.maximum_tracer_diffusivity ≥ 0 "maximum_tracer_diffusivity must be non-negative."
    @assert p.maximum_tke_diffusivity ≥ 0 "maximum_tke_diffusivity must be non-negative."
    return nothing
end

function validate_mixing_length(p)
    @assert p.C_surface > 0 "C_surface must be positive."
    @assert p.C_bottom > 0 "C_bottom must be positive."
    @assert p.C_sheared_plume ≥ 0 "C_sheared_plume must be non-negative."
    @assert p.Ri_width > 0 "Ri_width must be positive."

    for name in (
        :C_momentum_high,
        :C_momentum_low,
        :C_momentum_unstable,
        :C_momentum_convective,
        :C_momentum_entrainment,
        :C_tracer_high,
        :C_tracer_low,
        :C_tracer_unstable,
        :C_tracer_convective,
        :C_tracer_entrainment,
        :C_tke_high,
        :C_tke_low,
        :C_tke_unstable,
        :C_tke_convective,
        :C_tke_entrainment,
    )
        @assert getproperty(p, name) ≥ 0 "$name must be non-negative."
    end

    return nothing
end

function validate_tke_equation(p)
    @assert p.C_dissipation_high > 0 "C_dissipation_high must be positive."
    @assert p.C_dissipation_low > 0 "C_dissipation_low must be positive."
    @assert p.C_dissipation_unstable > 0 "C_dissipation_unstable must be positive."
    @assert p.C_dissipation_convective ≥ 0 "C_dissipation_convective must be non-negative."
    @assert p.C_dissipation_entrainment ≥ 0 "C_dissipation_entrainment must be non-negative."
    @assert p.C_surface_shear ≥ 0 "C_surface_shear must be non-negative."
    @assert p.C_surface_convection ≥ 0 "C_surface_convection must be non-negative."
    @assert p.C_bottom_dissipation ≥ 0 "C_bottom_dissipation must be non-negative."
    return nothing
end

function validate_background_diffusivities(K0)
    @assert K0.momentum ≥ 0 "background momentum diffusivity must be non-negative."
    @assert K0.tracer ≥ 0 "background tracer diffusivity must be non-negative."
    @assert K0.tke ≥ 0 "background TKE diffusivity must be non-negative."
    return nothing
end
