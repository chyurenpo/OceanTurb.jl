@inline KU0(m) = m.background_diffusivities.momentum
@inline KV0(m) = m.background_diffusivities.momentum
@inline KT0(m) = m.background_diffusivities.tracer
@inline KS0(m) = m.background_diffusivities.tracer
@inline Kk0(m) = m.background_diffusivities.tke
@inline Komega0(m) = m.background_diffusivities.omega

@inline KU(m, i) = KU0(m) + onface(m.state.nu_t, i)
@inline KV(m, i) = KV0(m) + onface(m.state.nu_t, i)
@inline KT(m, i) = KT0(m) + onface(m.state.kappa_t, i)
@inline KS(m, i) = KS0(m) + onface(m.state.kappa_t, i)

@inline Kk(m, i) =
    Kk0(m) + onface(m.state.Kk, i)

@inline Komega(m, i) =
    Komega0(m) + onface(m.state.Komega, i)
