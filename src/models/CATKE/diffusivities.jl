@inline KU(m, i) =
    m.background_diffusivities.momentum + m.state.nu_t[i]

@inline KV(m, i) =
    m.background_diffusivities.momentum + m.state.nu_t[i]

@inline KT(m, i) =
    m.background_diffusivities.tracer + m.state.kappa_t[i]

@inline KS(m, i) =
    m.background_diffusivities.tracer + m.state.kappa_t[i]

@inline Ke(m, i) =
    m.background_diffusivities.tke + m.state.K_e[i]
