@inline KU0(m) = m.background_diffusivities.momentum
@inline KV0(m) = m.background_diffusivities.momentum
@inline KT0(m) = m.background_diffusivities.tracer
@inline KS0(m) = m.background_diffusivities.tracer
@inline Ke0(m) = m.background_diffusivities.tke
@inline Kq2l0(m) = m.background_diffusivities.q2l

@inline KU(m, i) = KU0(m) + onface(m.state.νt, i)
@inline KV(m, i) = KV0(m) + onface(m.state.νt, i)
@inline KT(m, i) = KT0(m) + onface(m.state.κt, i)
@inline KS(m, i) = KS0(m) + onface(m.state.κt, i)

# MY transport coefficients q ℓ Sq and q ℓ Sl.
@inline Ke(m, i)   = Ke0(m)   + onface(m.state.Kk, i)
@inline Kq2l(m, i) = Kq2l0(m) + onface(m.state.Kq2l, i)
