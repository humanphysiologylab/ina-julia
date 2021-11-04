using LabelledArrays
using Sundials
const a_syms = (
    :tau_m,
    :tau_h,
    :tau_j,
    :m_inf,
    :h_inf,
    :v_cp,
    :I_leak,
    :I_Na,
    :I_p,
    :I_c,
    :I_comp,
    :I_in,
)
const reltol = 1e-3
const abstol = LVector((
    v_comp = 1e-2,
    v_p = 1e-2,
    v_m = 1e-2,
    m = 1e-4,
    h = 1e-4,
    j = 1e-4,
    I_out = 1e-2,
))

const solver = CVODE_BDF();
const dt = 1e-9  # initial
const tspan_initial = (0.0, 10.0)
const saveat = 5e-5 # tspan[1]: 5e-5: tspan[2]
