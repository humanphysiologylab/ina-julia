####################
# MODELS
####################

function calculate_model_1!(u, p, a)
    @unpack v_m = u
    a.tau_m = calculate_tau_m_a_b(v_m, p)
    a.tau_h = calculate_tau_h_a_b(v_m, p)
    a.tau_j = calculate_tau_j_a_b(v_m, p)

    a.m_inf = calculate_m_inf_v_half(v_m, p)
    a.h_inf = calculate_h_inf_v_half(v_m, p)
end

function calculate_model_2!(u, p, a)
    @unpack v_m = u
    a.tau_m = calculate_tau_m_a_b(v_m, p)
    a.tau_h = calculate_tau_h_a_b(v_m, p)
    a.tau_j = calculate_tau_j_a_b(v_m, p)

    a.m_inf = calculate_m_inf_a_b(v_m, p)
    a.h_inf = calculate_h_inf_a_b(v_m, p)
end

function calculate_model_3!(u, p, a)
    @unpack v_m = u
    a.tau_m = calculate_tau_m_v_half(v_m, p)
    a.tau_h = calculate_tau_h_v_half(v_m, p)
    a.tau_j = calculate_tau_j_v_half(v_m, p)

    a.m_inf = calculate_m_inf_v_half(v_m, p)
    a.h_inf = calculate_h_inf_v_half(v_m, p)
end

####################
# ALGEBRAIC
####################

function compute_algebraic!(u, p, a; model_type=:A)

    if model_type == :A
        calculate_model_1!(u, p, a)
    elseif model_type == :B
        calculate_model_2!(u, p, a)
    elseif model_type == :C
        calculate_model_3!(u, p, a)
    else
        error("no such model")
    end

    calculate_v_cp!(u, p, a)
    calculate_I_leak!(u, p, a)
    calculate_I_Na!(u, p, a)
    calculate_I_c!(u, p, a)
    calculate_I_p!(u, p, a)
    calculate_I_comp!(u, p, a)

    a.I_in = a.I_leak + a.I_Na + a.I_c - a.I_comp + a.I_p
end

####################
# RATES
####################
function compute_rates!(du, u, p, t, a; model_type=:A, safety_factor = 1e-8)
    @unpack v_comp, v_p, v_m, m, h, j, I_out = u
    # @show a
    @unpack m_inf, tau_m, h_inf, tau_h, tau_j = a
    compute_algebraic!(u, p, a; model_type = model_type)

    du.v_comp = calculate_d_v_comp(u, p, a)
    du.v_p = calculate_d_v_p(u, p, a)
    du.v_m = calculate_d_v_m(u, p, a)
    du.m = calculate_d_gate(m_inf, m, tau_m + safety_factor)
    du.h = calculate_d_gate(h_inf, h, tau_h + safety_factor)
    du.j = calculate_d_gate(h_inf, j, tau_j)
    du.I_out = calculate_d_I_out(u, p, a)
end

