function calculate_d_v_comp(u, p, a)
    @unpack v_c, x_c_comp, c_m, x_r_comp, R, alpha = p
    @unpack v_comp = u
    d_v_comp = (v_c - v_comp) / (x_c_comp * c_m * x_r_comp * R * (1 - alpha))
end

function calculate_d_v_p(u, p, a)
    @unpack c_p, R_f = p
    @unpack v_cp = a
    @unpack v_p = u
    d_v_p = (v_cp - v_p) / (c_p * R_f)
end

function calculate_d_v_m(u, p, a)
    @unpack v_off, R, c_m = p
    @unpack v_m, v_p  = u
    @unpack I_leak, I_Na = a
    d_v_m = (v_p + v_off - v_m) / (R * c_m) - 1e-9 * (I_leak + I_Na) / c_m
end

function calculate_d_I_out(u, p, a)
    @unpack tau_z= p
    @unpack I_in = a
    @unpack I_out = u
    d_I_out = (I_in - I_out) / tau_z
end

function calculate_d_gate(gate_inf, gate, tau_gate)
    d_gate = (gate_inf - gate) / tau_gate
end
