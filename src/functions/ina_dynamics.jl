function calculate_tau_m_v_half(v_m, p)
    @unpack a0_m, s_m, v_half_m, k_m = p
    tau_m = 1 / (a0_m * exp(v_m / s_m) * (1 + exp(-(v_m+v_half_m) / k_m)))
end

function calculate_tau_m_a_b(v_m, p)
    @unpack a0_m, s_m, b0_m, delta_m = p
    tau_m = 1 / (a0_m * exp(v_m / s_m) + b0_m * exp(-v_m / delta_m))
end

function calculate_tau_h_v_half(v_m, p)
    @unpack a0_h, s_h, v_half_h, k_h = p
    tau_h =  1 / (a0_h * exp(-v_m / s_h) * (1 + exp((v_m+v_half_h) / k_h)))
end

function calculate_tau_h_a_b(v_m, p)
    @unpack a0_h, s_h, b0_h, delta_h = p
    tau_h = 1 / (a0_h * exp(-v_m / s_h) + b0_h * exp(v_m / delta_h))
end

function calculate_tau_j_v_half(v_m, p)
    @unpack tau_j_const, a0_j, s_j, v_half_h, k_h = p
    tau_j = tau_j_const +  1 / (a0_j * exp(-v_m / s_j) * (1 + exp((v_m+v_half_h) / k_h)))
end

function calculate_tau_j_a_b(v_m, p)
    @unpack tau_j_const, a0_j, s_j, b0_j, delta_j = p
    tau_j = tau_j_const + 1 / (a0_j * exp(-v_m / s_j) + b0_j * exp(v_m / delta_j))
end

function calculate_m_inf_v_half(v_m, p)
    @unpack v_half_m, k_m = p
    m_inf = 1 / (1 + exp(-(v_half_m + v_m) / k_m))
end

function calculate_m_inf_a_b(v_m, p)
    @unpack a0_m, s_m, b0_m, delta_m = p
    m_inf = 1 / (1 + (b0_m/a0_m)*exp(-v_m *(1/delta_m + 1/s_m)))
end

function calculate_h_inf_v_half(v_m, p)
    @unpack v_half_h, k_h = p
    h_inf = 1 / (1 + exp((v_half_h + v_m) / k_h))
end

function calculate_h_inf_a_b(v_m, p)
    @unpack a0_h, s_h, b0_h, delta_h = p
    m_inf = 1 / (1 + (b0_h/a0_h)*exp(v_m *(1/delta_h + 1/s_h)))
end
