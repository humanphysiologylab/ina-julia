function calculate_v_cp!(u, p, a)
    @unpack v_c, alpha = p
    @unpack v_comp = u
    v_cp = v_c + (v_c - v_comp) * (1 / (1 - alpha) - 1)
    @pack! a = v_cp
end
function calculate_I_leak!(u, p, a)
    @unpack g_leak = p
    @unpack v_m = u
    I_leak = g_leak * v_m
    @pack! a = I_leak
end

function calculate_I_c!(u, p, a)
    @unpack c_m = p
    I_c = 1e9 * c_m * calculate_d_v_m(u, p, a)
    @pack! a = I_c
end

function calculate_I_p!(u, p, a)
    @unpack c_p = p
    I_p = 1e9 * c_p * calculate_d_v_p(u, p, a)
    @pack! a = I_p
end

function calculate_I_comp!(u, p, a)
    @unpack x_c_comp, c_m = p
    I_comp = 1e9 * x_c_comp * c_m * calculate_d_v_comp(u, p, a)
    @pack! a = I_comp
end

function calculate_I_Na!(u, p, a)
    @unpack g_max, v_rev = p
    @unpack v_m, m, h, j = u
    I_Na = g_max * h * m^3 * j * (v_m - v_rev)
    @pack! a =I_Na
end
