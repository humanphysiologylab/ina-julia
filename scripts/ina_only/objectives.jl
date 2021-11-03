using Parameters
using Statistics

find_step(t, protocol = protocol) = protocol.v[findfirst(x -> x >= t, protocol.t)]

function v_max_for_tau_a_b(a0, s, b0, delta)
    return log(a0 * delta / (b0 * s)) / (1 / s + 1 / delta)
end

function v_max_for_tau_v_half(k, s, v_half, sign = 1)
    return sign * k * log(k / (s - k)) - v_half
end

function sanity_checks(p, flag_tau_v_half = false)
    v_m_list = (-110, -20)
    if flag_tau_v_half
        if any((p.s_m, p.s_h, p.s_j) .<= (p.k_m, p.k_h, p.k_h))
            return [false]
        end

        v_max = (map(
            params -> v_max_for_tau_v_half(params...),
            [
                (p.k_m, p.s_m, p.v_half_m, -1),
                (p.k_h, p.s_h, p.v_half_h),
                (p.k_h, p.s_j, p.v_half_h),
            ],
        ))

        τ_m, τ_h, τ_j = [
            map(v_m -> calculate_τ(v_m, p), v_m_list) for calculate_τ ∈
            (INa.calculate_tau_m_v_half, INa.calculate_tau_h_v_half, INa.calculate_tau_j_v_half)
        ]

        τ_m_max, τ_h_max, τ_j_max = [
            INa.calculate_tau_m_v_half(v_max[1], p),
            INa.calculate_tau_h_v_half(v_max[2], p),
            INa.calculate_tau_j_v_half(v_max[3], p),
        ]
    else

        v_max = (map(
            params -> v_max_for_tau_a_b(params...),
            [
                (p.a0_m, -p.s_m, p.b0_m, -p.delta_m),
                (p.a0_h, p.s_h, p.b0_h, p.delta_h),
                (p.a0_j, p.s_j, p.b0_j, p.delta_j),
            ],
        ))

        τ_m, τ_h, τ_j = [
            map(v_m -> calculate_τ(v_m, p), v_m_list) for calculate_τ ∈
            (INa.calculate_tau_m_a_b, INa.calculate_tau_h_a_b,INa.calculate_tau_j_a_b)
        ]

        τ_m_max, τ_h_max, τ_j_max = [
            INa.calculate_tau_m_a_b(v_max[1], p),
            INa.calculate_tau_h_a_b(v_max[2], p),
            INa.calculate_tau_j_a_b(v_max[3], p),
        ]
    end

    flags = [
        all(1e-9 .<= τ_m .<= 1e-3),
        all(1e-8 .<= τ_h .<= 1e-2),
        all(5e-8 .<= τ_j .<= 2.0),
        all(v_m_list[1] .<= v_max .<= v_m_list[2]),
        τ_m_max >= 4e-5,
        τ_h_max >= 1e-3,
        τ_j_max >= 1e-2,
    ]
end

function prepare_p(x, p_keys_opt, p_dict, mask_log)#, kwargs_loss)
    @assert length(x) == length(p_keys_opt)
    p = deepcopy(p_dict)
    for (k, v, is_log) in zip(p_keys_opt, x, mask_log)
        p[k] = mask_log[k] ? exp(v) : v
    end
    #p["α"] = kwargs_loss.α
    return p
end

function solve_model(prob, remake_kwargs, solve_kwargs)
    prob_remade = remake(prob; remake_kwargs...)
    sol = solve(prob_remade; solve_kwargs...)
end

function calculate_u∞!(u∞, p, prob, solve_kwargs_default)::Symbol  # returns status
    @unpack reltol, abstol, solver, saveat, dt = solve_kwargs_default

    t_eq_end = 5.0
    remake_kwargs = (p = p, tspan = (0.0, t_eq_end), callback = nothing)
    solve_kwargs = (saveat = [t_eq_end], reltol, abstol, solver, dt)
    sol_eq = solve_model(prob, remake_kwargs, solve_kwargs)
    u∞[:] = sol_eq.u[1]
    return sol_eq.retcode
end


function calculate_data_segments_serial(
    data_segmented,
    u∞,
    p,
    prob,
    solve_kwargs_default,
    n_steps = 20,
)::Symbol  # returns status

    tspan_end = prob.tspan[2]
    step_size = solve_kwargs_default.saveat
    all_steps = 20
    segment_size_sec = tspan_end / all_steps
    segment_size_timestamps = Int(tspan_end / all_steps / step_size)

    for i ∈ 1:n_steps
        # Threads.@threads for i ∈ 1: n_steps
        tspan_segment = (segment_size_sec * (i - 1), segment_size_sec * i - step_size)
        remake_kwargs = (u0 = deepcopy(u∞), p = deepcopy(p), tspan = tspan_segment)
        sol_segment =
            solve_model(deepcopy(prob), remake_kwargs, deepcopy(solve_kwargs_default))

        if sol_segment.retcode == :Success
            i_start = (i - 1) * segment_size_timestamps + 1
            i_end = i * segment_size_timestamps
            data_segmented[i_start:i_end] = sol_segment[:I_out]
        else
            sol_segment.retcode
        end
    end

    return :Success
end


function calculate_data_segmented!(
    data_segmented,
    p,
    prob,
    solve_kwargs_default,
    n_steps = 20,
)::Symbol

    u∞ = deepcopy(prob.u0)
    retcode = calculate_u∞!(u∞, p, prob, solve_kwargs_default)

    if retcode != :Success
        return retcode
    end

    retcode = calculate_data_segments_serial(
        data_segmented,
        u∞,
        p,
        prob,
        solve_kwargs_default,
        n_steps,
    )

    if retcode != :Success
        return retcode
    end

    return :Success

end

function calculate_loss_segments!(
    data_segmented::Vector{Float64},
    data_true::Vector{Float64},
    x::Vector{Float64},
    prob::ODEProblem,
    solve_kwargs_default,
    weight = ones(length(data_segmented)),
    flag_tau_v_half = false,
    n_steps = 20,
)::Float64

    p = prepare_p(x, p_kwargs...)
    is_ok = all(sanity_checks(p, flag_tau_v_half))

    if !is_ok
        return Inf
    end

    retcode =
        calculate_data_segmented!(data_segmented, p, prob, solve_kwargs_default, n_steps)

    if retcode == :Success
        #residuals = data_segmented - data_true
        #@unpack c, α = p
        #loss = calculate_loss_robust.(residuals, α=α, c=c)
        loss = RMSE(data_segmented, data_true, weight)
        loss = mean(loss)
    else
        loss = Inf
    end
    return loss

end

function create_ode_problem(du, u, a, config, tspan, callback_name)
    u₀ = config["runtime"]["states"]
    protocol = config["runtime"]["protocol"]
    p_initial = config["runtime"]["constants"]
    model_type = config["runtime"]["model_type"]  #ex: A
    @info "model type $model_type"
    
    rhs = ODEFunction((du, u, p, t) -> INa.compute_rates!(du, u, p, t, a; model_type = model_type))
    prob = ODEProblem(rhs, u₀, tspan, p_initial, callback = callback_name)
    return prob
end

function change_step_v!(integrator)
    t = integrator.t
    v_c = find_step(t)
    integrator.p.v_c = v_c
    set_proposed_dt!(integrator, 1e-9)
end

function change_du!(integrator)
    integrator.u = u₀
end