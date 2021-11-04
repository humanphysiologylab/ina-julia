using DifferentialEquations

cb_step_v = PresetTimeCallback(protocol.t, change_step_v!, save_positions = (false, false))

t_steady_state = (0.0:0.25:5.0)
t_steady_state = (0.0:0.6:9.0)
cb_steady_state =
    PresetTimeCallback(t_steady_state, change_du!, save_positions = (false, false))
callbackset = CallbackSet(cb_step_v, cb_steady_state)

fitness_progress_history = Vector{Vector{Float64}}()
cb_bbo = oc -> push!(fitness_progress_history, [num_func_evals(oc), best_fitness(oc)])

fitness_progress_history_NM = []
cb_NM = x -> push!(fitness_progress_history_NM, calculate_loss_segments!(
                                                                            data_optim,
                                                                            data_true,
                                                                            x,
                                                                            prob,
                                                                            solve_kwargs_default,
                                                                            weight,
                                                                            flag_tau_v_half))
