println(pwd())
# import Pkg
# using DataFrames
# using BlackBoxOptim
# using Queryverse
# using BlackBoxOptim: num_func_evals
# using LabelledArrays
# using NamedTupleTools
#
# include("csvs.jl")
# include("odes.jl")
# include("objectives.jl")
# include("../../src/losses.jl")
#
# p_initial = LVector((
#     g_max = 6392.4,
#     s_m = 12.8321,
#     b0_j = 11358.31,
#     tau_z = 0.0001,
#     k_m = 7.66677,
#     b0_h = 58131.407,
#     v_half_h = 59.69117,
#     tau_j_const = 0.00132,
#     c_m = 2.47e-11,
#     s_j = 69.417,
#     v_off = -2.0,
#     R_f = 326421.4,
#     c_p = 4.5e-13,
#     s_h = 21.62,
#     k_h = 5.507771,
#     x_c_comp = 0.05,
#     alpha = 0.75,
#     v_rev = 18.0,
#     v_half_m = 28.35802,
#     delta_h = 10.27583,
#     x_r_comp = 0.2,
#     v_c = -80.0,
#     g_leak = 0.4812443,
#     delta_m = 22.823,
#     a0_j = 0.6462,
#     a0_m = 16529.06,
#     delta_j = 7.445,
#     R = 2.495e7,
#     a0_h = 6.26,
#     b0_m = 386.7,
# ))
#
#
#
# u₀ = LVector((
#     v_comp = -80.0,
#     v_p = -80.0,
#     v_m = -80.0,
#     m = 0.0,
#     h = 1.0,
#     j = 1.0,
#     I_out = 0.0,
# ))
#
# a = @LVector Real a_syms
#
#
# calculate_model_ina! = calculate_model_1!
# compute_algebraic!(u₀, p_initial, a)
#
# du = similar(u₀)
# u = deepcopy(u₀)
#
# compute_rates!(du, u, p_initial, 0.0, a)
#
# rhs = ODEFunction((du, u, p, t) -> compute_rates!(du, u, p, t, a))
#
# cb_step_v = PresetTimeCallback(protocol.t, change_step_v!, save_positions = (false, false))
# t_steady_state = (0.0:0.25:5.0)
# cb_steady_state =
#     PresetTimeCallback(t_steady_state, change_du!, save_positions = (false, false))
# callbackset = CallbackSet(cb_step_v, cb_steady_state)
#
# u₁ = deepcopy(u₀)
# prob = ODEProblem(rhs, u₁, tspan, p_initial, callback = callbackset)
# solve_kwargs_default = (; reltol, abstol, solver, saveat, dt)
# sol = solve_model(prob, (;), solve_kwargs_default);
#
# u₁ = deepcopy(u₀)
# ###
# p_start = LVector(namedtuple(Symbol.(legend_constants_1.name), legend_constants_1.value));
# ###
# println(calculate_u∞!(u₁, p_start, prob, solve_kwargs_default))
#
# data_segmented = zeros(100001);
#
# prob_full = ODEProblem(rhs, u₀, tspan, p_start, callback = cb_step_v)
# @time trace = calculate_data_segments_serial(
#     data_segmented,
#     u₁,
#     p_start,
#     prob_full,
#     solve_kwargs_default,
# )
#
# data_segmented_u_0 = zeros(100001)
# @time calculate_data_segmented!(
#     data_segmented_u_0,
#     p_start,
#     prob_full,
#     solve_kwargs_default,
#     20,
# )
#
# loss_1 = sum((data_segmented .- data_segmented_u_0) .^ 2)
# loss_2 = calculate_rmse(data_segmented, data_segmented_u_0)
# data_segmented_u_1 = zeros(100001)
# @time loss_3 = calculate_loss_segments!(
#     data_segmented_u_1,
#     data_segmented_u_0,
#     p_start,
#     prob_full,
#     solve_kwargs_default,
# )
#
# println(loss_1, " ", loss_2, " ", loss_3)
