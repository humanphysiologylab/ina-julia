using BlackBoxOptim
using BlackBoxOptim: num_func_evals, bboptimize
using Dates
using PyCall
using DataFrames
using Queryverse
using Logging
using DifferentialEquations
using Plots
#import Pkg; Pkg.add("ArgParse")
ScipyOptimize = pyimport("scipy.optimize")
@info "packages are loaded"

include("../../modules/INa.jl")
include("constants.jl")
include("io_utils.jl")
include("objectives.jl")
include("../../src/losses.jl")
@info "files are loaded"

parsed_args = parse_commandline()
config_name = abspath(parsed_args["config"])
@info "config $config_name"
config = make_config(config_name)
output_dir_name = make_output_dirs(config)
protocol = config["runtime"]["protocol"]

find_step(t, protocol = protocol) = protocol.v[findfirst(x -> x >= t, protocol.t)]

include("callbacks.jl")
@info "callbacks are loaded"

a = LVector(namedtuple(a_syms, zeros(length(a_syms))))
u₀ = config["runtime"]["states"]
u = deepcopy(u₀)
du = similar(u₀)
p_dict = config["runtime"]["constants"]
model_type = config["runtime"]["model_type"]
INa.compute_algebraic!(u₀, p_dict, a, model_type=model_type)

flag = config["runtime"]["flag"]
if flag == :act
    const tspan = (0.0, 5.0)
    const array_length = 100000
    const n_steps = 20
elseif flag == :inact
    const tspan = (0.0, 9.0)
    const array_length = 180000
    const n_steps = 15
else
    error("no such flag, choose only 'act' or 'inact'")
end

prob = create_ode_problem(du, u, a, config, tspan, cb_step_v)

solve_kwargs_default = (; reltol, abstol, solver, saveat, dt)
sol = solve_model(prob, (;), solve_kwargs_default);

p_dict = config["runtime"]["constants"]
p_keys_opt = keys(config["runtime"]["gene_keys_opt"])
mask_log = config["runtime"]["mask_multipliers"]
constants_bounds = config["runtime"]["constants_bounds"]
bounds = config["runtime"]["bounds"]
p_kwargs = (; p_keys_opt, p_dict, mask_log)

x₀ = [mask_log[k] ? log.(p_dict[k]) : p_dict[k] for k in p_keys_opt]

remade_bounds = [
    mask_log[k] ? Tuple{Float64,Float64}(log.(constants_bounds[k])) : #(bounds[k] * p_initial[k])) :
    Tuple{Float64,Float64}(constants_bounds[k]) for k in p_keys_opt
]

p = prepare_p(x₀, p_kwargs...)

data_true = config["initial"]["experimental_conditions"]["trace_1"]["phenotype"].I_out
flag_tau_v_half = config["initial"]["flag_tau_v_half"]
weight = config["initial"]["experimental_conditions"]["trace_1"]["sample_weight"]

popsize = config["initial"]["popsize"]#512
maxtime = config["initial"]["maxtime"]#3 * 60 * 60.
traceint = config["initial"]["traceint"]#60*15.
NM_maxiter = config["initial"]["NM_iters"]

fitness_progress_history = Vector{Vector{Float64}}()
callback_bbo = oc -> push!(fitness_progress_history, [num_func_evals(oc), best_fitness(oc)])

data_optim = zeros(array_length)
@info "BlackBoxOptim started"
res = bboptimize(
    x -> calculate_loss_segments!(
        data_optim,
        data_true,
        x,
        prob,
        solve_kwargs_default,
        weight,
        flag_tau_v_half,
        n_steps
        );
    SearchRange = remade_bounds,
    MaxTime = maxtime,
    TraceInterval = traceint,
    PopulationSize = popsize,
    # NThreads=Threads.nthreads()-1,
    CallbackFunction = callback_bbo,
    CallbackInterval = 0.5,
)

x_opt = best_candidate(res)
save_progress_history(fitness_progress_history, output_dir_name, "loss_BBO.csv")

res_all = DataFrame(name = [p_keys_opt...], real = x₀, BBO = x_opt)
res_all |> save(output_dir_name * "/res_BBO.csv", header = true)
@info "Nelder - Mead started"
res_n_m = ScipyOptimize.minimize(
    x -> calculate_loss_segments!(
        data_optim,
        data_true,
        x,
        prob,
        solve_kwargs_default,
        weight,
        flag_tau_v_half,
        n_steps
        ),
    x_opt,
    bounds = remade_bounds,
    callback = cb_NM,
    method = "Nelder-Mead",
    options = Dict(:maxiter => NM_maxiter),
)

x_opt_NM = res_n_m["x"]
df_fitness_progress_history_NM = DataFrame(numb = [k for k = 1:length(fitness_progress_history_NM)],
                                        loss = fitness_progress_history_NM)
df_fitness_progress_history_NM |> save(joinpath(output_dir_name, "loss_NM.csv"), header = true)

res_all = DataFrame(name = [p_keys_opt...], real = x₀, BBO = x_opt, NM = x_opt_NM)
res_all |> save(output_dir_name * "/res_NM.csv", header = true)