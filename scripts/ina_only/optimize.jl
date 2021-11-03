using BlackBoxOptim
using BlackBoxOptim: num_func_evals, bboptimize
using Dates
using PyCall
using DataFrames
using Queryverse
using Logging
#import Pkg; Pkg.add("ArgParse")
using ArgParse
ScipyOptimize = pyimport("scipy.optimize")

include("io_utils.jl")
include("objectives.jl")
include("constants.jl")
include("../../src/losses.jl")

@info "packages are loaded"

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
            help = "a config name"
            required = true
    end
    return parse_args(s)
end


parsed_args = parse_commandline()

@debug "Parsed args:"
for (arg, val) in parsed_args
    @debug "  $arg  =>  $val"
end

config_name = abspath(parsed_args["config"])
config = make_config(config_name)

time = (Dates.format(Dates.now(), "yy_mm_dd_HHMMSS"))
output_dir = abspath(normpath(config["initial"]["output_folder_name"]))
output_dir_name = output_dir * time
if !isdir(output_dir)
    println("NO")
    mkdir(output_dir)
end
if !isdir(output_dir_name)
    println("NO")
    mkdir(output_dir_name)
end

calculate_model_ina! = config["runtime"]["model"]
protocol = config["runtime"]["protocol"]

a = @LVector Real a_syms
u₀ = config["runtime"]["states"]

du = similar(u₀)
u = deepcopy(u₀)

p_initial = config["runtime"]["constants"]

compute_algebraic!(u₀, p_initial, a)
compute_rates!(du, u, p_initial, 0., a)

find_step(t, protocol = protocol) = protocol.v[findfirst(x -> x >= t, protocol.t)]

cb_step_v = PresetTimeCallback(protocol.t, change_step_v!, save_positions = (false, false))

t_steady_state = (0.0:0.25:5.0)
cb_steady_state =
    PresetTimeCallback(t_steady_state, change_du!, save_positions = (false, false))
callbackset = CallbackSet(cb_step_v, cb_steady_state);

rhs = ODEFunction((du, u, p, t) -> compute_rates!(du, u, p, t, a))
prob = ODEProblem(rhs, u₀, tspan, p_initial,  callback = cb_step_v)

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

data_true = similar(zeros(100000))
calculate_data_segmented!(data_true, p, prob, solve_kwargs_default)
# data_true = config["initial"]["experimental_conditions"]["trace_1"]["phenotype"].I_out
flag_tau_v_half = config["initial"]["flag_tau_v_half"]
weight = ones(100000)#config["initial"]["experimental_conditions"]["trace_1"]["sample_weight"]
function loss_for_optim(x)
    data_optim = zeros(100000)
    return calculate_loss_segments!(
        data_optim,
        data_true,
        x,
        prob,
        solve_kwargs_default,
        weight,
        flag_tau_v_half,

    )
end

popsize = config["initial"]["popsize"]#512
maxtime = config["initial"]["maxtime"]#3 * 60 * 60.
traceint = config["initial"]["traceint"]#60*15.
NM_maxiter = config["initial"]["NM_iters"]

fitness_progress_history = Vector{Vector{Float64}}()
callback_bbo = oc -> push!(fitness_progress_history, [num_func_evals(oc), best_fitness(oc)])

res = bboptimize(
    loss_for_optim;
    SearchRange = remade_bounds,
    MaxTime = maxtime,
    TraceInterval = traceint,
    PopulationSize = popsize,
    #     NThreads=Threads.nthreads()-1,
    CallbackFunction = callback_bbo,
    CallbackInterval = 0.5,
)
x_opt = best_candidate(res)

df_fitness_progress_history = DataFrame(
    numb = [fitness_progress_history[k][1] for k = 1:length(fitness_progress_history)],
    loss = [fitness_progress_history[k][2] for k = 1:length(fitness_progress_history)],
)

df_fitness_progress_history |> save(output_dir_name * "/loss_BBO.csv", header = true)

res_all = DataFrame(name = [p_keys_opt...], real = x₀, BBO = x_opt)
res_all |> save(output_dir_name * "/res_BBO.csv", header = true)

function CallbackNM(x)
    push!(FitnessProgressHistoryNM, loss_for_optim(x))
end
FitnessProgressHistoryNM = []
res_n_m = ScipyOptimize.minimize(
    (x, args) -> loss_for_optim(x),
    x_opt,
    constants_bounds,
    callback = CallbackNM,
    method = "Nelder-Mead",
    options = Dict(:maxiter => NM_maxiter),
)
x_opt_NM = res_n_m["x"]
df_FitnessProgressHistoryNM = DataFrame(loss = FitnessProgressHistoryNM)
df_FitnessProgressHistoryNM |> save(joinpath(output_dir_name, "loss_NM.csv"), header = true)

res_all = DataFrame(name = [p_keys_opt...], real = x₀, BBO = x_opt, NM = x_opt_NM)
res_all |> save(output_dir_name * "/res_NM.csv", header = true)
