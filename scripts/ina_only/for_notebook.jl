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

config_name = abspath("scripts/configs/config_all.json")
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

INa.compute_algebraic!(u₀, p_dict, a)
INa.compute_rates!(du, u, p_dict, 0.0, a)

prob = create_ode_problem(du, u, a, config, tspan, cb_step_v)

solve_kwargs_default = (; reltol, abstol, solver, saveat, dt)
sol = solve_model(prob, (;), solve_kwargs_default);


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

array_length = 100000

data_true = similar(zeros(array_length))
calculate_data_segmented!(data_true, p, prob, solve_kwargs_default)
plot(data_true)


