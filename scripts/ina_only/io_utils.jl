# import Pkg
# Pkg.add("JSON")
using JSON
using DataFrames: DataFrame
using CSV: File as CSVFile
using LabelledArrays
using NamedTupleTools
using ArgParse
using Logging

read_csv(filename) = DataFrame(CSVFile(filename))

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
        help = "a config name"
        required = true
    end
    return parse_args(s)
end

function create_genes_dict_from_config(config)
    genes_dict = Dict()
    for (ec_name, ec) in config["initial"]["experimental_conditions"]
        genes_dict[ec_name] = ec["params"]
    end
    return genes_dict
end

function generate_bounds_gammas_mask_multipliers(genes_dict)
    gene_keys_opt, bounds, mask_multipliers = [], [], []
    for genes in values(genes_dict)
        for (gene_name, gene) in genes
            #println(gene_name)
            bound, mask_multiplier = gene["bounds"], gene["is_multiplier"]
            push!(mask_multipliers, mask_multiplier)
            push!(bounds, bound)
            push!(gene_keys_opt, gene_name)
        end
    end
    LVector(namedtuple(Symbol.(gene_keys_opt), gene_keys_opt)),
    LVector(namedtuple(Symbol.(gene_keys_opt), bounds)),
    LVector(namedtuple(Symbol.(gene_keys_opt), mask_multipliers))
end

function make_config(config_filename)
    config_path = abspath(config_filename)
    config_initial = JSON.parsefile(config_path)
    # config = Dict("filenames" => Dict(name => config_initial[name] for name in ["filename_legend_states",
    #                                                                 "filename_legend_constants",
    #                                                                 "filename_legend_algebraic",
    #                                                                 "filename_protocol",
    #                                                                 ]))
    config = Dict("initial" => JSON.parsefile(config_path))
    config["initial"]["config_path"] = config_path

    runtime = Dict()
    config["runtime"] = runtime
    config["runtime"]["genes_dict"] = create_genes_dict_from_config(config) # how about gamma?
    #     config["runtime"]["constants_dict"] = make_constants_dict_from_config(config)

    df_constants = read_csv(abspath(config["initial"]["filename_legend_constants"]))
    config["runtime"]["constants"] =
        LVector(namedtuple(Symbol.(df_constants.name), df_constants.value))
    config["runtime"]["constants_bounds"] = LVector(
        namedtuple(
            Symbol.(df_constants.name),
            [df_constants.bound_1[k], df_constants.bound_2[k]] for
            k = 1:length(df_constants.bound_2)
        ),
    )

    df_states = read_csv(abspath(config["initial"]["filename_legend_states"]))
    config["runtime"]["states"] =
        LVector(namedtuple(Symbol.(df_states.name), df_states.value))


    for (exp_cond_name, exp_cond) in config["initial"]["experimental_conditions"]
        if exp_cond_name == "common"
            continue
        end

        filename_phenotype = abspath(exp_cond["filename_phenotype"])
        exp_cond["filename_phenotype"] = filename_phenotype
        exp_cond["phenotype"] = read_csv(filename_phenotype)

        if haskey(exp_cond, "filename_sample_weight")
            filename_weight = abspath(exp_cond["filename_sample_weight"])
            exp_cond["filename_sample_weight"] = filename_weight
            exp_cond["sample_weight"] = read_csv(filename_weight).w
        end
    end

    protocol = read_csv(abspath(config["initial"]["filename_protocol"]))
    config["runtime"]["protocol"] = protocol

    model_type = config["initial"]["model_type"]
    config["runtime"]["model_type"] = Symbol(model_type) # Meta.parse(Main.eval(model_name))  # 


    gene_keys_opt, bounds, mask_multipliers =
        generate_bounds_gammas_mask_multipliers(config["runtime"]["genes_dict"])
    config["runtime"]["bounds"] = bounds
    config["runtime"]["mask_multipliers"] = mask_multipliers
    config["runtime"]["gene_keys_opt"] = gene_keys_opt

    return config
end

function make_output_dirs(config)
    time = (Dates.format(Dates.now(), "yy_mm_dd_HHMMSS"))
    output_main_dir = abspath(normpath(config["initial"]["output_folder_name"]))
    output_dir = joinpath(output_main_dir, time)

    if !isdir(output_main_dir)
        @info "Making $output_main_dir"
        mkdir(output_main_dir)
    else
        @info "Output directory $output_main_dir already exist!"
    end

    if !isdir(output_dir)
        @info "Making $output_dir"
        mkdir(output_dir)
    else
        @info "Output directory $output_dir already exist!"
    end
    return output_dir
end

function save_progress_history(fitness_progress_history, output_dir_name, name)
    df_fitness_progress_history = DataFrame(
        numb = [fitness_progress_history[k][1] for k = 1:length(fitness_progress_history)],
        loss = [fitness_progress_history[k][2] for k = 1:length(fitness_progress_history)],
    )

    df_fitness_progress_history |> save(joinpath(output_dir_name, name), header = true)
    return nothing
end
