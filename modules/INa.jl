module INa
    using Parameters
    include("../src/functions/ina_dynamics.jl")
    include("../src/functions/algebraic.jl")
    include("../src/functions/rates.jl")
    
    using DifferentialEquations
    using Sundials
    include("../src/problems/ina_current.jl")
end