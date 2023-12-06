using DifferentialEquations,Statistics, Distributions, Plots, StatsBase
using LinearAlgebra, NLsolve
using DataFrames, CSV

#include("code/example_1/ex1.jl")
#include("code/Trait_Plot/Trait_Plot.jl")
#include("code/example_2/ex2.jl")
#include("code/example_3/ex3.jl")
#include("code/example_4/ex4.jl")
#include("./code/example_5/ex5.jl")
#include("./code/example_6/ex6.jl")
#include("./code/example_7/ex7.jl")
#include("./code/example_9/ex9.jl")
#include("./code/AD_example_1/code.jl")
#include("./code/AD_example_1/PIP_code.jl")
#include("./code/AD_example_1/AD_hill_climb.jl")
#include("./code/AD_example_1/PIP_code_alt.jl")
#include("./code/AD_example_1/within_host_PIP_code.jl")
include("code/AD_2_Species/code.jl")
#include("code/AD_2_Species/draw_PIP.jl")
#include("code/AD_2_Species/analysis.jl")
#include("code/AD_2_Species/draw_TEP.jl")
include("code/AD_2_Species/draw_SelectionGrad.jl")