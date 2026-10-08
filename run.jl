using DifferentialEquations,Statistics, Distributions, Plots, StatsBase
using LinearAlgebra, NLsolve, Roots, Sobol
using DataFrames, CSV, JLD2, PoissonRandom

#--Function and struct definitions--
include("code/auxiliary_funs/bisection.jl")
include("code/solve_cubic_code/solve_cubic.jl")
include("code/AD_2_Species/code.jl")
include("code/AD_2_Species/draw_PIP.jl")
include("code/AD_2_Species/find_host_shifts/find_host_shifts.jl")


#--Scripts to run simulations and draw figures--
#include("code/AD_2_Species/draw_TEP.jl")
include("code/AD_2_Species/stochastic_individual_based_model/continuous_model/approximate_model.jl") #Code to run stochastic simulations and create sample path figures
#include("code/AD_2_Species/PIP_parameter_check/PIP_parameter_check.jl") #Code to run parameter sweep
#include("code/AD_2_Species/interspecific_contact_influence/interspecific_contact_influence.jl") #Code to draw the effect of interspecific contact intensity and host-species similarity
