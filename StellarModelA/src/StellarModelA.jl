module StellarModelA

using OrdinaryDiffEq
using LinearAlgebra

include("ModelA.jl")
include("Integrate.jl")
include("Sweep.jl")       # <- antes de Plotting
include("Plotting.jl")    # <- los docstrings aquí son los que pueden romper

export PhysParams, ICParams, StopParams, SolveParams, ModelAResult
export modelA!, integrate_modelA, surface_callback
export sweep_modelA
export plot_profiles, plot_overlaid   # usa SIN bang (coherente con el helper que te di)

export _mk_solver

function _mk_solver(name::AbstractString)
    name == "Tsit5"  && return Tsit5()
    name == "DP5"    && return DP5()
    name == "Vern7"  && return Vern7()
    name == "Rodas5" && return Rodas5()
    name == "TRBDF2" && return TRBDF2()
    error("Unknown solver algo='$name'")
end
end # module
