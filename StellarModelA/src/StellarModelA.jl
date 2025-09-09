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

end # module
