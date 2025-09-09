# examples/demo_modelA.jl
Pkg.add("StellarModelA")
using StellarModelA, OrdinaryDiffEq   # the module is found via the project

# Define parameters
phys = PhysParams(6.6743e-8, 2.48e14, 5/3)    # G, K, γ
ic   = ICParams(1e4, 50.0)                    # r0, ρc
stop = StopParams(1e-6, 1.5e11)               # ε, rmax
num  = SolveParams(Tsit5(), 1e-10, 1e-8, nothing)  # no saveat → post-sample

# Run a single integration
res = integrate_modelA(phys, ic, stop, num)

# Print summary
println("R ≈ ", res.R, " cm")
println("M* ≈ ", res.Mstar, " g")

# Optional: plot (only if you included Plotting.jl)
using Plots
plot(res.r, res.ρ; xlabel="r [cm]", ylabel="ρ [g cm⁻³]", yscale=:log10, lw=2, label="ρ(r)")
plot!(res.r, res.M; xlabel="r [cm]", ylabel="M [g]", lw=2, label="M(r)")
