#!/usr/bin/env julia
# =============================================================
# Milestone 1 – Baseline run for Model A (single model)
# Reads params/default.toml, runs once, saves CSV + plots.
# Output: ./datasets/m1_YYYYMMDD-HHMMSS
# =============================================================

using TOML, Dates, CSV, DataFrames, Printf, Plots

# --- Include module from src ---
include(joinpath(@__DIR__, "..", "src", "StellarModelA.jl"))
using .StellarModelA
using .StellarModelA: integrate_modelA, PhysParams, ICParams, StopParams, SolveParams

# --- Read parameters from TOML ---
params_path = joinpath(@__DIR__, "..", "params", "default.toml")
cfg = TOML.parsefile(params_path)

G      = cfg["physical"]["G"]
K      = cfg["eos"]["K"]
γ      = cfg["eos"]["gamma"]
ρc     = cfg["ic"]["rho_c"]
r0     = cfg["ic"]["r0"]
ε      = cfg["stop"]["epsilon"]
rmax   = cfg["stop"]["rmax"]

algo   = cfg["solver"]["algo"]
abstol = cfg["solver"]["abstol"]
reltol = cfg["solver"]["reltol"]
saveat = cfg["solver"]["saveat"] == "auto" ? nothing : cfg["solver"]["saveat"]

solver = StellarModelA._mk_solver(algo)

# Build parameter structs
phys = PhysParams(G, K, γ)
ic   = ICParams(r0, ρc)
stop = StopParams(ε, rmax)
num  = SolveParams(solver, abstol, reltol, saveat)

# Integrate single run (no extra kwargs)
res = integrate_modelA(phys, ic, stop, num)

# Save outputs 
tag    = Dates.format(now(), "yyyymmdd-HHMMSS")
outdir = joinpath(@__DIR__, "..", "datasets", "m1_" * tag)
mkpath(joinpath(outdir, "plots"))

# Profile CSV
CSV.write(joinpath(outdir, "profile_single.csv"),
          DataFrame(r = res.r, M = res.M, rho = res.ρ))

# Summary CSV
R, Mstar = res.R, res.Mstar
rho_bar  = 3Mstar / (4π * R^3)
ratio    = ρc / rho_bar
CSV.write(joinpath(outdir, "summary_single.csv"),
          DataFrame([(; K=K, gamma=γ, rho_c=ρc, R=R, Mstar=Mstar,
                       rho_bar=rho_bar, ratio=ratio)]))

# Plots (just two, readable)
rhat = res.r ./ res.R
plot(rhat, res.ρ; lw=2, yaxis=:log10,
     xlabel="r̂ = r/R", ylabel="ρ [g cm⁻³]",
     title="Density profile (log)")
savefig(joinpath(outdir, "plots", "rho_single.png"))

plot(rhat, res.M; lw=2, xaxis=:log10,
     xlabel="r̂ = r/R (log)", ylabel="M [g]",
     title="Mass profile (semilog-x)")
savefig(joinpath(outdir, "plots", "mass_single.png"))

println("✅ Milestone 1 COMPLETED (single run): ", outdir)
println("   • summary_single.csv")
println("   • profile_single.csv")
println("   • plots/rho_single.png, plots/mass_single.png")
