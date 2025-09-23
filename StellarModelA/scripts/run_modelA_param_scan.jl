# Runs Model A parameter scan (γ, K), computes R, M⋆, ρc/ρ̄,
# runs r0/tolerance sensitivity, and saves comparison plots.

#!/usr/bin/env julia
# Milestone 2: Parameter scan (γ,K) + diagnostics (r0,tols) + comparison plots
using Plots
using TOML, Dates, Printf, CSV, DataFrames
using OrdinaryDiffEq
include(joinpath(@__DIR__, "..", "src", "StellarModelA.jl"))
using .StellarModelA
using .StellarModelA: PhysParams, ICParams, StopParams, SolveParams, ModelAResult
using .StellarModelA: integrate_modelA, sweep_modelA, plot_overlaid

# -------------------- 0) Leer config base (igual a M1) --------------------
params_path = joinpath(@__DIR__, "..", "params", "default.toml")
cfg = TOML.parsefile(params_path)

G     = cfg["physical"]["G"]
K0    = cfg["eos"]["K"]
γ0    = cfg["eos"]["gamma"]
r0    = cfg["ic"]["r0"]
ρc0   = cfg["ic"]["rho_c"]
ε     = cfg["stop"]["epsilon"]
rmax  = cfg["stop"]["rmax"]

algo_str = cfg["solver"]["algo"]
abstol   = cfg["solver"]["abstol"]
reltol   = cfg["solver"]["reltol"]
saveat   = cfg["solver"]["saveat"] == "auto" ? nothing : cfg["solver"]["saveat"]
ppd      = cfg["output_grid"]["points_per_decade"]

# Grillas para el SCAN (puedes sobreescribir desde TOML si quieres)
Ks         = get(cfg["sweep"], "Ks",     [K0, 2K0, 5K0])
gammas     = get(cfg["sweep"], "gammas", [1.3, 1.5, 5/3])
rho_cs     = get(cfg["sweep"], "rho_cs", [ρc0])   # si quieres varias ρc, añade aquí
use_threads = get(cfg["sweep"], "use_threads", false)

# -------------------- 1) Factory de solver --------------------

solver = StellarModelA._mk_solver(algo_str)

# -------------------- 2) Ejecutar SCAN (γ,K,ρc) --------------------
results = sweep_modelA(Ks, gammas, rho_cs;
    G=G, r0=r0, ε=ε, rmax=rmax,
    solver=solver, abstol=abstol, reltol=reltol,
    saveat=saveat, points_per_decade=ppd,
    use_threads=use_threads
)

# -------------------- 3) Armar salida + métricas M2 --------------------
tag    = Dates.format(now(), "yyyymmdd-HHMMSS")
outdir = joinpath(@__DIR__, "..", "datasets", "m2_"*tag)
mkpath(outdir)

# Tabla resumen con (R, M*, rho_bar, ratio)
rows = NamedTuple[]
for res in results
    phys, ic = res.params.phys, res.params.ic
    R, Mstar = res.R, res.Mstar
    rho_bar = 3*Mstar / (4π * R^3)          # ρ̄
    ratio   = ic.ρc / rho_bar               # ρc/ρ̄
    push!(rows, (; K=phys.K, gamma=phys.γ, rho_c=ic.ρc,
                  R=R, Mstar=Mstar, rho_bar=rho_bar, ratio=ratio))
end
df_summary = DataFrame(rows)
sort!(df_summary, [:gamma, :K, :rho_c])
CSV.write(joinpath(outdir, "summary_scan.csv"), df_summary)

# -------------------- 4) Figuras comparativas (ρ y M) --------------------
# Etiquetas compactas por (γ,K,ρc)
lblfun = r -> @sprintf("γ=%.3g, K=%.2e, ρc=%.2e", r.params.phys.γ, r.params.phys.K, r.params.ic.ρc)

# Densidad: log y r normalizado
p_rho  = plot_overlaid(results; which=:rho,  logρ=false,  norm_r=true, lbl=lblfun)
p_mass = plot_overlaid(results; which=:mass, logρ=false, norm_r=true, lbl=lblfun)

# Guardar
plotsdir = joinpath(outdir, "plots"); mkpath(plotsdir)
Base.invokelatest(Plots.savefig, p_rho,  joinpath(plotsdir, "rho_overlaid.png"))
Base.invokelatest(Plots.savefig, p_mass, joinpath(plotsdir, "mass_overlaid.png"))

# -------------------- 5) Sensibilidad: r0 y tolerancias --------------------
# Elige un (K,γ,ρc) “representativo” para la prueba de estabilidad
K_ref  = first(Ks); γ_ref = first(gammas); ρc_ref = first(rho_cs)

r0_grid     = [r0/10, r0, 10r0]                 # puedes ampliar: [1e3, 3e3, 1e4, 3e4, 1e5]
abstol_grid = [1e-8, 1e-10, 1e-12]
reltol_grid = [1e-6, 1e-8, 1e-10]

sens_rows = NamedTuple[]
for r0v in r0_grid, ab in abstol_grid, re in reltol_grid
    # Construir structs explícitos para integrate_modelA
    phys = PhysParams(G, K_ref, γ_ref)
    ic   = ICParams(r0v, ρc_ref)
    stop = StopParams(ε, rmax)
    # si no pasamos saveat, Integrate.jl crea una malla log posmuestreo
    local_saveat = (saveat === nothing) ? nothing : saveat
    num = SolveParams(solver, ab, re, local_saveat)

    res = integrate_modelA(phys, ic, stop, num)
    R, Mstar = res.R, res.Mstar
    rho_bar  = 3*Mstar / (4π * R^3)
    ratio    = ρc_ref / rho_bar

    push!(sens_rows, (; r0=r0v, abstol=ab, reltol=re, R=R, Mstar=Mstar,
                       rho_bar=rho_bar, ratio=ratio, npts=length(res.r)))
end
df_sens = DataFrame(sens_rows)
sort!(df_sens, [:r0, :abstol, :reltol])
CSV.write(joinpath(outdir, "sensitivity_r0_tols.csv"), df_sens)

# -------------------- 6) Mensaje final --------------------
println("✅ Milestone 2 listo en: ", outdir)
println("   • summary_scan.csv           → (K, γ, ρc) ↦ (R, M*, ρ̄, ρc/ρ̄)")
println("   • sensitivity_r0_tols.csv    → estabilidad en R, M* frente a r0 y tolerancias")
println("   • plots/rho_overlaid.png     → ρ(r) comparativo (log, r̂=r/R)")
println("   • plots/mass_overlaid.png    → M(r) comparativo (r̂=r/R)")
