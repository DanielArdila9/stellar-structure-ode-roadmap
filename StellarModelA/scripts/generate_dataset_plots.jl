# scripts/make_dataset_from_toml_rhat.jl
using TOML, Dates, Printf, CSV, DataFrames, Plots
include(joinpath(@__DIR__, "..", "src/StellarModelA.jl"))
using .StellarModelA  # exports sweep_modelA, ModelAResult, etc.
using OrdinaryDiffEq
# ---------- 0) Controls ----------

# ---------- 1) Load TOML ----------
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

Nhat = get(cfg["resample"], "Nhat", 512)  # default to 512 if not in TOML

ppd = cfg["output_grid"]["points_per_decade"]

Ks         = get(cfg["sweep"], "Ks",     [K0])
gammas     = get(cfg["sweep"], "gammas", [γ0])
rho_cs     = get(cfg["sweep"], "rho_cs", [ρc0])
use_threads = get(cfg["sweep"], "use_threads", false)

# ---------- 2) Solver factory ----------
function _mk_solver(name::AbstractString)
    name == "Tsit5"  && return Tsit5()
    name == "DP5"    && return DP5()
    name == "Vern7"  && return Vern7()
    name == "Rodas5" && return Rodas5()
    name == "TRBDF2" && return TRBDF2()
    error("Unknown solver algo='$name'")
end
solver = _mk_solver(algo_str)

# ---------- 3) Run sweep ----------
# If saveat === nothing, sweep_modelA builds a log(r) grid using points_per_decade. 
results = sweep_modelA(Ks, gammas, rho_cs;
    G=G, r0=r0, ε=ε, rmax=rmax,
    solver=solver, abstol=abstol, reltol=reltol,
    saveat=saveat, points_per_decade=ppd,
    use_threads=use_threads
)  # returns Vector{ModelAResult} with r, M(r), ρ(r).  

# ---------- 4) Output dir ----------
tag    = Dates.format(now(), "yyyymmdd-HHMMSS")
outdir = joinpath(@__DIR__, "..", "datasets", "set_" * tag)
mkpath(outdir)

# ---------- 5) Simple 1D linear interpolation ----------
function _interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::AbstractVector{<:Real})
    yi = similar(Float64.(xq))
    i = 1
    @inbounds for j in eachindex(xq)
        xv = xq[j]
        while i < length(x) && x[i+1] < xv
            i += 1
        end
        if xv <= x[1]
            yi[j] = y[1]
        elseif xv >= x[end]
            yi[j] = y[end]
        else
            x1, x2 = x[i], x[i+1]
            y1, y2 = y[i], y[i+1]
            t = (xv - x1) / (x2 - x1)
            yi[j] = (1 - t) * y1 + t * y2
        end
    end
    return yi
end

# ---------- 6) Resample each run to exactly Nhat on r̂ ∈ (0,1] ----------
function resample_rhat(res; Nhat::Int=Nhat)
    R, Mstar = res.R, res.Mstar
    r, M, ρ  = res.r, res.M, res.ρ
    ρc       = res.params.ic.ρc

    rhat_q = collect(range(1/Nhat, 1.0; length=Nhat))  # avoid r̂=0
    rq     = rhat_q .* R
    Mq     = _interp1(r, M, rq)
    ρq     = _interp1(r, ρ, rq)

    return (; R, Mstar, ρc,
             r=rq, M=Mq, rho=ρq,
             r_hat=rhat_q, M_hat=Mq ./ Mstar, rho_hat=ρq ./ ρc)
end

# ---------- 7) Build big profiles + summary ----------
rows = DataFrame[]
srows = DataFrame[]  # summary rows

for res in results
    phys, ic = res.params.phys, res.params.ic
    s = resample_rhat(res; Nhat=Nhat)

    push!(rows, DataFrame(;
        K = phys.K, gamma = phys.γ, rho_c = ic.ρc,
        R = s.R, Mstar = s.Mstar,
        r = s.r, M = s.M, rho = s.rho,
        r_hat = s.r_hat, M_hat = s.M_hat, rho_hat = s.rho_hat
    ))

    push!(srows, DataFrame(;
        K = phys.K, gamma = phys.γ, rho_c = ic.ρc,
        R = s.R, Mstar = s.Mstar, npts = Nhat,
        status = get(res.params, :status, "unknown")
    ))
end

profiles = vcat(rows...)
CSV.write(joinpath(outdir, "profiles_resampled.csv"), profiles)

# unique summary per (K,γ,ρc)
summary = unique(vcat(srows...))
CSV.write(joinpath(outdir, "summary.csv"), summary)

println("✅ Dataset listo en: ", outdir)
println("   • profiles_resampled.csv  (one big table, each run has exactly Nhat = $Nhat rows)")
println("   • summary.csv             (one row per (K, γ, ρc))")

# --- Overlaid plots ---
p_rho  = plot_overlaid(results; which=:rho,  logρ=true, norm_r=true, lbl=r->"ρc=$(r.params.ic.ρc)")
p_mass = plot_overlaid(results; which=:mass, logρ=false, norm_r=true, lbl=r->"ρc=$(r.params.ic.ρc)")

plotsdir = joinpath(outdir, "plots")
mkpath(plotsdir)
Base.invokelatest(Plots.savefig, p_rho,  joinpath(plotsdir, "all_density.png"))
Base.invokelatest(Plots.savefig, p_mass, joinpath(plotsdir, "all_mass.png"))