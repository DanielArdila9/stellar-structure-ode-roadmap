# run_sweep.jl

include("src/StellarModelA.jl")
using .StellarModelA
using OrdinaryDiffEq
using Plots

# Define parameter grids
Ks   = [1e14, 2e14, 5e14]       # try three values of K
γs   = [5/3]                    # fix γ
ρcs  = [50.0, 100.0, 200.0]     # vary central density

# Run sweep
results = sweep_modelA(Ks, γs, ρcs)

println("Number of models: ", length(results))
for (i, res) in enumerate(results)
    println("Run $i: R ≈ $(res.R) cm, M* ≈ $(res.Mstar) g, ρc = $(res.params.ic.ρc)")
end

# Plot all density profiles
plot_overlaid(results; which=:rho, logρ=true,
              lbl = res -> "ρc=$(res.params.ic.ρc), K=$(res.params.phys.K)")

# Plot all mass profiles
plot_overlaid(results; which=:mass,
              lbl = res -> "ρc=$(res.params.ic.ρc), K=$(res.params.phys.K)")

R      = [res.R      for res in results]
Mstar  = [res.Mstar  for res in results]

scatter(R, Mstar;
    xscale=:log10, yscale=:log10,
    xlabel="R [cm]", ylabel="M★ [g]",
    markerstrokewidth=0, label="sweep",
    title="Mass–Radius relation")

using Plots

Ks = unique([res.params.phys.K for res in results])

p = plot(title="Mass–Radius by K",
         xscale=:log10, yscale=:log10,
         xlabel="R [cm]", ylabel="M★ [g]")

for K in Ks
    sel = [res for res in results if res.params.phys.K == K]
    sel = sort(sel, by = r -> r.R)

    scatter!(p, [r.R for r in sel], [r.Mstar for r in sel];
             markersize=6, markerstrokewidth=0, label="K=$(K)")
    plot!(p, [r.R for r in sel], [r.Mstar for r in sel]; lw=2, label="")
end

display(p)          # <- needed in scripts
# savefig(p, "MR_by_K.png")  # optional

const Msol = 1.98847e33   # gram
const Rsol = 6.957e10     # cm

using Plots

Ks = unique([res.params.phys.K for res in results])

p = plot(title="Mass–Radius by K",
         xscale=:log10, yscale=:log10,
         xlabel="R [R☉]", ylabel="M★ [M☉]")

for K in Ks
    sel = [res for res in results if res.params.phys.K == K]
    sel = sort(sel, by = r -> r.R)

    scatter!(p, [r.R / Rsol for r in sel],
                [r.Mstar / Msol for r in sel];
                markersize=6, markerstrokewidth=0,
                label="K=$(K)")
    plot!(p, [r.R / Rsol for r in sel],
             [r.Mstar / Msol for r in sel];
             lw=2, label="")
end

display(p)



##NEW

using Plots
include("src/StellarModelA.jl")
using .StellarModelA

# --- helpers for grids ---
loggrid(a, b; n) = exp10.(range(log10(a), log10(b), length=n))
lingrid(a, b; n) = collect(range(a, b, length=n))

# --- parameter grids ---
# K spans decades → log spacing
Ks  = loggrid(5e13, 5e14; n=7)

# Common polytropic gammas (γ = 1 + 1/n):
#   n=3   → γ=4/3  (relativistic degenerate / radiation-dominated)
#   n=2.5 → γ=1.4
#   n=2   → γ=1.5
#   n=1.5 → γ=5/3  (ideal monoatomic / non-relativistic degenerate)
#   n=1   → γ=2
γs  = Float64[4//3, 7//5, 3//2, 5//3, 2//1]  # == [1.333..., 1.4, 1.5, 1.666..., 2.0]

# Central densities — widen a bit
ρcs = loggrid(20.0, 400.0; n=10) # or: loggrid(20.0, 400.0; n=10)

# --- run sweep (dense radial sampling, threaded) ---
results = sweep_modelA(
    Ks, γs, ρcs;
    solver = Tsit5(), abstol=1e-10, reltol=1e-8,
    saveat = nothing, points_per_decade = 100,
    use_threads = true
)

@info "Runs completed" length(results)

# --- plot overlays ---
plot_overlaid(results; which=:rho, logρ=true,
              lbl = res -> "γ=$(round(res.params.phys.γ, digits=3)), K=$(res.params.phys.K)")
plot_overlaid(results; which=:mass,
              lbl = res -> "γ=$(round(res.params.phys.γ, digits=3)), K=$(res.params.phys.K)")

# --- Mass–Radius grouped by γ in solar units ---
const Msol = 1.98847e33
const Rsol = 6.957e10

function plot_MR_by_gamma(results)
    p = plot(title="Mass–Radius by γ",
             xscale=:log10, yscale=:log10,
             xlabel="R [R☉]", ylabel="M★ [M☉]")

    # unique γ with stable sorting
    gammas = sort(unique([res.params.phys.γ for res in results]))
    for γ in gammas
        sel = sort([r for r in results if isapprox(r.params.phys.γ, γ; atol=1e-12)],
                   by = r -> r.R)
        scatter!(p, [r.R/Rsol for r in sel], [r.Mstar/Msol for r in sel];
                 markersize=6, markerstrokewidth=0,
                 label="γ=$(round(γ, digits=3))")
        plot!(p, [r.R/Rsol for r in sel], [r.Mstar/Msol for r in sel];
              lw=2, label="")
    end
    display(p)
    return p
end

plot_MR_by_gamma(results)
