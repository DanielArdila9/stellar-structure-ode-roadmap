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



# === Grid helpers ===
loggrid(a, b; n) = exp10.(range(log10(a), log10(b), length=n))
lingrid(a, b; n) = collect(range(a, b, length=n))

# === Choose parameter grids ===
# K often varies over decades → use log spacing
Ks  = loggrid(5e13, 5e14; n=7)       # 7 values from 5e13 to 5e14
γs  = [5/3]                          # keep gamma fixed
# ρc can be linear or log; start wider than before
ρcs = lingrid(20.0, 400.0; n=10)     # 10 values: 20 → 400 g/cm^3
# (or) ρcs = loggrid(20.0, 400.0; n=10)

# === Run the sweep ===
results = sweep_modelA(
    Ks, γs, ρcs;
    # solver/tolerances as before:
    solver = Tsit5(), abstol=1e-10, reltol=1e-8,
    # dense radial sampling (points per decade in r):
    saveat = nothing, points_per_decade = 100,
    # parallel across parameter combos:
    use_threads = true
)

# === Quick sanity: how many runs, any failures? ===
@info "runs completed" length(results)
# Optional: view a compact table

# === Plots ===
# Density overlays (log y):
plot_overlaid(results; which=:rho, logρ=true,
              lbl = res -> "ρc=$(res.params.ic.ρc), K=$(res.params.phys.K)")
# Mass overlays:
plot_overlaid(results; which=:mass,
              lbl = res -> "ρc=$(res.params.ic.ρc), K=$(res.params.phys.K)")

# Mass–Radius by K (solar units, clearer physically)
const Msol = 1.98847e33; const Rsol = 6.957e10
Ksuniq = unique([res.params.phys.K for res in results])
p = plot(title="Mass–Radius by K", xscale=:log10, yscale=:log10,
         xlabel="R [R☉]", ylabel="M★ [M☉]")
for K in Ksuniq
    sel = sort([res for res in results if res.params.phys.K == K], by = r -> r.R)
    scatter!(p, [r.R/Rsol for r in sel], [r.Mstar/Msol for r in sel];
             markersize=6, markerstrokewidth=0, label="K=$(K)")
    plot!(p, [r.R/Rsol for r in sel], [r.Mstar/Msol for r in sel]; lw=2, label="")
end
display(p)
