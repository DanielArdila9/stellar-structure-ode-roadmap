using OrdinaryDiffEq

# ---------- preset cgs + politropo n=1.5 (γ=5/3) ----------
G  = 6.6743e-8
γ  = 5/3
K  = 2.48e14      # dyn/cm^2 / (g/cm^3)^(5/3)  (calibrado para R~Rsun)
r0 = 1e4         # cm (pequeño, evita r=0)

# ---------- RHS de Model A: u = [M, ρ], p = (G, K, γ) ----------
function modelA!(du, u, p, r)
    M, ρ = u
    G, K, γ = p
    du[1] = 4π * r^2 * ρ
    du[2] = -(G * M) / (K * γ) * ρ^(2 - γ) / r^2
end

# CI regulares en r0
ρc = 8.45         # g/cm^3  (valor "de juguete" que cuadra M~Msun con R~Rsun)
M0 = (4π/3) * ρc * r0^3
u0 = [M0, ρc]   
p  = (G, K, γ)

ε = 1e-6 
condition(u, r, integ) = u[2] - ε    # ρ - ε
affect!(integ) = terminate!(integ)
stop_surface = ContinuousCallback(condition, affect!; rootfind=true)

# Problema e integración (techo rmax holgado; el evento corta en R)
rmax = 1.5e11                   # ~2 R_sun como techo
prob = ODEProblem(modelA!, u0, (r0, rmax), p)
sol = solve(prob, Tsit5(); callback=stop_surface,
            abstol=1e-10, reltol=1e-8,
            save_everystep=false,
            saveat=range(r0, rmax, length=2000))  # <-- dense sampling
# Radio y masa finales

R     = sol.t[end]
Mstar = sol.u[end][1]
println("R ≈ ", R, " cm")
println("M⋆ ≈ ", Mstar, " g")


using Plots

r_vals = sol.t
M_vals = getindex.(sol.u, 1)   # M(r)
ρ_vals = getindex.(sol.u, 2)   # ρ(r)

# Densidad
plot(r_vals, ρ_vals; lw=2, xlabel="r [cm]", ylabel="ρ(r) [g cm⁻³]",
     label="densidad", legend=:topright)

# Masa
plot(r_vals, M_vals; lw=2, xlabel="r [cm]", ylabel="M(r) [g]",
     label="masa acumulada", legend=:bottomright)
