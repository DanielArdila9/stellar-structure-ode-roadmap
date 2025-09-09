# Physics, state definitions, parameters, and RHS for "Model A".
# No `module` here; this file is included into `module StellarModelA`.

############################
# Parameter containers
############################

struct PhysParams
    G::Float64     # gravitational constant (cgs)
    K::Float64     # polytropic constant in P = K * ρ^γ
    γ::Float64     # adiabatic index (e.g., 5/3 for n=1.5 polytrope)
end

struct ICParams
    r0::Float64    # small start radius to avoid r=0 singularity (cm)
    ρc::Float64    # central density (g/cm^3)
end

struct StopParams
    ε::Float64     # density threshold to define the surface (≈0)
    rmax::Float64  # safety ceiling for the domain (cm)
end

struct SolveParams
    solver         # e.g., Tsit5()
    abstol::Float64
    reltol::Float64
    saveat         # either `nothing` or a vector/range of r at which to save
end

############################
# Result container
############################

struct ModelAResult
    R::Float64                # surface radius where ρ ≈ ε
    Mstar::Float64            # M(R)
    params::NamedTuple        # metadata (phys, ic, stop, numerics)
    r::Vector{Float64}        # sampled radii
    M::Vector{Float64}        # M(r) at `r`
    ρ::Vector{Float64}        # ρ(r) at `r`
end

############################
# Regularized initial conditions (helper)
############################

"""
    regularized_u0(ic::ICParams)

Returns the initial state `u0 = [M(r0), ρ(r0)]` using the regularity expansion:
M(r0) ≈ (4π/3)*ρc*r0^3, ρ(r0) ≈ ρc
"""
function regularized_u0(ic::ICParams)
    M0 = (4π/3) * ic.ρc * ic.r0^3
    return [M0, ic.ρc]
end

############################
# RHS for Model A
############################

"""
    modelA!(du, u, phys::PhysParams, r)

State u = [M, ρ]. Equations:
M'(r)  = 4π r^2 ρ
ρ'(r)  = -(G M)/(Kγ r^2) * ρ^(2-γ)
"""
function modelA!(du, u, phys::PhysParams, r)
    @inbounds begin
        M, ρ = u
        G, K, γ = phys.G, phys.K, phys.γ
        du[1] = 4π * r^2 * ρ
        # Guard: avoid division by ~0 if r is extremely small (shouldn't happen if r0>0)
        du[2] = - (G * M) / (K * γ) * ρ^(2 - γ) / (r^2)
    end
    return nothing
end
