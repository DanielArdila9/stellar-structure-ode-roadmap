# Build ODEProblem, callbacks, and run the solve.
# No `module` here; this file is included into `module StellarModelA`.

############################
# Surface event (callback)
############################

"""
    surface_callback(ε::Float64)

Returns a ContinuousCallback that stops integration when ρ(r) crosses ε.
Assumes the state is u = [M, ρ].
"""
function surface_callback(ε::Float64)
    condition(u, r, integ) = u[2] - ε
    affect!(integ) = terminate!(integ)
    return ContinuousCallback(condition, affect!; rootfind=true)
end

############################
# Build ODEProblem
############################

"""
    build_problem(phys::PhysParams, ic::ICParams, stop::StopParams)

Builds an ODEProblem for Model A on (r0, rmax) with initial u0.
"""
function build_problem(phys::PhysParams, ic::ICParams, stop::StopParams)
    u0    = regularized_u0(ic)
    rspan = (ic.r0, stop.rmax)
    # Wrap RHS so it matches OrdinaryDiffEq signature f!(du, u, p, r)
    f!(du, u, p, r) = modelA!(du, u, phys, r)
    prob = ODEProblem(f!, u0, rspan, nothing;
                      isoutofdomain = (u,p,r)-> (u[2] < 0.0)) # forbid negative density
    return prob
end

############################
# Integrate once (single run)
############################

"""
    integrate_modelA(phys::PhysParams, ic::ICParams,
                     stop::StopParams, num::SolveParams)

Integrates outward from r0 until ρ ≈ ε (surface), or until rmax.
Returns ModelAResult with (R, M(R)) and sampled profiles.

Saving policy:
- If `num.saveat === nothing`, solve lightly, then post-sample on a uniform grid.
- Else, use `saveat` directly for output sampling.
"""
function integrate_modelA(phys::PhysParams, ic::ICParams,
                          stop::StopParams, num::SolveParams)::ModelAResult

    # Sanity checks
    if ic.ρc <= stop.ε
        error("ic.ρc (central density) must be > stop.ε. Got ρc=$(ic.ρc), ε=$(stop.ε)")
    end

    prob = build_problem(phys, ic, stop)
    cb   = surface_callback(stop.ε)

    # Choose saving policy
    has_saveat = !(num.saveat === nothing)

    sol = if has_saveat
        solve(prob, num.solver; callback=cb,
              abstol=num.abstol, reltol=num.reltol,
              save_everystep=false, saveat=num.saveat)
    else
        solve(prob, num.solver; callback=cb,
              abstol=num.abstol, reltol=num.reltol,
              save_everystep=false)
    end

    # Did we stop by surface, or hit rmax?
    R = sol.t[end]
    uR = sol.u[end]
    ρR = uR[2]

    # If we used no saveat, post-sample on a clean grid
    # If we used no saveat, post-sample on a LOG grid (mejor cerca del centro)
    r_samp, M_samp, ρ_samp = if has_saveat
        r_vec = sol.t
        M_vec = getindex.(sol.u, 1)
        ρ_vec = getindex.(sol.u, 2)
        (r_vec, M_vec, ρ_vec)
    else
        npts = 800
        r_vec = exp10.(range(log10(ic.r0), log10(R), length=npts))
        M_vec = [sol(r)[1] for r in r_vec]
        ρ_vec = [sol(r)[2] for r in r_vec]
        (collect(r_vec), M_vec, ρ_vec)
    end

    # Summaries
    Mstar = uR[1]

    meta = (
        phys = phys,
        ic   = ic,
        stop = stop,
        num  = num,
        status = (ρR ≈ stop.ε ? "surface" : "ended"),
        solver = string(typeof(num.solver))
    )

    return ModelAResult(R, Mstar, meta, r_samp, M_samp, ρ_samp)
end
