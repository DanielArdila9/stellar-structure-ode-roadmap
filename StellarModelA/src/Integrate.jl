# Build ODEProblem, callbacks, and run the solve.
# This file is included inside `module StellarModelA`

############################
# Surface event (callback)
############################

"""
    surface_callback(ε::Float64)

ContinuousCallback that stops the integration when ρ(r) crosses ε (>0).
Assumes state vector u = [M, ρ].
"""
function surface_callback(ε::Float64)
    @assert ε > 0 "ε must be > 0"
    condition(u, r, integ) = u[2] - ε
    affect!(integ) = terminate!(integ)
    return ContinuousCallback(condition, affect!; rootfind=true)
end

############################
# Build ODEProblem
############################

"""
    build_problem(phys::PhysParams, ic::ICParams, stop::StopParams)

Builds an `ODEProblem` for Model A over (r0, rmax) with a regularized u0.
Passes `phys` as `p` to avoid closures.
"""
function build_problem(phys::PhysParams, ic::ICParams, stop::StopParams)
    # quick guards (fail fast if something is invalid)
    @assert phys.K  > 0                "K must be > 0"
    @assert phys.γ  > 1                "γ must be > 1"
    @assert ic.r0   > 0                "r0 must be > 0"
    @assert stop.rmax > ic.r0          "rmax must be > r0"
    @assert ic.ρc   > stop.ε           "ρc must be > ε"

    u0    = regularized_u0(ic)
    rspan = (ic.r0, stop.rmax)

    # Directly use modelA!(du,u,p,r) with p=phys
    prob = ODEProblem(modelA!, u0, rspan, phys;
                      # prevent the solver from stepping into unphysical states
                      isoutofdomain = (u,p,r)->(u[2] < 0.0))
    return prob
end

############################
# Integrate once (single run)
############################

"""
    integrate_modelA(phys::PhysParams, ic::ICParams,
                     stop::StopParams, num::SolveParams) :: ModelAResult

Integrates outward from r0 until ρ ≈ ε (surface) or rmax.
Returns a `ModelAResult` with (R, M(R)) and sampled profiles.

Sampling policy:
- If `num.saveat === nothing`, solve lightly and **post-sample on a log grid**.
- If `num.saveat` is a vector, use it directly during integration.
"""
function integrate_modelA(phys::PhysParams, ic::ICParams,
                          stop::StopParams, num::SolveParams)::ModelAResult

    prob = build_problem(phys, ic, stop)
    cb   = surface_callback(stop.ε)

    # Does it have an explicit saveat?
    has_saveat = !(num.saveat === nothing)

    # Solve
    sol = has_saveat ?
        solve(prob, num.solver; callback=cb,
              abstol=num.abstol, reltol=num.reltol,
              save_everystep=false, saveat=num.saveat) :
        solve(prob, num.solver; callback=cb,
              abstol=num.abstol, reltol=num.reltol,
              save_everystep=false)

    # Final radius and state
    R  = sol.t[end]
    uR = sol.u[end]
    Mstar, ρR = uR[1], uR[2]

    # Profiles: either use saved steps or log post-sampling
    if has_saveat
        r_samp = sol.t
        M_samp = getindex.(sol.u, 1)
        ρ_samp = getindex.(sol.u, 2)
    else
        # logarithmic adaptive mesh: points_per_decade * decades
        decades = max(1.0, log10(R/ic.r0))
        points_per_decade = 80                 # adjustable
        npts = Int(clamp(round(decades * points_per_decade), 200, 4000))
        r_samp = exp10.(range(log10(ic.r0), log10(R), length=npts))
        # interpolate solution on the log grid
        M_samp = [sol(r)[1] for r in r_samp]
        ρ_samp = [sol(r)[2] for r in r_samp]
    end

    # Status/metadata
    status_str = try
        # ReturnCode.Terminated if it stopped by callback
        (sol.retcode == ReturnCode.Terminated) ? "surface" : string(sol.retcode)
    catch
        # fallback if ReturnCode is not in scope
        ρR ≈ stop.ε ? "surface" : "ended"
    end

    meta = (
        phys   = phys,
        ic     = ic,
        stop   = stop,
        num    = num,
        status = status_str,
        solver = string(typeof(num.solver))
    )

    return ModelAResult(R, Mstar, meta, r_samp, M_samp, ρ_samp)
end
