# Build ODEProblem, callbacks, and run the solve.
# This file is included inside `module StellarModelA`

############################
# Surface event (callback)
############################

"""
    surface_callback(ε::Float64)

ContinuousCallback que detiene la integración cuando ρ(r) cruza ε (>0).
Asume estado u = [M, ρ].
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

Construye un `ODEProblem` para Model A en (r0, rmax) con u0 regularizado.
Pasa `phys` como `p` para evitar closures.
"""
function build_problem(phys::PhysParams, ic::ICParams, stop::StopParams)
    # guardas rápidas (ayudan a fallar temprano)
    @assert phys.K  > 0                "K must be > 0"
    @assert phys.γ  > 1                "γ must be > 1"
    @assert ic.r0   > 0                "r0 must be > 0"
    @assert stop.rmax > ic.r0          "rmax must be > r0"
    @assert ic.ρc   > stop.ε           "ρc must be > ε"

    u0    = regularized_u0(ic)
    rspan = (ic.r0, stop.rmax)

    # Usa directamente modelA!(du,u,p,r) con p=phys
    prob = ODEProblem(modelA!, u0, rspan, phys;
                      # evita que el solver avance a estados físicamente inválidos
                      isoutofdomain = (u,p,r)->(u[2] < 0.0))
    return prob
end

############################
# Integrate once (single run)
############################

"""
    integrate_modelA(phys::PhysParams, ic::ICParams,
                     stop::StopParams, num::SolveParams) :: ModelAResult

Integra hacia afuera desde r0 hasta que ρ ≈ ε (superficie) o rmax.
Devuelve `ModelAResult` con (R, M(R)) y perfiles muestreados.

Política de muestreo:
- Si `num.saveat === nothing`, se resuelve ligero y luego se **post-muestrea en malla log**.
- Si `num.saveat` es un vector, se usa directamente durante la integración.
"""
function integrate_modelA(phys::PhysParams, ic::ICParams,
                          stop::StopParams, num::SolveParams)::ModelAResult

    prob = build_problem(phys, ic, stop)
    cb   = surface_callback(stop.ε)

    # ¿Trae saveat explícito?
    has_saveat = !(num.saveat === nothing)

    # Resolver
    sol = has_saveat ?
        solve(prob, num.solver; callback=cb,
              abstol=num.abstol, reltol=num.reltol,
              save_everystep=false, saveat=num.saveat) :
        solve(prob, num.solver; callback=cb,
              abstol=num.abstol, reltol=num.reltol,
              save_everystep=false)

    # Radio y estado final
    R  = sol.t[end]
    uR = sol.u[end]
    Mstar, ρR = uR[1], uR[2]

    # Perfilería: usa directamente lo guardado o post-muestreo log
    if has_saveat
        r_samp = sol.t
        M_samp = getindex.(sol.u, 1)
        ρ_samp = getindex.(sol.u, 2)
    else
        # malla logarítmica adaptativa: puntos_por_decada * décadas
        decades = max(1.0, log10(R/ic.r0))
        puntos_por_decada = 80                 # ajustable
        npts = Int(clamp(round(decades * puntos_por_decada), 200, 4000))
        r_samp = exp10.(range(log10(ic.r0), log10(R), length=npts))
        # usa interpolación del sol
        M_samp = [sol(r)[1] for r in r_samp]
        ρ_samp = [sol(r)[2] for r in r_samp]
    end

    # Estado/metadata
    status_str = try
        # ReturnCode.Terminated si cortó por callback
        (sol.retcode == ReturnCode.Terminated) ? "surface" : string(sol.retcode)
    catch
        # compat si ReturnCode no está en scope
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
