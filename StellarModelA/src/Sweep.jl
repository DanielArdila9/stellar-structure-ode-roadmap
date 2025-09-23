# Parameter grid orchestration (sweeps over K, γ, ρc, etc.)
# Included inside `module StellarModelA`

raw"""
    sweep_modelA(Ks, γs, ρcs;
                 G=6.6743e-8, r0=1e4, ε=1e-6, rmax=1.5e11,
                 solver=Tsit5(), abstol=1e-10, reltol=1e-8,
                 saveat=nothing, points_per_decade=80,
                 use_threads::Bool=false)

Cartesian sweep over Ks × γs × ρcs. If `saveat === nothing`, a log-spaced
radial grid is generated automatically (points_per_decade).
If `use_threads=true`, the sweep runs in parallel and preserves the same
lexicographic order as the input grids. Returns Vector{ModelAResult}.
"""
function sweep_modelA(Ks, γs, ρcs;
    G=6.6743e-8, r0=1e4, ε=1e-6, rmax=1.5e11,
    solver=Tsit5(), abstol=1e-10, reltol=1e-8,
    saveat=nothing, points_per_decade=80,
    use_threads::Bool=false)

    Ks = collect(Ks); γs = collect(γs); ρcs = collect(ρcs)
    @assert !isempty(Ks); @assert !isempty(γs); @assert !isempty(ρcs)
    @assert 0 < r0 < rmax
    @assert ε > 0

    # build save grid if needed
    local_saveat = saveat
    if local_saveat === nothing
        decades = max(1.0, log10(rmax/r0))
        npts = Int(clamp(round(decades*points_per_decade), 200, 4000))
        local_saveat = exp10.(range(log10(r0), log10(rmax), length=npts))
    end
    num = SolveParams(solver, abstol, reltol, local_saveat)

    combos  = collect(Iterators.product(Ks, γs, ρcs))
    results = ModelAResult[]  # collect only successes

    function run_one(K, γ, ρc)
        ρc <= ε && throw(ArgumentError("ρc=$(ρc) must be > ε=$(ε)"))
        phys = PhysParams(G, K, γ)
        ic   = ICParams(r0, ρc)
        stop = StopParams(ε, rmax)
        integrate_modelA(phys, ic, stop, num)
    end

    if use_threads && Base.Threads.nthreads() > 1
        # thread-local buffers to avoid push! contention
        bufs = [ModelAResult[] for _ in 1:Threads.nthreads()]
        Threads.@threads for i in eachindex(combos)
            K, γ, ρc = combos[i]
            tid = Threads.threadid()
            try
                push!(bufs[tid], run_one(K, γ, ρc))
            catch err
                @warn "Sweep case failed" K γ ρc err
            end
        end
        for b in bufs; append!(results, b); end
    else
        for i in eachindex(combos)
            K, γ, ρc = combos[i]
            try
                push!(results, run_one(K, γ, ρc))
            catch err
                @warn "Sweep case failed" K γ ρc err
            end
        end
    end

    return results
end


raw"""
    sweep_summary(results::Vector{ModelAResult})

Return a vector of NamedTuples with key metrics per run:
(:R, :Mstar, :K, :γ, :ρc, :status).
"""
function sweep_summary(results::Vector{ModelAResult})
    return [(;
        R      = res.R,
        Mstar  = res.Mstar,
        K      = res.params.phys.K,
        γ      = res.params.phys.γ,
        ρc     = res.params.ic.ρc,
        rho_mean = 3*res.Mstar / (4π*res.R^3),
        ratio    = res.params.ic.ρc / (3*res.Mstar / (4π*res.R^3)),
        status = get(res.params, :status, "unknown")
    ) for res in results]
end

