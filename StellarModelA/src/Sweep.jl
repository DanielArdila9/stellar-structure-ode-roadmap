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

    # Normalize inputs to vectors (accept any iterable)
    Ks  = collect(Ks);  γs = collect(γs);  ρcs = collect(ρcs)
    @assert !isempty(Ks)  "Ks must be non-empty"
    @assert !isempty(γs)  "γs must be non-empty"
    @assert !isempty(ρcs) "ρcs must be non-empty"
    @assert r0 > 0 && rmax > r0 "Require 0 < r0 < rmax"
    @assert ε > 0 "ε must be > 0"

    # Build default log saveat if none provided
    local_saveat = saveat
    if local_saveat === nothing
        decades = max(1.0, log10(rmax / r0))
        npts = Int(clamp(round(decades * points_per_decade), 200, 4000))
        local_saveat = exp10.(range(log10(r0), log10(rmax), length=npts))
    end
    num = SolveParams(solver, abstol, reltol, local_saveat)

    # Cartesian product & prealloc (deterministic order)
    combos  = collect(Iterators.product(Ks, γs, ρcs))
    results = Vector{ModelAResult}(undef, length(combos))

    # Single run (with guard for ε >= ρc)
    function run_one(K, γ, ρc)
        if ρc <= ε
            throw(ArgumentError("ρc=$(ρc) must be > ε=$(ε)"))
        end
        phys = PhysParams(G, K, γ)
        ic   = ICParams(r0, ρc)
        stop = StopParams(ε, rmax)
        integrate_modelA(phys, ic, stop, num)
    end

    if use_threads && Base.Threads.nthreads() > 1
        Threads.@threads for i in eachindex(combos)
            K, γ, ρc = combos[i]
            try
                results[i] = run_one(K, γ, ρc)
            catch err
                @warn "Sweep case failed" K γ ρc err
                results[i] = nothing
            end
        end
    else
        for i in eachindex(combos)
            K, γ, ρc = combos[i]
            try
                results[i] = run_one(K, γ, ρc)
            catch err
                @warn "Sweep case failed" K γ ρc err
                results[i] = nothing
            end
        end
    end

    # Drop failed entries (if any)
    return [res for res in results if res !== nothing]
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
        status = get(res.meta, :status, "unknown")
    ) for res in results]
end
