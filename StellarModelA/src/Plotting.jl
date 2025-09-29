# Lazy-load Plots only when needed
function _ensure_plots()
    if !isdefined(@__MODULE__, :Plots)
        @eval using Plots
    end
    return nothing
end

"""
    plot_profiles(res::ModelAResult; logρ=false, norm_r=false)

Plot the density ρ(r) and enclosed mass M(r) from a single ModelAResult.

Arguments:
- `logρ`    : if true, plot density on a log10 scale
- `norm_r`  : if true, normalize radius by the stellar radius R (so r̂ = r/R)
"""
function plot_profiles(res::ModelAResult; logρ::Bool=false, norm_r::Bool=false)
    _ensure_plots()
    rr   = norm_r ? res.r ./ res.R : res.r
    MM   = res.M
    rrho = res.ρ

    xlabel_r = norm_r ? "r̂ = r/R" : "r [cm]"

    # Density profile
    pρ = Plots.plot(rr, rrho;
        lw=2, label="ρ(r)",
        xlabel=xlabel_r, ylabel="ρ [g cm⁻³]",
        yscale=(logρ ? :log10 : :identity),
        title="Density profile"
    )

    # Enclosed mass profile
    pM = Plots.plot(rr, MM;
        lw=2, label="M(r)",
        xlabel=xlabel_r, ylabel="M [g]",
        title="Enclosed mass"
    )

    display(pρ); display(pM)
    return pρ, pM
end


"""
    plot_overlaid(results::Vector{ModelAResult};
                  which=:rho, logρ=false, norm_r=false, lbl=...)

Overlay multiple profiles (density or mass) from a vector of ModelAResult.

Arguments:
- `which`   : choose `:rho` for density or `:mass` for enclosed mass
- `logρ`    : if true, use log10 y-scale for density plots
- `norm_r`  : if true, normalize radius by R (so r̂ = r/R)
- `lbl`     : function that maps a result to a curve label
"""
function plot_overlaid(results::Vector{ModelAResult};
                       which::Symbol=:rho, logρ::Bool=false,
                       norm_r::Bool=false,
                       lbl = res -> "ρc=$(res.params.ic.ρc)")
    _ensure_plots()

    xlabel_r = norm_r ? "r̂ = r/R" : "r [cm]"

    if which === :rho
        # Overlay density profiles
        p = Plots.plot(title="Density profiles",
                       xlabel=xlabel_r, ylabel="ρ [g cm⁻³]",
                       yscale=(logρ ? :log10 : :identity))
        for res in results
            rr   = norm_r ? res.r ./ res.R : res.r
            rrho = res.ρ
            Plots.plot!(p, rr, rrho; lw=2, label=lbl(res))
        end
    elseif which === :mass
        # Overlay mass profiles
        p = Plots.plot(title="Mass profiles",
                       xlabel=xlabel_r, ylabel="M [g]")
        for res in results
            rr, MM = norm_r ? (res.r ./ res.R, res.M) : (res.r, res.M)
            Plots.plot!(p, rr, MM; lw=2, label=lbl(res))
        end
    else
        error("which must be :rho or :mass")
    end

    display(p)
    return p
end
