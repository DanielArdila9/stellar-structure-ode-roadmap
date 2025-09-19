# Tiny plotting helpers (optional). Included inside `module StellarModelA`.

# Tiny plotting helpers (optional). Included inside `module StellarModelA`.

function _ensure_plots()
    if !isdefined(@__MODULE__, :Plots)
        @eval using Plots
    end
    return nothing
end
function plot_profiles(res::ModelAResult; logρ::Bool=false, norm_r::Bool=false)
    _ensure_plots()
    rr   = norm_r ? res.r ./ res.R : res.r
    MM   = res.M
    rrho = res.ρ

    xlabel_r = norm_r ? "r̂ = r/R" : "r [cm]"

    pρ = Plots.plot(rr, rrho;
        lw=2, label="ρ(r)",
        xlabel=xlabel_r, ylabel="ρ [g cm⁻³]",
        yscale=(logρ ? :log10 : :identity),
        title="Density profile"
    )

    pM = Plots.plot(rr, MM;
        lw=2, label="M(r)",
        xlabel=xlabel_r, ylabel="M [g]",
        title="Enclosed mass"
    )

    display(pρ); display(pM)
    return pρ, pM
end


function plot_overlaid(results::Vector{ModelAResult};
                       which::Symbol=:rho, logρ::Bool=false,
                       norm_r::Bool=false,
                       lbl = res -> "ρc=$(res.params.ic.ρc)")
    _ensure_plots()

    xlabel_r = norm_r ? "r̂ = r/R" : "r [cm]"

    if which === :rho
        p = Plots.plot(title="Density profiles",
                       xlabel=xlabel_r, ylabel="ρ [g cm⁻³]",
                       yscale=(logρ ? :log10 : :identity))
        for res in results
            rr   = norm_r ? res.r ./ res.R : res.r
            rrho = res.ρ
            Plots.plot!(p, rr, rrho; lw=2, label=lbl(res))
        end
    elseif which === :mass
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
