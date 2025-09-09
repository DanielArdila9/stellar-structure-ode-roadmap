# Tiny plotting helpers (optional). Included inside `module StellarModelA`.

# Tiny plotting helpers (optional). Included inside `module StellarModelA`.

function _ensure_plots()
    if !isdefined(@__MODULE__, :Plots)
        @eval using Plots
    end
    return nothing
end

raw"""
    plot_profiles(res::ModelAResult; logρ=false) -> (pρ, pM)
"""
function plot_profiles(res::ModelAResult; logρ::Bool=false)
    _ensure_plots()
    rr, MM, rrho = res.r, res.M, res.ρ

    pρ = Plots.plot(rr, rrho;
        lw=2, label="ρ(r)",
        xlabel="r [cm]", ylabel="ρ [g cm⁻³]",
        yscale = (logρ ? :log10 : :identity),
        title = "Density profile"
    )

    pM = Plots.plot(rr, MM;
        lw=2, label="M(r)",
        xlabel="r [cm]", ylabel="M [g]",
        title = "Enclosed mass"
    )

    display(pρ); display(pM)
    return pρ, pM
end

raw"""
    plot_overlaid(results::Vector{ModelAResult};
                  which=:rho, logρ=false,
                  lbl = res -> "ρc=$(res.params.ic.ρc)") -> p
"""
function plot_overlaid(results::Vector{ModelAResult};
                       which::Symbol=:rho, logρ::Bool=false,
                       lbl = res -> "ρc=$(res.params.ic.ρc)")
    _ensure_plots()

    if which === :rho
        p = Plots.plot(title="Density profiles",
                       xlabel="r [cm]", ylabel="ρ [g cm⁻³]",
                       yscale = (logρ ? :log10 : :identity))
        for res in results
            rr, rrho = res.r, res.ρ
            Plots.plot!(p, rr, rrho; lw=2, label=lbl(res))
        end
    elseif which === :mass
        p = Plots.plot(title="Mass profiles",
                       xlabel="r [cm]", ylabel="M [g]")
        for res in results
            rr, MM = res.r, res.M
            Plots.plot!(p, rr, MM; lw=2, label=lbl(res))
        end
    else
        error("which must be :rho or :mass")
    end

    display(p)
    return p
end
