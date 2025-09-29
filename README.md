# Stellar Structure ODE Roadmap

Numerical experiments with **stellar structure equations** using ODE solvers and a polytropic equation of state (EOS).

We study a simplified star model (**Model A**) where mass and density evolve as:

$$
\frac{dM}{dr} = 4\pi r^2 \rho,\qquad
\frac{d\rho}{dr} = - \frac{G M}{K \gamma r^2}\, \rho^{2-\gamma}.
$$

Integration starts at a small radius r0 with regularized central conditions and stops automatically at the stellar surface ($\rho \to 0$).

---

## Project Structure

```
stellar-structure-ode-roadmap/
├── StellarModelA/          # Julia module
│   ├── src/
│   │   ├── ModelA.jl       # Physics: parameters, RHS for Model A
│   │   ├── Integrate.jl    # ODEProblem setup, event callback, integration
│   │   ├── Sweep.jl        # Parameter sweeps (K, γ, ρc) + summaries
│   │   ├── Plotting.jl     # Plotting helpers (ρ(r), M(r))
│   │   └── StellarModelA.jl# Main module wrapper
│   ├── datasets/           # Generated CSVs (ignored by git)
│   ├── params/             # TOML parameter files (default configs)
│   └── scripts/
│       ├── run_modelA_baseline.jl      # Milestone 1: baselines
│       └── run_modelA_param_scan.jl    # Milestone 2: scans + diagnostics
├── README.md
└── .gitignore
```

---

## How to Run

### Milestone 1 – Baselines and Data
Run a few test cases to check the solver and generate synthetic datasets:
```bash
julia --project=. scripts/run_modelA_baseline.jl
```
Outputs:
- `profiles_resampled.csv` (fixed-grid profiles of ρ(r), M(r))
- `summary.csv` (one row per parameter set)
- Overlay plots of density and mass

### Milestone 2 – Parameter Scans and Diagnostics
Perform systematic parameter scans and sensitivity tests:
```bash
julia --project=. scripts/run_modelA_param_scan.jl
```
Outputs:
- `summary_scan.csv` → Table with (K, γ, ρc) → (R, M*, ρ̄, ρc/ρ̄)
- `sensitivity_r0_tols.csv` → Stability results vs. r0, abstol, reltol
- `plots/` → Comparison plots of ρ(r) and M(r) across parameter sets

---

## Milestones

**Milestone 1 (Weeks 1–2): Baselines and Data**
- Implemented Model A with event-driven stop.
- Produced clean runs and datasets.
- Verified monotone M(r) growth and physically plausible ρ(r).

**Milestone 2 (Weeks 3–4): Parameter Scans and Diagnostics**
- Scanned γ ∈ {1.3, 1.5, 5/3} and multiple K values.
- Tabulated surface radius R, stellar mass M*, mean density ρ̄, and ratio ρc/ρ̄.
- Performed sensitivity studies on r0 and solver tolerances.
- Generated comparison plots across parameter sets.

---

## Dependencies

This project uses [Julia](https://julialang.org) with:

- [OrdinaryDiffEq.jl](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) (ODE solvers: Tsit5, DP5, Vern7, etc.)
- [Plots.jl](http://docs.juliaplots.org/latest/) (for figures)
- [CSV.jl](https://csv.juliadata.org/) and [DataFrames.jl](https://dataframes.juliadata.org/) (for output tables)
- [TOML.jl](https://github.com/JuliaLang/TOML.jl) (for parameter configs)

Install them in the project environment:
```bash
julia --project=.
julia> ]instantiate
```

---

## Reproducibility

- Default parameters are stored in `params/default.toml`.
- All outputs (CSVs, plots) are timestamped and saved under `StellarModelA/datasets/`.
- Run scripts directly from `scripts/` to reproduce each milestone’s deliverables.

---

Next milestone (3) will add **Model B** (convective, 3 variables: M(r), P(r), T(r)) and comparisons with Model A.
