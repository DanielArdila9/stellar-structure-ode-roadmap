# Stellar Structure ODE Roadmap

Numerical experiments with **stellar structure equations** using ODE solvers and a polytropic equation of state (EOS).  
The goal is to generate CSV profiles and plots for density, pressure, and enclosed mass, stopping integration at the stellar surface (ρ → 0).

---

## Project Structure

```
stellar-structure-ode-roadmap/
├── StellarModelA/       # Julia module implementing Model A
│   ├── modelA.jl        # Core ODE system and integration
│   ├── plots.jl         # Plotting helpers
│   └── utils.jl         # Utility functions (callbacks, etc.)
├── run_modelA.jl        # Julia driver script (to be added)
├── run_modelA.py        # Python wrapper (optional; to be added)
├── params.toml          # Default parameters (K, γ, ρc, tolerances; to be added)
├── output/              # Generated CSVs and plots
├── .gitignore
└── README.md
```

---

## Features (in progress)

- Numerical integration of stellar structure ODEs with polytropic EOS.  
- Termination at the stellar surface using callbacks.  
- Generation of radial profiles:  
  - Mass M(r)  
  - Density ρ(r)  
  - Pressure P(r)  
- Export to CSV and plots (via `Plots.jl`).  
- Parameter grid runs for different values of K, γ, and central density ρc.  

---

## Next Steps

- [ ] Add `params.toml` with default physical and numerical parameters.  
- [ ] Implement Python driver for reproducibility.  
- [ ] Extend to other EOS or models (e.g., Model B).  
