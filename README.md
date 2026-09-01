# MathBio Warsaw

Teaching materials for **“How Does a Hydra Decide Where to Grow a Head?”**, a lecture on chemical and mechanical routes to biological pattern formation. The repository combines a LaTeX/Beamer slide deck, Julia scripts used to generate numerical illustrations, and an interactive reaction–diffusion solver.

The materials were prepared for the [MAT-BIO Exchange Summer School](https://matbio.natur.cuni.cz/mat-bio-exchange), a five-day international school held in Warsaw on 31 August–4 September 2026. The school introduces students to mathematical modelling in biology through lectures, discussion, collaboration, and hands-on exploration.

## Lecture contents

1. **The biological puzzle** — Hydra as a model organism, regeneration, symmetry breaking, and grafting experiments.
2. **Reaction–diffusion approach** — chemical patterning, diffusion, local reaction dynamics, phase portraits, and stability.
3. **Turing instability** — diffusion-driven instability, linearization, spatial eigenmodes, unstable bands, and mode selection.
4. **Gierer–Meinhardt model** — steady states, local stability, activation thresholds, and a worked Turing-instability example.
5. **Growth, injury, and source density** — live numerical experiments on domain size, wounding, and spatial tissue heterogeneity.
6. **Mechanochemical pattern formation** — mechanical feedback, an elastic spring model, stationary states, bifurcations, and oscillations.

## Repository structure

```text
Latex/
├── main.tex             # Beamer entry point
├── preamble.tex         # packages, fonts, colours, and slide theme
├── macros.tex           # reusable slide and navigation macros
├── references.bib       # bibliography
├── sections/            # lecture sections, A–G in presentation order
├── figures/             # one folder per section file, plus common/
└── docs/                # supplementary source material

Julia/
├── competition_phase_portrait.jl
├── stability_phase_portraits.jl
├── gierer_meinhardt_example.jl
└── Solver/              # interactive 1D reaction–diffusion application
```

`Latex/main.tex` controls the lecture order through `\input{...}` statements. Individual sections can be excluded temporarily by commenting out the corresponding input line.

## LaTeX slides

The presentation uses Beamer with `fontspec` and `unicode-math`, so it must be compiled with LuaLaTeX. From the repository root:

```bash
cd Latex
latexmk -lualatex main.tex
```

The generated deck is written to `Latex/main.pdf`. The configured fonts are Arial and Fira Math.

## Julia figure scripts

The standalone scripts use `Plots.jl` and write their output directly to the corresponding directories under `Latex/figures/`:

- `competition_phase_portrait.jl` — phase portrait for a two-species competition system;
- `stability_phase_portraits.jl` — stable and unstable trajectory comparisons;
- `gierer_meinhardt_example.jl` — nullclines, phase portraits, and spatial-mode growth rates for the Gierer–Meinhardt model.

Run a script from the repository root, for example:

```bash
julia Julia/gierer_meinhardt_example.jl
```

## Interactive solver

`Julia/Solver/` is an interactive solver for one-dimensional reaction–diffusion and mechanochemical models. It discretizes space with a sparse finite-difference Laplacian and integrates the resulting method-of-lines system with the adaptive `TRBDF2` method from DifferentialEquations.jl. The application supports:

- Gierer–Meinhardt and mechanochemical model families loaded from `Julia/Solver/Models/`;
- Neumann and periodic boundary conditions;
- live parameter and time-step control;
- localized or random perturbations of the solution;
- saving and restoring the current simulation state;
- splitting, merging, reordering, and deleting domain segments.

Install the solver dependencies and start the QML interface:

```bash
cd Julia/Solver
julia install_dependencies.jl
julia --threads=auto --project=. main_qml.jl
```

The previous GLMakie interface remains available:

```bash
julia --threads=auto --project=. main.jl
```

New models can be added as `.jl` files under `Julia/Solver/Models/`; each file must evaluate to a valid `ModelSpec` or `RDModel` definition.
