# MathBio Warsaw

Materials for the presentation **"How Does a Hydra Decide Where to Grow a Head?"** — chemical and mechanical routes to pattern formation, covering regeneration, Turing instability, and mechanochemical models.

## Structure

- `Latex/` — Beamer slide deck source
  - `main.tex` — entry point; sections are toggled via `\input{...}` lines
  - `sections/` — individual sections (hydra biology, reaction-diffusion, local dynamics, Turing instability, Gierer–Meinhardt, simulations, mechanochemistry)
  - `figures/` — images and diagrams used in the slides
  - `macros.tex`, `preamble.tex` — shared LaTeX macros and preamble
  - `references.bib` — bibliography
- `Julia/` — supporting numerical simulations
  - `competition_phase_portrait.jl`
  - `gierer_meinhardt_example.jl`
  - `stability_phase_portraits.jl`

## Building

Compile `Latex/main.tex` with a standard LaTeX toolchain (e.g. `latexmk` or `pdflatex`).
