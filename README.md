# Calcium Oscillation Model (Aspatial)

This repository simulates cytosolic and ER calcium (Ca²⁺) dynamics using an aspatial ODE model based on the Li-Rinzel framework, extended to include ER buffering, SOCE entry, cytoplasmic clearance, and volume scaling between compartments. Parameters are set to match those used in the manuscript. The code reproduces Ca²⁺ oscillations under both wild-type and RTN3OE ER conditions, with ER volumes optionally derived from precomputed 3D network geometries.

To get started, run `model_run.m`, which sets parameters, initializes state variables, solves the ODEs, and plots cytosolic and ER Ca²⁺ over time. The core ODE system is defined in `calcium_ode_cyto_V.m`.

## References
[1] Li, Y. X., & Rinzel, J. (1994). Equations for InsP₃ receptor-mediated [Ca²⁺]ᵢ oscillations derived from a detailed kinetic model: a Hodgkin-Huxley like formalism. *Journal of Theoretical Biology*, 166(4), 461–473. https://doi.org/10.1006/jtbi.1994.1041


