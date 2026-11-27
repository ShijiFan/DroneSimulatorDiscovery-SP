# Reproducibility Package — Near-Axis Rotor Scattering for Open-Set 77 GHz Drone Detection

## Contents
- `simulators/O-b1-dist.py`: Exact distilled six-parameter three-component rotor scattering model ($F = 0.71$ dB).
- `simulators/O-b2.py`: Exact validation-selected near-axis rotor scattering model ($F = 1.14$ dB, track AUROC 0.989).
- `data/W14_summary.json`, `data/W14b_summary.json`: Spectral fitness and cross-taper simulator evaluations.
- `data/W18_summary.json`, `data/W18_tables.json`: Recording-clustered bootstrap open-set metrics and ablations.
- `data/frozen_extensions.json`: Disjoint two-fold threshold transfer and validation-set size results.
- `DATA_DICTIONARY.md`: Field definitions and mapping to manuscript tables.

## Dataset
Measured 77 GHz FMCW recordings are publicly available from the SAAB SIRS benchmark:
Karlsson et al., Zenodo 5896641 (doi:10.5281/zenodo.5896641).
Full synthesis pipeline, training scripts, and per-seed evaluation logs will be released publicly upon acceptance.
