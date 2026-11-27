# Data Dictionary

- `F`: Root-mean-square residual spectral error (dB) over 0.25 < |f| < 8.5 kHz against pooled drone spectrum.
- `AUROC`: Track-level open-set AUROC averaged equally across Gull, Human, and Grey Heron (K = 20 retained dwells).
- `Pd@1%`: Track detection probability at 1% pooled non-drone false-alarm rate.
- `Pfa@0`: Pooled non-drone false-alarm rate at uncalibrated logit zero.
- `5 ms / 9 ms`: Gaussian two-way amplitude FWHM beamwidth taper used during synthetic rendering.

In the repository copy of `data/frozen_extensions.json`, source artifact keys use filenames with unique suffixes in place of machine-local paths. All recorded values are preserved.
