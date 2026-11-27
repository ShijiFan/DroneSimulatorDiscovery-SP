import numpy as np

RMIN, RMAX = 0.0, 0.10
NSC = 4
AMP = 0.11
TAPER = 4.5
ASYM0, ASYM_SLOPE, ASYM_MAX = 0.35, 5.0, 0.85


def rotor_signal(t, L, f0, n_b, n_rot, beta, lam, rng):
    """Rigid line blades (Martin-Mulgrew) + inboard hub/blade-root scatterers.
    Inboard scatterers at small random radii and random azimuths on the hub give a
    low-Doppler hump; their aspect-dependent RCS (stronger when receding, more so
    further out along the blade root) gives the negative-Doppler asymmetry."""
    out = np.zeros(len(t), dtype=complex)
    k0 = 4 * np.pi / lam * np.cos(beta)
    for k in range(n_rot):
        f = f0 * (1 + 0.05 * rng.uniform(-1, 1))
        phi = rng.uniform(0, 2 * np.pi)
        rho = RMIN + (RMAX - RMIN) * rng.uniform(0, 1, NSC) ** 1.5
        dth = rng.uniform(0, 2 * np.pi, NSC)
        asym = np.minimum(ASYM_MAX, ASYM0 + ASYM_SLOPE * rho)
        for b in range(n_b):
            psi = 2 * np.pi * f * t + phi + 2 * np.pi * b / n_b
            x = k0 * L * np.cos(psi)
            out += np.exp(0.5j * x) * np.sinc(x / (2 * np.pi))
            for r, a, d in zip(rho, asym, dth):
                p = psi + d
                out += AMP * (1 - TAPER * r) * (1 + a * np.sin(p)) * np.exp(1j * k0 * r * L * np.cos(p))
    return out
