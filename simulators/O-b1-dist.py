"""Parametric family for the slow-time return of a small multirotor's rotating parts.

The signal is the coherent sum of three physical scatterer classes, for each of n_rot rotors
(rate f = f0 * (1 +/- 5 %) random per rotor, random start angle) and each of n_b blades at
angle psi_b(t) = 2 pi f t + phi + 2 pi b / n_b.  Projected line-of-sight phase uses
k0 = 4 pi / lam * cos(beta) (radar elevation beta foreshortens the in-plane rotation).

1. Blades (reference power 1): each blade is a uniform line scatterer from the hub to the tip
   at radius L.  Integrating exp(j k0 r cos psi) over r in [0, L] gives
   exp(j x/2) sinc(x / 2pi), x = k0 L cos psi: a strong specular "flash" when the blade is
   broadside (cos psi ~ 0) plus a weak continuum out to the tip Doppler +/- 2 f0 k0 L / 2.
   No free parameter: this is the fixed backbone of the family.

2. Symmetric hub features (power p_sym relative to blades): motor bell, hub cap and screw
   heads, i.e. a few (3 per blade) point scatterers at very small radius r = r_sym * L * u^3
   (u ~ U(0,1), so most sit close to the axis), each with a random scattering phase.  They are
   visible over the whole revolution, so they fill a narrow, symmetric band of micro-Doppler
   around zero (half-width ~ 2 pi f k0 r_sym L / 2pi).

3. One-sided blade-root / pitch-link features (power p_one relative to blades): 2 point
   scatterers per blade at somewhat larger radius r = r_one * L * (0.1 + 0.9 u^2).  The hub
   housing shadows them over part of the revolution, modelled as an amplitude visibility
   1 + occ * sin(psi + occ_phase): occ = 0 is fully visible, occ -> 1 nearly blanks one
   half-turn.  Because visibility correlates with the sign of the radial velocity
   (d/dt cos psi = -sin psi), the energy lands preferentially on one side of zero Doppler:
   this is the only source of Doppler-spectrum asymmetry, and occ_phase sets which side and
   how strongly (max asymmetry when the shadowed half-turn is the approaching one).

Each class is normalised to its own mean power before scaling, so p_sym and p_one are
power ratios to the blade return and are independent of the random draws.
"""
import numpy as np

PARAMS = [("p_sym", 0.0, 10.0), ("r_sym", 0.005, 0.30), ("p_one", 0.0, 10.0), ("r_one", 0.01, 0.40), ("occ", 0.0, 1.0), ("occ_phase", -3.1416, 3.1416)]

THETA = {'p_sym': 2.156365917241298, 'r_sym': 0.23968545730567517, 'p_one': 6.671447063905026, 'r_one': 0.10958238015687433, 'occ': 1.0, 'occ_phase': 0.8257539242318159}


def rotor_signal(t, L, f0, n_b, n_rot, beta, lam, rng):
    th = THETA
    t = np.asarray(t, dtype=float)
    n = len(t)
    k0 = 4 * np.pi / lam * np.cos(beta)
    out_b = np.zeros(n, dtype=complex)
    out_s = np.zeros(n, dtype=complex)
    out_o = np.zeros(n, dtype=complex)
    for _ in range(int(n_rot)):
        f = f0 * (1 + 0.05 * rng.uniform(-1, 1))
        phi = rng.uniform(0, 2 * np.pi)
        for b in range(int(n_b)):
            psi = 2 * np.pi * f * t + phi + 2 * np.pi * b / n_b
            c = np.cos(psi)
            # 1. blade line scatterer
            x = k0 * L * c
            out_b += np.exp(0.5j * x) * np.sinc(x / (2 * np.pi))
            # 2. symmetric hub points (always visible)
            r_s = th["r_sym"] * L * rng.uniform(0, 1, 3) ** 3
            ph_s = rng.uniform(0, 2 * np.pi, 3)
            out_s += np.exp(1j * (k0 * np.outer(c, r_s) + ph_s)).sum(axis=1)
            # 3. one-sided root points (shadowed by hub housing)
            r_o = th["r_one"] * L * (0.1 + 0.9 * rng.uniform(0, 1, 2) ** 2)
            ph_o = rng.uniform(0, 2 * np.pi, 2)
            vis = 1 + th["occ"] * np.sin(psi + th["occ_phase"])
            out_o += vis * np.exp(1j * (k0 * np.outer(c, r_o) + ph_o)).sum(axis=1)
    pb = np.mean(np.abs(out_b) ** 2) + 1e-30
    ps = np.mean(np.abs(out_s) ** 2) + 1e-30
    po = np.mean(np.abs(out_o) ** 2) + 1e-30
    return out_b + out_s * np.sqrt(th["p_sym"] * pb / ps) + out_o * np.sqrt(th["p_one"] * pb / po)
