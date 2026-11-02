#!/usr/bin/env python3
"""P14 / WP-D0': cogenese por conservacao de N_psi (ver P14D0b_DEFINITION.md).
Censo dos 39 estados copiado de proofs/proofs/P_yield_ratio.py (passo D.6)."""
import math
from fractions import Fraction as F
from itertools import combinations_with_replacement as cwr

Q9 = [1, -5, -5, 10, 4, 4]
LAM_REPO, LAM_PAPER = 1.2662e12, 1.2609e12
B_E, A, P = 325 / 1152, F(13, 12), F(5, 6)
YP, M_N = 0.245, 0.9389
NU_VIS = 3 + 3 * (1 - YP / 2)
OBS = 0.1200 / 0.02237


def q9s(t): return sum(Q9[k] for k in t)
def tri(t): return sum(t) % 6
def mRJ(t): return abs(sum(F(2 * k - 5, 2) for k in t if k >= 3))


def stable(t):
    n = len(t); reach = {(0, 0): 1}
    for k in t:
        new = dict(reach)
        for (q, tr), mask in reach.items():
            key = (q + Q9[k], (tr + k) % 6)
            new[key] = new.get(key, 0) | (mask << 1)
        reach = new
    mid = ((1 << n) - 1) & ~1
    return not any(q == 0 and tr in (0, 3) and (mask & mid) for (q, tr), mask in reach.items())


def mass(t, lam=LAM_REPO):
    return lam * math.exp(-(7 + len(t) - float(A) * float(mRJ(t)) ** float(P)) / B_E)


states = []
for N in (9, 12, 15, 18, 21, 27):
    for t in cwr(range(6), N):
        if q9s(t) == 0 and tri(t) in (0, 3) and stable(t):
            states.append((N, t, mass(t)))
assert len(states) == 39, len(states)

print("##### Censo: 39 estados escuros estaveis (H4) #####")
print(f"  {'N':>3} {'m [GeV]':>12} {'m/N [GeV]':>11}  rotulos (contagem por k)")
for N, t, m in sorted(states, key=lambda s: s[2] / s[0]):
    cnt = [t.count(k) for k in range(6)]
    print(f"  {N:3d} {m:12.4e} {m/N:11.4e}  {cnt}  I={tri(t)} mRJ={mRJ(t)}")

print(f"\n  nu_vis = {NU_VIS:.4f} psions visiveis por barion ; m_N = {M_N} GeV")


def omega(m_over_N):
    return NU_VIS * m_over_N / M_N


r1 = min(states, key=lambda s: s[2] / s[0])
r3 = min(states, key=lambda s: s[2])
mN_r2 = sum(s[2] for s in states) / sum(s[0] for s in states)
deg_r1 = [s for s in states if abs(s[2] / s[0] / (r1[2] / r1[0]) - 1) < 1e-9]

print("\n##### Regras de reliquia (fixadas antes) #####")
print(f"  R1 min(m/N): N={r1[0]}, m={r1[2]:.4e} GeV, m/N={r1[2]/r1[0]:.4e}  (degenerescencia: {len(deg_r1)} estados)")
print(f"  R2 democratico: m/N efetivo = {mN_r2:.4e}")
print(f"  R3 min(m): N={r3[0]}, m={r3[2]:.4e} GeV, m/N={r3[2]/r3[0]:.4e}")

print("\n##### Predicao Omega_DM/Omega_b = nu_vis (m/N)/m_N  (comparacao so agora) #####")
print(f"  observado (Planck 2018) = {OBS:.3f}")
for name, mn in (("R1", r1[2] / r1[0]), ("R2", mN_r2), ("R3", r3[2] / r3[0])):
    w = omega(mn)
    print(f"  {name}: {w:.4e}  ({100*(w/OBS-1):+.1f}%)  {'COMPATIVEL (+-10%)' if abs(w/OBS-1) <= 0.10 else ''}")
w1p = omega(r1[2] / r1[0] * LAM_PAPER / LAM_REPO)
print(f"  R1 com Lambda do paper (1.2609e12): {w1p:.4e} ({100*(w1p/OBS-1):+.1f}%)")
print(f"\n  m/N exigido pelo observado: {OBS*M_N/NU_VIS:.4f} GeV")
print(f"  sensibilidade: d ln m / d ln b_e = (7+N-A mRJ^p)/b_e para o estado R1 = "
      f"{(7 + r1[0] - float(A)*float(mRJ(r1[1]))**float(P))/B_E:.1f}")
