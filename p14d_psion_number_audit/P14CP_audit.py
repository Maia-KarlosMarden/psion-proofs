#!/usr/bin/env python3
"""P14 / WP-CP: auditoria de CP (ver P14CP_DEFINITION.md)."""
import cmath
import math
from fractions import Fraction as F
import numpy as np

N = 6
j = 2.5
Q9 = [1, -5, -5, 10, 4, 4]


def spin():
    Jz = np.diag([k - j for k in range(N)]).astype(complex)
    Jp = np.zeros((N, N), complex)
    for k in range(N - 1):
        m = k - j
        Jp[k + 1, k] = math.sqrt(j * (j + 1) - m * (m + 1))
    Jm = Jp.conj().T
    return [(Jp + Jm) / 2, (Jp - Jm) / 2j, Jz]


J = spin()
PL = np.diag([1, 1, 1, 0, 0, 0]).astype(complex)
PR = np.diag([0, 0, 0, 1, 1, 1]).astype(complex)

print("##### A1: realidade do tensor interno da eq. (1) #####")
V_S = np.kron(np.eye(N), np.eye(N))
A = [Ja @ PL for Ja in J]
B = [Ja @ PR for Ja in J]
V_D1 = 0.5 * sum(np.kron(a, b) + np.kron(a.conj().T, b.conj().T) for a, b in zip(A, B))


def graded(M, q):
    out = np.zeros_like(M)
    for r in range(N):
        for c in range(N):
            if abs(M[r, c]) > 0 and Q9[r] - Q9[c] == q:
                out[r, c] = M[r, c]
    return out


V_D1p = np.zeros_like(V_D1)
for a, b in zip(A, B):
    for q in range(-30, 31):
        aq, bq = graded(a, q), graded(b, -q)
        V_D1p += 0.5 * (np.kron(aq, bq) + np.kron(aq.conj().T, bq.conj().T))
for name, V in (("(psi-bar psi)^2", V_S), ("J.Sigma D1", V_D1), ("J.Sigma D1'", V_D1p)):
    print(f"  {name:16s}: max|Im| = {np.abs(V.imag).max():.2e} ; hermitiano: {np.allclose(V, V.conj().T)}")
# o que importa e J_y (x) S_y: (iA)(x)(iB) real
print(f"  J_y e imaginaria pura: {np.allclose(J[1].real, 0)} ; J_y (x) J_y real: "
      f"{np.abs(np.kron(J[1], J[1]).imag).max():.1e}")

print("\n##### A2: eps_psi com hopping complexo t = |t| e^{i phi} #####")
t_abs = 0.5
m = [abs(k - j) for k in range(N)]


def eps_psi(phi):
    """Bloco da ligacao |-1/2,L> (k=2) <-> |+1/2,R> (k=3): H2 = [[1/2, t],[t*, 1/2]] + desdobramento
    pela energia de flip j (J.6: tan theta = t/j). eps = sin^2(2 theta), theta do autovetor de
    H = [[0, t],[t*, j]] (desdobramento j = E_flip)."""
    t = t_abs * cmath.exp(1j * phi)
    H = np.array([[0, t], [np.conj(t), j]])
    w, v = np.linalg.eigh(H)
    c = abs(v[0, 0]); s = abs(v[1, 0])
    return (2 * s * c) ** 2


for phi in (0.0, 0.3, math.pi / 3, -math.pi / 3, 1.0, -1.0):
    print(f"  phi = {phi:+.4f}: eps_psi = {eps_psi(phi):.10f}")
print(f"  referencia 25/169 = {25/169:.10f} (tan theta = 1/5)")
print("  eps_psi depende so de |t|: par em phi -> CP-par")

print("\n##### A3: vacuos Z6 do condensado de enrolamento, w_k = e^{2 pi i k/6} #####")
for k in range(N):
    w = cmath.exp(2j * math.pi * k / 6)
    ims = [(w ** n).imag for n in (1, 2, 3)]
    tag = "CP-invariante" if all(abs(x) < 1e-12 for x in ims) else f"conjugado de k={(6-k)%6}"
    print(f"  k={k}: Im(w), Im(w^2), Im(w^3) = {ims[0]:+.4f}, {ims[1]:+.4f}, {ims[2]:+.4f}   {tag}")

print("\n##### A4: vies entre vacuos (K.10): Delta C(s) #####")
mj = [F(2 * k - 5, 2) for k in range(N)]
C0 = sum(x * x for x in mj)
dC = [C0 - sum(mj[k] * mj[(k + s) % N] for k in range(N)) for s in range(N)]
print(f"  Delta C(s), s=0..5 = {[str(x) for x in dC]}  (paper: 0, 15, 24, 27, 24, 15)")
sym = all(dC[s] == dC[(6 - s) % 6] for s in range(N))
print(f"  Delta C(s) = Delta C(6-s): {sym}  -> vacuos CP-conjugados (k, 6-k) degenerados")
for label, f in (("Boltzmann T=dC(1)", lambda E: math.exp(-float(E) / float(dC[1]))),
                 ("Boltzmann T=dC(1)/10", lambda E: math.exp(-10 * float(E) / float(dC[1]))),
                 ("uniforme", lambda E: 1.0),
                 ("so o menor nao trivial", lambda E: 1.0 if E == dC[1] else 0.0)):
    p = [f(dC[k]) for k in range(N)]
    Z = sum(p)
    net = [sum(p[k] / Z * (cmath.exp(2j * math.pi * k / 6) ** n).imag for k in range(N)) for n in (1, 2, 3)]
    print(f"  populacao {label:24s}: sum p_k Im(w_k^n), n=1,2,3 = {net[0]:+.2e}, {net[1]:+.2e}, {net[2]:+.2e}")

print("\n##### A5: kappa_CP #####")
me, LAM = 0.51099895e-3, 1.2609e12
print(f"  (m_e/Lambda_meta)^2 = {(me/LAM)**2:.3e}  (real, positivo: magnitude, nao fase)")

ok_A1 = np.abs(V_D1.imag).max() < 1e-12 and np.abs(V_D1p.imag).max() < 1e-12
ok_A2 = abs(eps_psi(0.7) - eps_psi(-0.7)) < 1e-12
print(f"\n  VEREDITO: A1 real={ok_A1}, A2 par={ok_A2}, A4 simetrico={sym}  ->  "
      f"{'SEM FONTE LIQUIDA DE CP' if (ok_A1 and ok_A2 and sym) else 'HA FONTE (ver acima)'}")
