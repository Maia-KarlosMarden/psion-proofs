#!/usr/bin/env python3
"""P14 / WP-J1: assimetria materia/antimateria de J.1 nas conjugacoes C_A e C_B."""
import math
from itertools import combinations
from fractions import Fraction as F
import numpy as np

N, j = 6, 2.5
Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]
mk = [k - j for k in range(N)]
NAMES = {(0, 1, 2): "e-", (0, 4, 5): "p", (0, 1, 5): "n/nu_tau", (0, 2, 4): "nu_mu/S0",
         (1, 2, 3): "nu_e/L0", (1, 3, 5): "S+", (2, 3, 4): "S+'", (3, 4, 5): "exotico Q=+2"}

CONJ = {"C_A: k->6-k": [(6 - k) % 6 for k in range(N)],
        "C_B: k->k+3 (Z2)": [(k + 3) % 6 for k in range(N)]}


def spin():
    Jz = np.diag(mk).astype(complex)
    Jp = np.zeros((N, N), complex)
    for k in range(N - 1):
        Jp[k + 1, k] = math.sqrt(j * (j + 1) - mk[k] * (mk[k] + 1))
    Jm = Jp.conj().T
    return [(Jp + Jm) / 2, (Jp - Jm) / 2j, Jz]


J = spin()
PL = np.diag([1, 1, 1, 0, 0, 0]).astype(complex)
PR = np.eye(N) - PL
A = [Ja @ PL for Ja in J]
B = [Ja @ PR for Ja in J]
h = np.diag([abs(x) for x in mk]).astype(complex)
h[2, 3] = h[3, 2] = 0.5


def graded(M, q):
    out = np.zeros_like(M)
    for r in range(N):
        for c in range(N):
            if abs(M[r, c]) > 0 and 9 * (Q9[r] - Q9[c]) == q:
                out[r, c] = M[r, c]
    return out


PAIRS = {"D1": list(zip(A, B)),
         "D1'": [(graded(a, q), graded(b, -q)) for a, b in zip(A, B) for q in range(-30, 31)
                 if np.abs(graded(a, q)).max() > 0 and np.abs(graded(b, -q)).max() > 0]}

# ---------- Fock space (Jordan-Wigner) ----------
D = 2 ** N
def c_op(k):
    m = np.array([[1.0]])
    for l in range(N):
        if l < k:
            f = np.diag([1.0, -1.0])
        elif l == k:
            f = np.array([[0.0, 1.0], [0.0, 0.0]])
        else:
            f = np.eye(2)
        m = np.kron(m, f)
    return m.astype(complex)


C = [c_op(k) for k in range(N)]
Cd = [x.conj().T for x in C]


def one_body(M):
    return sum(M[a, b] * Cd[a] @ C[b] for a in range(N) for b in range(N) if abs(M[a, b]) > 0)


def basis_index(t):
    # estado c†_{t0} c†_{t1} c†_{t2} |0>, |0> = todos desocupados = indice 0 com convencao abaixo
    v = np.zeros(D, complex)
    v[0] = 1.0
    for k in reversed(t):
        v = Cd[k] @ v
    return v


TRIPLES = list(combinations(range(N), 3))
VEC = {t: basis_index(t) for t in TRIPLES}
Nop = sum(Cd[k] @ C[k] for k in range(N))
idx3 = [i for i in range(D) if abs(Nop[i, i].real - 3) < 1e-9]

H0 = one_body(h)
OJS = {}
for name, pairs in PAIRS.items():
    O = sum(one_body(a) @ one_body(b) for a, b in pairs)
    OJS[name] = 0.5 * (O + O.conj().T)


def perm_matrix(p):
    P = np.zeros((N, N))
    for k in range(N):
        P[p[k], k] = 1
    return P


def norm(X):
    return np.linalg.norm(X)


print("##### Passo 2: paridade sob as conjugacoes #####")
for cname, p in CONJ.items():
    P = perm_matrix(p)
    hp = P @ h @ P.T
    print(f"  {cname}")
    print(f"    h (um psion)          : |impar|/|total| = {norm((h - hp)/2)/norm(h):.4f}")
    for vname, pairs in PAIRS.items():
        V = sum(np.kron(a, b) for a, b in pairs)
        V = 0.5 * (V + V.conj().T)
        PP = np.kron(P, P)
        Vp = PP @ V @ PP.T
        print(f"    J.Sigma {vname:3s} (tensor) : |par|/|total| = {norm((V + Vp)/2)/norm(V):.4f} ; "
              f"|impar|/|total| = {norm((V - Vp)/2)/norm(V):.4f}")

print("\n##### Passo 3: energia da tupla vs conjugada (setor de 3 psions) #####")
for cname, p in CONJ.items():
    print(f"\n  === {cname} ===")
    for vname in ("D1", "D1'"):
        for g in (0.1, 1.0):
            H = H0 + g * OJS[vname]
            Hs = H[np.ix_(idx3, idx3)]
            w, U = np.linalg.eigh(Hs)
            print(f"   vertice {vname}, g = {g}:")
            print(f"     {'tupla':>10s} {'nome':>13s} {'Q':>4s} -> {'C(t)':>10s} {'Q':>4s}   "
                  f"{'dE0':>7s} {'dE1':>7s} {'dE_exato':>9s}  favorecida")
            for t in [(0, 1, 2), (0, 4, 5), (0, 1, 5), (0, 2, 4), (1, 2, 3), (1, 3, 5), (2, 3, 4), (3, 4, 5)]:
                ct = tuple(sorted(p[k] for k in t))
                vt, vc = VEC[t], VEC[ct]
                e0 = (vt.conj() @ H0 @ vt).real - (vc.conj() @ H0 @ vc).real
                e1 = (vt.conj() @ H @ vt).real - (vc.conj() @ H @ vc).real
                # autoestado exato de maior sobreposicao
                def ex(v):
                    ov = np.abs(U.conj().T @ v[idx3]) ** 2
                    return w[int(np.argmax(ov))]
                ee = ex(vt) - ex(vc)
                qt, qc = sum(Q9[k] for k in t), sum(Q9[k] for k in ct)
                fav = "propria" if ee < -1e-9 else ("conjugada" if ee > 1e-9 else "empate")
                same = " (auto-conjugada)" if ct == t else ""
                print(f"     {str(t):>10s} {NAMES[t]:>13s} {str(qt):>4s} -> {str(ct):>10s} {str(qc):>4s}   "
                      f"{e0:+7.3f} {e1:+7.3f} {ee:+9.4f}  {fav}{same}")

print("\n##### Criterio (iii): conjugacao de carga Q(C(t)) = -Q(t) #####")
for cname, p in CONJ.items():
    ok = []
    for t in [(0, 1, 2), (0, 4, 5)]:
        ct = tuple(sorted(p[k] for k in t))
        ok.append(sum(Q9[k] for k in ct) == -sum(Q9[k] for k in t))
    print(f"  {cname}: eletron {ok[0]}, proton {ok[1]}")
