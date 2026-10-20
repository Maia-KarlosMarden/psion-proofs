#!/usr/bin/env python3
"""P14 / WP-SIM: rede de Boltzmann de tuplas de psions (ver P14SIM_DEFINITION.md)."""
import json
import math
from itertools import combinations
import numpy as np
from scipy.integrate import solve_ivp

N, j = 6, 2.5
mk = [k - j for k in range(N)]
Q9 = [1, -5, -5, 10, 4, 4]
TR = list(combinations(range(N), 3))
TI = {t: i for i, t in enumerate(TR)}
NT = len(TR)
conjA = lambda t: tuple(sorted((6 - k) % 6 for k in t))
CJ = [TI[conjA(t)] for t in TR]
LEVELS = (8, 10)

# ---------- operadores internos (como WP-J1) ----------
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
h1 = np.diag([abs(x) for x in mk]).astype(complex)
h1[2, 3] = h1[3, 2] = 0.5
D = 2 ** N


def c_op(k):
    m = np.array([[1.0]])
    for l in range(N):
        f = np.diag([1.0, -1.0]) if l < k else (np.array([[0.0, 1.0], [0.0, 0.0]]) if l == k else np.eye(2))
        m = np.kron(m, f)
    return m.astype(complex)


C = [c_op(k) for k in range(N)]
Cd = [x.conj().T for x in C]
one = lambda M: sum(M[a, b] * Cd[a] @ C[b] for a in range(N) for b in range(N) if abs(M[a, b]) > 0)
H0 = one(h1)
O = sum(one(a) @ one(b) for a, b in zip(A, B))
OJS = 0.5 * (O + O.conj().T)


def vec(t):
    v = np.zeros(D, complex); v[0] = 1.0
    for k in reversed(t):
        v = Cd[k] @ v
    return v


V3 = np.array([vec(t) for t in TR])          # 20 x 64


def Hmat(g):
    H = H0 + g * OJS
    return (V3.conj() @ H @ V3.T)              # 20 x 20 na base de triplas


# ---------- rearranjos 2->2 (mesmo nivel) ----------
REARR = []
for a, b in combinations(range(NT), 2):
    ms = sorted(TR[a] + TR[b])
    for c, d in combinations(range(NT), 2):
        if (c, d) > (a, b) and {c, d} != {a, b} and sorted(TR[c] + TR[d]) == ms:
            REARR.append((a, b, c, d))
REARR = np.array(REARR)

LAM, MPL, GSTAR = 1.2609e12, 1.2209e19, 106.75
G4 = (2 * math.pi ** 2 / 3) ** 2
T_RH, T_END, K_EQ = 1.38e11, 1e8, 10.0
Y0 = 0.208 * 12 / GSTAR


def gamma_over_H(T, cG):
    return cG * G4 * T ** 3 * MPL / (1.66 * math.sqrt(GSTAR) * LAM ** 4)


def build(g, law, symmetrize=False, u=(0.5, 0.5)):
    H = Hmat(g)
    diag = H.diagonal().real
    W = np.abs(H) ** 2
    np.fill_diagonal(W, 0.0)
    if law == "L1":
        kb = diag ** 2
    else:
        kb = np.exp(-(diag - diag.min()))
    kb = kb / kb.mean()
    if symmetrize:
        kb = 0.5 * (kb + kb[CJ])
        W = 0.5 * (W + W[np.ix_(CJ, CJ)])
    return kb, W, np.array(u)


def rhs_factory(kb, W, u, cG):
    idx = np.array(TR)
    def f(x, Y):
        T = T_RH * math.exp(-x)
        r = gamma_over_H(T, cG)
        Yk = Y[:N]
        Yt = Y[N:].reshape(2, NT)
        dk = np.zeros(N)
        dt = np.zeros((2, NT))
        prod = Yk[idx[:, 0]] * Yk[idx[:, 1]] * Yk[idx[:, 2]]
        for li in range(2):
            Jb = kb * u[li] * (prod - Yt[li] / K_EQ)
            dt[li] += Jb
            np.subtract.at(dk, idx[:, 0], Jb); np.subtract.at(dk, idx[:, 1], Jb); np.subtract.at(dk, idx[:, 2], Jb)
            dt[li] += W @ Yt[li] - W.sum(axis=0) * Yt[li]
            a, b, c, d = REARR.T
            Jr = Yt[li, a] * Yt[li, b] - Yt[li, c] * Yt[li, d]
            np.subtract.at(dt[li], a, Jr); np.subtract.at(dt[li], b, Jr)
            np.add.at(dt[li], c, Jr); np.add.at(dt[li], d, Jr)
        return r * np.concatenate([dk, dt.ravel()])
    return f


def run(cG, g, law, symmetrize=False, u=(0.5, 0.5)):
    kb, W, uu = build(g, law, symmetrize, u)
    Y = np.concatenate([np.full(N, 1 / N), np.zeros(2 * NT)])
    xend = math.log(T_RH / T_END)
    sol = solve_ivp(rhs_factory(kb, W, uu, cG), [0, xend], Y, method="LSODA", rtol=1e-10, atol=1e-16)
    Yf = sol.y[:, -1]
    Yt = Yf[N:].reshape(2, NT)
    p, pbar = TI[(0, 4, 5)], TI[(0, 1, 2)]
    dp = Yt[0, p] - Yt[0, pbar]
    de = Yt[1, pbar] - Yt[1, p]
    Npsi0 = 1.0
    Npsi = Yf[:N].sum() + 3 * Yt.sum()
    Qnet = sum(Q9[k] for k in range(N) for _ in [0]) * 0  # placeholder
    Qc = sum(sum(Q9[k] for k in TR[i]) / 9 * (Yt[0, i] + Yt[1, i]) for i in range(NT)) + \
         sum(Q9[k] / 9 * Yf[k] for k in range(N))
    return dict(dp=dp, de=de, Npsi_err=abs(Npsi - Npsi0), lev_diff=np.abs(Yt[0] - Yt[1]).max(),
                Q=Qc, eta=7.04 * Y0 * dp, comp=Yt.sum(), ok=sol.success)


if __name__ == "__main__":
    print("##### Controles #####")
    c1 = run(1e-3, 1.0, "L1", symmetrize=True)
    print(f"  C1 pesos simetrizados: Delta_p = {c1['dp']:.2e}, Delta_e = {c1['de']:.2e}  "
          f"{'OK' if max(abs(c1['dp']), abs(c1['de'])) < 1e-12 else 'FALHA'}")
    res = {}
    print("\n##### Grade (12 corridas) #####")
    print(f"  {'c_G':>6s} {'g':>4s} {'lei':>3s} {'Delta_p':>11s} {'Delta_e':>11s} {'eta(p)':>10s} "
          f"{'Q liquida':>10s} {'N_psi err':>9s} {'|Y8-Y10|':>9s} {'compostos':>9s}")
    for cG in (1e-5, 1e-3, 1.0):
        for g in (0.1, 1.0):
            for law in ("L1", "L2"):
                r = run(cG, g, law)
                res[f"{cG}_{g}_{law}"] = {k: float(v) for k, v in r.items() if k != "ok"}
                print(f"  {cG:6.0e} {g:4.1f} {law:>3s} {r['dp']:+11.3e} {r['de']:+11.3e} {r['eta']:+10.2e} "
                      f"{r['Q']:+10.2e} {r['Npsi_err']:9.1e} {r['lev_diff']:9.1e} {r['comp']:9.3f}")
    allc2 = max(v["Npsi_err"] for v in res.values())
    allc3 = max(v["lev_diff"] for v in res.values())
    print(f"\n  C2 conservacao de N_psi (pior): {allc2:.1e}  {'OK' if allc2 < 1e-8 else 'FALHA'}")
    print(f"  C3 identidade dos niveis (pior): {allc3:.1e}  {'OK' if allc3 < 1e-12 else 'FALHA'}")
    hyd = [k for k, v in res.items() if v["dp"] > 0 and v["de"] > 0]
    print(f"\n  Materia de hidrogenio (Delta_p > 0 e Delta_e > 0): {hyd if hyd else 'NENHUMA corrida'}")
    print(f"  comparacao (so agora): eta observado 6.1e-10")
    json.dump(res, open("P14SIM_results.json", "w"), indent=1)
