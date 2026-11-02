#!/usr/bin/env python3
"""P14 / WP-COMP: bifurcacao por criacao de pares de compositos (ver P14COMP_DEFINITION.md)."""
import io, contextlib, json, math, runpy
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp

with contextlib.redirect_stdout(io.StringIO()):
    S = runpy.run_path('P14NOW_sim.py', run_name='lib')
TR, REARR, weights = S['TR'], S['REARR'], S['weights']
SI = [i for i, t in enumerate(TR) if sum(t) % 6 in (0, 3)]
n = len(SI)
pos = {g: j for j, g in enumerate(SI)}
RE = np.array([[pos[a], pos[b], pos[c], pos[d]] for a, b, c, d in REARR if all(x in pos for x in (a, b, c, d))])
Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]
QS = np.array([float(sum(Q9[k] for k in TR[g])) for g in SI])
YEQ = 0.208 * 2 / 106.75
TEQ = (1.66 * math.sqrt(106.75) * 1.2609e12 ** 4 / ((2 * math.pi ** 2 / 3) ** 2 * 1.2209e19)) ** (1 / 3)
be = 325 / 1152; LAM = 1.2609e12
EB = {"inf": math.inf, "m1": LAM * math.exp(-1 / be), "m2": LAM * math.exp(-2 / be)}
ETA = 6.1e-10
th = 2 * math.pi / 6
SHW = np.array([math.sin((7 if QS[j] == 0 else 2) * th) for j in range(n)]) / math.sin(2 * th)
ra, rb, rc, rd = RE.T if len(RE) else (np.array([], int),) * 4


def rhs_factory(law, eps, s, Eb, rate, Tof):
    kb, W = weights(1.0, law, 0.0, 1)
    k = kb[SI]; k = k / k.mean()
    Ws = W[np.ix_(SI, SI)]

    def trans(T):
        return Ws @ T - Ws.sum(axis=1) * T

    def rear(T):
        d = np.zeros(n)
        if len(RE):
            Jr = T[ra] * T[rb] - T[rc] * T[rd]
            np.subtract.at(d, ra, Jr); np.subtract.at(d, rb, Jr)
            np.add.at(d, rc, Jr); np.add.at(d, rd, Jr)
        return d / YEQ

    def f(x, Y):
        r = rate(x); T = Tof(x)
        tB, tL, aB, aL = Y[:n], Y[n:2 * n], Y[2 * n:3 * n], Y[3 * n:]
        J = lambda u, v: k * (YEQ ** 2 - u * v) / YEQ
        JBB, JLL, JBL, JLB = J(tB, aB), J(tL, aL), J(tB, aL), J(tL, aB)
        d = 0.5 * eps * s * (JBB + JLL + JBL + JLB)
        dtB = JBB + JBL + d + trans(tB) + rear(tB)
        daL = JLL + JBL + d + trans(aL) + rear(aL)
        dtL = JLL + JLB - d + trans(tL) + rear(tL)
        daB = JBB + JLB - d + trans(aB) + rear(aB)
        if Eb != math.inf:
            w = math.exp(-Eb / T)
            dtB -= 0.5 * w * (tB - tL); dtL += 0.5 * w * (tB - tL)
            daB -= 0.5 * w * (aB - aL); daL += 0.5 * w * (aB - aL)
        return r * np.concatenate([dtB, dtL, daB, daL])
    return f


def run(law, eps, s, Eb, TRH):
    Tof = lambda x: TRH * math.exp(-x)
    f = rhs_factory(law, eps, s, Eb, lambda x: (Tof(x) / TEQ) ** 3, Tof)
    sol = solve_ivp(f, [0, math.log(1e3)], np.zeros(4 * n), method='LSODA', rtol=1e-10, atol=1e-26)
    Y = sol.y[:, -1]
    tB, tL, aB, aL = Y[:n], Y[n:2 * n], Y[2 * n:3 * n], Y[3 * n:]
    B, L = (tB - aB).sum(), (tL - aL).sum()
    Q = (QS * (tB - aB + tL - aL)).sum()
    return dict(B=B, L=L, Q=Q, eta=7.04 * B, ok=sol.success)


if __name__ == "__main__":
    print(f"singletos: {[TR[g] for g in SI]}; rearranjos entre singletos: {len(RE)}; T_eq = {TEQ:.3e}; "
          f"E_b(m1) = {EB['m1']:.2e}, E_b(m2) = {EB['m2']:.2e}")
    print("##### Controles #####")
    c1 = max(abs(run(l, 0.0, np.ones(n), math.inf, 1.4e10)['eta']) for l in ("L1", "L2"))
    print(f"  C1 eps = 0: |eta| = {c1:.1e}  {'OK' if c1 < 1e-18 else 'FALHA'}")
    f = rhs_factory("L1", 0.1, np.ones(n), math.inf, lambda x: 1e3, lambda x: 1e10)
    sol = solve_ivp(f, [0, 20], np.full(4 * n, YEQ), method='LSODA', rtol=1e-10, atol=1e-26)
    Y = sol.y[:, -1]; c4 = abs((Y[:n] - Y[2 * n:3 * n]).sum())
    print(f"  C4 equilibrio (eps = 0,1): |B| = {c4:.1e}  {'OK' if c4 < 1e-16 else 'FALHA'}")
    grid = (1.38e11, 4e10, 2e10, 1.4e10, 7e9, 4e9, 1.4e9, 7e8, 3e8, 1.5e8)
    res, c2, c3 = {}, 0.0, 0.0
    kap = 2 * math.pi ** 2 / 3 / (16 * math.pi)
    for bias, s, eps in (("U", np.ones(n), 1e-3), ("HW", SHW, kap * math.sin(2 * th))):
        print(f"\n##### vies {bias} (eps = {eps:.3e}) #####")
        for law in ("L1", "L2"):
            for ebk, Eb in EB.items():
                row = []
                for TRH in grid:
                    r = run(law, eps, s, Eb, TRH)
                    if abs(r['B']) > 1e-20:
                        c2 = max(c2, abs(r['B'] + r['L']) / abs(r['B']))
                    c2 = max(c2, abs(r['Q']) / max(abs(r['B']), 1e-30) if abs(r['B']) > 1e-20 else 0)
                    row.append((TRH, r['eta']))
                    res[f"{bias}_{law}_{ebk}_{TRH}"] = float(r['eta'])
                row.sort()
                cross = []
                for (T1, e1), (T2, e2) in zip(row, row[1:]):
                    a, b = abs(e1), abs(e2)
                    if a > 0 and b > 0 and (a - ETA) * (b - ETA) < 0:
                        cross.append(math.exp(math.log(T1) + (math.log(ETA) - math.log(a)) * math.log(T2 / T1) / math.log(b / a)))
                tag = "eta/eps" if bias == "U" else "eta"
                fac = 1 / eps if bias == "U" else 1.0
                print(f"  {law} E_b={ebk:3s}: {tag} = " + ", ".join(f"{T:.1e}:{e*fac:+.1e}" for T, e in row)
                      + (f"\n      T_RH* = {[f'{c:.2e}' for c in cross]}" if bias == "HW" else ""))
        rp = run("L1", eps, s, math.inf, 1.4e10); rn = run("L1", -eps, s, math.inf, 1.4e10)
        c3 = max(c3, abs(rn['eta'] / rp['eta'] + 1))
    print(f"\n  C2 max(|B+L|/|B|, |Q|/|B|) = {c2:.1e}  {'OK' if c2 < 1e-12 else 'FALHA'}")
    print(f"  C3 eps -> -eps: desvio {c3:.1e}  {'OK' if c3 < 1e-6 else 'FALHA'}")
    json.dump(res, open('P14COMP_results.json', 'w'), indent=1)
