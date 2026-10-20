#!/usr/bin/env python3
"""P14 / WP-BIF: bifurcacao de nivel na formacao (ver P14BIF_DEFINITION.md)."""
import io, contextlib, json, math, runpy
import numpy as np
from scipy.integrate import solve_ivp

with contextlib.redirect_stdout(io.StringIO()):
    S = runpy.run_path('P14NOW_sim.py', run_name='lib')
TR, NT, REARR, weights = S['TR'], S['NT'], S['REARR'], S['weights']
gamma_over_H, T_RH, T_END, K_EQ = S['gamma_over_H'], S['T_RH'], S['T_END'], S['K_EQ']
idx = np.array(TR)
ra, rb, rc, rd = REARR.T
YEQ, YEQT = 1 / 12, K_EQ * (1 / 12) ** 3
Q9 = np.array([1, -5, -5, 10, 4, 4]) / 9
QT = np.array([Q9[list(t)].sum() for t in TR])
SW = {"S1": np.ones(NT), "S2": np.array([(sum(k <= 2 for k in t) - sum(k >= 3 for k in t)) / 3 for t in TR])}
Y0 = 0.208 * 12 / 106.75
ETA = 6.1e-10
RTOL = 1e-11


def sector_ops(W):
    def trans(T):
        return W @ T - W.sum(axis=1) * T

    def rear(T):
        Jr = T[ra] * T[rb] - T[rc] * T[rd]
        d = np.zeros(NT)
        np.subtract.at(d, ra, Jr); np.subtract.at(d, rb, Jr)
        np.add.at(d, rc, Jr); np.add.at(d, rd, Jr)
        return d
    return trans, rear


def make_rhs(law, eps, sw, rate):
    kb, W = weights(1.0, law, 0.0, 1)
    trans, rear = sector_ops(W)
    s = SW[sw]

    def f(x, Y):
        r = rate(x)
        psi, apsi = Y[:6], Y[6:12]
        tB, tL, aB, aL = (Y[12 + 20 * i:32 + 20 * i] for i in range(4))
        out = []
        dps = []
        for P, B_, L_, sg in ((psi, tB, tL, +1), (apsi, aB, aL, -1)):
            P3 = P[idx[:, 0]] * P[idx[:, 1]] * P[idx[:, 2]]
            J = kb * (P3 - (B_ + L_) / K_EQ)
            dP = np.zeros(6)
            for col in range(3):
                np.subtract.at(dP, idx[:, col], J)
            src = 0.5 * sg * eps * s * J
            dB = 0.5 * kb * P3 - kb * B_ / K_EQ + src + trans(B_) + rear(B_)
            dL = 0.5 * kb * P3 - kb * L_ / K_EQ - src + trans(L_) + rear(L_)
            out.append((dB, dL)); dps.append(dP)
        (dtB, dtL), (daB, daL) = out
        JB = tB * aB - (YEQT / 2) ** 2
        JL = tL * aL - (YEQT / 2) ** 2
        dtB -= JB; daB -= JB; dtL -= JL; daL -= JL
        Ja = psi * apsi - YEQ ** 2
        return r * np.concatenate([dps[0] - Ja, dps[1] - Ja, dtB, dtL, daB, daL])
    return f


def summarize(Y):
    psi, apsi = Y[:6], Y[6:12]
    tB, tL, aB, aL = (Y[12 + 20 * i:32 + 20 * i] for i in range(4))
    Bv, Lv = tB - aB, tL - aL
    B, L = Bv.sum(), Lv.sum()
    dfree = (psi - apsi).sum()
    return dict(B=B, L=L, dfree=dfree, Npsi=dfree + 3 * (B + L), Q=(QT * (Bv + Lv)).sum(),
                Bv=Bv, Lv=Lv)


def run(law, eps, sw, cG):
    f = make_rhs(law, eps, sw, lambda x: gamma_over_H(T_RH * math.exp(-x), cG))
    Y0v = np.concatenate([np.full(12, 1 / 12), np.zeros(80)])
    sol = solve_ivp(f, [0, math.log(T_RH / T_END)], Y0v, method="LSODA", rtol=RTOL, atol=1e-22)
    return summarize(sol.y[:, -1]), sol.success


if __name__ == "__main__":
    print("##### Controles #####")
    c1 = max(max(abs(run(l, 0.0, sw, 1e-3)[0][k]) for k in ("B", "L", "Q")) for l in ("L1", "L2") for sw in ("S1", "S2"))
    print(f"  C1 eps = 0: maior |B|, |L|, |Q| = {c1:.1e}  {'OK' if c1 < 1e-12 else 'FALHA'}")
    c5 = 0.0
    for l in ("L1", "L2"):
        f = make_rhs(l, 0.1, "S1", lambda x: 1e3)
        Yeq = np.concatenate([np.full(12, YEQ), np.full(80, YEQT / 2)])
        sol = solve_ivp(f, [0, 20], Yeq, method="LSODA", rtol=RTOL, atol=1e-22)
        r = summarize(sol.y[:, -1])
        c5 = max(c5, abs(r['B']), abs(r['L']))
    print(f"  C5 equilibrio com eps = 0,1: maior |B|, |L| = {c5:.1e}  {'OK' if c5 < 1e-12 else 'FALHA'}")

    print("\n##### Grade, eps = 1e-3 #####")
    eps = 1e-3
    res, c2, c3, c4 = {}, 0.0, 0.0, 0.0
    best = None
    for sw in ("S1", "S2"):
        for law in ("L1", "L2"):
            for cG in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0):
                r, ok = run(law, eps, sw, cG)
                c2 = max(c2, abs(r['Npsi']))
                if abs(r['B']) > 0:
                    c3 = max(c3, abs(r['B'] + r['L']) / abs(r['B']))
                c4 = max(c4, abs(r['Q']))
                eta_per_eps = 7.04 * Y0 * r['B'] / eps
                res[f"{sw}_{law}_{cG}"] = dict(B=float(r['B']), L=float(r['L']), dfree=float(r['dfree']),
                                               eta_per_eps=float(eta_per_eps), ok=ok)
                print(f"  {sw} {law} cG={cG:6.0e}: B/eps = {r['B']/eps:+.3e}  L/eps = {r['L']/eps:+.3e}  "
                      f"(B+L)/|B| = {(r['B']+r['L'])/abs(r['B']) if r['B'] else 0:+.1e}  Q = {r['Q']:+.1e}  "
                      f"eta/eps = {eta_per_eps:+.3e}  ok={ok}")
                if best is None or abs(eta_per_eps) > abs(best[0]):
                    best = (eta_per_eps, sw, law, cG)
    print(f"\n  C2 |N_psi| maior = {c2:.1e}  {'OK' if c2 < 1e-12 else 'FALHA'}")
    print(f"  C3 (B+L)/|B| maior = {c3:.1e}  {'OK' if c3 < 1e-12 else 'FALHA'}")
    print(f"  C4 |Q| maior = {c4:.1e}  {'OK' if c4 < 1e-12 else 'FALHA'}")

    e1, sw, law, cG = best
    rp = run(law, eps, sw, cG)[0]; rn = run(law, -eps, sw, cG)[0]; r10 = run(law, 10 * eps, sw, cG)[0]
    c6 = abs(rn['B'] / rp['B'] + 1)
    lin = abs(r10['B'] / (10 * rp['B']) - 1)
    print(f"  C6 eps -> -eps: desvio {c6:.1e}  {'OK' if c6 < 1e-6 else 'FALHA'}; linearidade (1e-2 vs 1e-3) {lin:.1e}")

    print(f"\n##### Engenharia reversa (melhor ponto: {sw} {law} cG={cG:.0e}, eta/eps = {e1:+.3e}) #####")
    for lab, D in (("D = 1", 1.0), ("xi* = 3,6 (D = 2,1e-2)", 3.6 ** -3), ("xi* = 39 (D = 1,7e-5)", 39.0 ** -3)):
        need = ETA / (abs(e1) * D)
        print(f"  {lab:24s}: eps necessario = {need:.2e}  {'(< 1, possivel)' if need < 1 else '(> 1, impossivel)'}"
              f"   com eps = eps_psi = 0,148: eta = {abs(e1)*0.148*D:.2e}")

    print("\n##### Contabilidade final (B, L, Q conservados; estados mais leves p, e, nu) #####")
    B, L = rp['B'], rp['L']
    Np = B; Ne = Np; nub = -L - Ne
    print(f"  por barion: N_e/N_p = {Ne/Np:.6f}, (N_nubar - N_nu)/N_p = {nub/Np:.6f}  (estrutura prevista: 1 e 2)")
    json.dump(res, open("P14BIF_results.json", "w"), indent=1)
