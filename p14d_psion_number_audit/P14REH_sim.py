#!/usr/bin/env python3
"""P14 / WP-REH: reaquecimento gradual, degraus consistentes, normalizacao (ver P14REH_DEFINITION.md)."""
import io, contextlib, json, math, runpy
import numpy as np
from scipy.integrate import solve_ivp

with contextlib.redirect_stdout(io.StringIO()):
    C = runpy.run_path('P14COMP_sim.py', run_name='lib')
rhs_factory, n, QS, SI, TR, TEQ0 = C['rhs_factory'], C['n'], C['QS'], C['SI'], C['TR'], C['TEQ']
MPL, GS = 1.2209e19, 106.75
ETA = 6.1e-10
th = 2 * math.pi / 6
KAP = 2 * math.pi ** 2 / 3 / (16 * math.pi)
DL = [2 if TR[g] == (0, 1, 2) else (7 if QS[j] == 0 else -1) for j, g in enumerate(SI)]
SW = np.array([math.sin(d * th) for d in DL]); EPS = KAP * np.abs(SW).max(); SW = SW / np.abs(SW).max()
Hrad = lambda T: 1.66 * math.sqrt(GS) * T * T / MPL
Tof = lambda rhoR: (30 * rhoR / (math.pi ** 2 * GS)) ** 0.25


def run(law, TRH, ratio, cG, eps=EPS):
    Gphi = Hrad(TRH)
    HI = ratio * Gphi
    g = rhs_factory(law, eps, SW, math.inf, lambda x: 1.0, lambda x: 1.0)
    Teq = TEQ0 * cG ** (-1 / 3)
    rho0 = 3 * HI ** 2 * MPL ** 2 / (8 * math.pi)
    st = {'Tmax': 0.0, 'Ttr': None}

    def f(N, y):
        lphi, lR = y[0], y[1]
        rphi, rR = math.exp(lphi), math.exp(lR)
        H = math.sqrt(8 * math.pi / (3 * MPL ** 2) * (rphi + rR))
        T = Tof(rR)
        src = (Gphi / H) * rphi / rR
        dphi = -3 - Gphi / H
        dR = -4 + src
        D = 0.75 * src
        Y = y[2:]
        GH = (T / Teq) ** 3 * Hrad(T) / H
        dY = GH * g(0.0, Y) - Y * D
        return np.concatenate([[dphi, dR], dY])

    NRH = math.log(ratio) / 1.5
    Nend = NRH + 16
    y0 = np.concatenate([[math.log(rho0), math.log(rho0 * 1e-30)], np.zeros(4 * n)])
    sol = solve_ivp(f, [0, Nend], y0, method='LSODA', rtol=1e-9, atol=1e-24, dense_output=False, max_step=0.05)
    lR = sol.y[1]; Ts = np.array([Tof(math.exp(v)) for v in lR])
    Tmax = Ts.max()
    # T_RH efetivo: onde H = Gamma_phi
    Hs = np.sqrt(8 * math.pi / (3 * MPL ** 2) * (np.exp(sol.y[0]) + np.exp(lR)))
    i = np.argmin(abs(np.log(Hs / Gphi)))
    Y = sol.y[2:, -1]
    tB, tL, aB, aL = Y[:n], Y[n:2 * n], Y[2 * n:3 * n], Y[3 * n:]
    B, L = (tB - aB).sum(), (tL - aL).sum()
    Q = (QS * (tB - aB + tL - aL)).sum()
    return dict(eta=7.04 * B, B=B, L=L, Q=Q, Tmax=Tmax, TRHeff=Ts[i], ok=sol.success,
                phi_end=math.exp(sol.y[0, -1] - lR[-1]), Tend=Ts[-1])


if __name__ == "__main__":
    print(f"Delta-l (singletos {[TR[g] for g in SI]}): {DL}; pesos {np.round(SW, 3).tolist()}; eps = {EPS:.4f}")
    r0 = run('L1', 1e10, 1e4, 1.0, eps=0.0)
    print(f"C1 eps = 0: |eta| = {abs(r0['eta']):.1e}  {'OK' if abs(r0['eta']) < 1e-15 else 'FALHA'}")
    grid = (1e8, 3e8, 1e9, 3e9, 1e10, 2e10, 4e10, 1e11)
    res, c2, c3, inband, tot = {}, 0.0, 0.0, 0, 0
    for cG in (0.1, 1.0, 10.0):
        for ratio in (1e2, 1e4, 1e6, 1e8):
            for law in ('L1', 'L2'):
                row = []
                for TRH in grid:
                    r = run(law, TRH, ratio, cG)
                    tot += 1
                    if 1e-10 <= abs(r['eta']) <= 1e-9:
                        inband += 1
                    if abs(r['B']) > 1e-12:
                        c2 = max(c2, abs(r['B'] + r['L']) / abs(r['B']), abs(r['Q']) / abs(r['B']))
                    c3 = max(c3, abs(r['TRHeff'] / TRH - 1))
                    row.append((TRH, r['eta'], r['Tmax']))
                    res[f"{cG}_{ratio}_{law}_{TRH}"] = {k: float(v) for k, v in r.items() if k != 'ok'}
                roots = []
                for (T1, e1, _), (T2, e2, _) in zip(row, row[1:]):
                    a, b = abs(e1), abs(e2)
                    if a > 0 and b > 0 and (a - ETA) * (b - ETA) < 0:
                        roots.append(math.exp(math.log(T1) + (math.log(ETA) - math.log(a)) * math.log(T2 / T1) / math.log(b / a)))
                print(f"  c_G={cG:4.1f} H_I/G={ratio:.0e} {law}: Tmax/T_RH~{row[3][2]/3e9:.0f}  eta: "
                      + ", ".join(f"{T:.0e}:{e:+.1e}" for T, e, _ in row)
                      + f"\n      raizes T_RH* = {[f'{x:.1e}' for x in roots]}")
    print(f"\nC2 max(|B+L|/|B|, |Q|/|B|) = {c2:.1e}  {'OK' if c2 < 1e-9 else 'FALHA'}")
    print(f"C3 desvio do T_RH efetivo = {c3:.2f}  {'OK' if c3 < 0.3 else 'FALHA'}")
    print(f"fracao da grade com 1e-10 <= |eta| <= 1e-9: {inband}/{tot}")
    json.dump(res, open('P14REH_results.json', 'w'), indent=1)
