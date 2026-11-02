#!/usr/bin/env python3
"""P14 / WP-RUNG: selecao por proto-cor e degrau como enrolamento (ver P14RUNG_DEFINITION.md)."""
import io, contextlib, json, math, runpy
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp

with contextlib.redirect_stdout(io.StringIO()):
    B = runpy.run_path('P14BIF_sim.py', run_name='lib')
TR, NT = B['TR'], B['NT']
mk, summ, goH, Y0 = B['make_rhs'], B['summarize'], B['gamma_over_H'], B['Y0']
Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]
NAMES = {(0, 1, 2): "e- / Delta- (B)", (0, 4, 5): "p", (0, 1, 5): "n / nu_tau", (0, 2, 4): "Sigma0 / nu_mu",
         (1, 2, 3): "Lambda0 / nu_e", (1, 3, 5): "Sigma+-like", (2, 3, 4): "Sigma+-like", (3, 4, 5): "Delta++-like (Q=+2)"}
ETA = 6.1e-10


def Dval(t):
    M = sum(F(2 * k - 5, 2) for k in t); Q = sum(Q9[k] for k in t)
    return (2 * M - 6 * Q + 3) % 4


singlet = [t for t in TR if sorted(k % 3 for k in t) == [0, 1, 2]]
Qt = {t: sum(Q9[k] for k in t) for t in TR}


def weights(k):
    th = 2 * math.pi * k / 6
    s = np.zeros(NT)
    for i, t in enumerate(TR):
        if t in singlet:
            dl = 7 if Qt[t] == 0 else 2
            s[i] = math.sin(dl * th)
    m = np.abs(s).max()
    return (s / m if m > 0 else s), m


def eta_run(law, k, kappa, TRH):
    s, m = weights(k)
    key = f"W{k}"
    mk.__globals__['SW'][key] = s
    eps = kappa * m
    f = mk(law, eps, key, lambda x: goH(TRH * math.exp(-x), 1.0))
    y0 = np.concatenate([np.full(12, 1 / 12), np.zeros(80)])
    sol = solve_ivp(f, [0, math.log(1e3)], y0, method='LSODA', rtol=1e-11, atol=1e-24)
    r = summ(sol.y[:, -1])
    return 7.04 * Y0 * r['B'], r


if __name__ == "__main__":
    print("##### Parte 1: tripas singleto de proto-cor #####")
    print(f"  singletos ({len(singlet)}):")
    for t in singlet:
        print(f"    {t}  Q = {str(Qt[t]):>3s}  D = {Dval(t)}  {NAMES.get(t, '?')}")
    rest = [t for t in TR if t not in singlet]
    frac = [t for t in rest if Qt[t].denominator != 1]
    neu = [t for t in rest if Qt[t] == 0]
    print(f"  nao singletos: {len(rest)}; de carga fracionaria: {len(frac)}; inteiras: {[ (t, str(Qt[t]), Dval(t)) for t in rest if Qt[t].denominator == 1]}")
    print(f"  as duas neutras com D = 2 estao entre os nao singletos: {all(t in rest for t in [(0,1,4),(0,2,5)])}")

    print("\n##### Parte 3: eta com vies de enrolamento (H-W) #####")
    for k in range(6):
        s, m = weights(k)
        print(f"  k={k}: sin(2 theta) = {math.sin(2*2*math.pi*k/6):+.3f}, sin(7 theta) = {math.sin(7*2*math.pi*k/6):+.3f}")
    grid = (1.38e10, 4e9, 1.4e9, 7e8, 3e8, 1.5e8, 7e7)
    extra = (2e10, 3e10, 4e10, 5e10)
    res = {}
    c1 = c2 = c3 = 0.0
    for law in ("L1", "L2"):
        for kappa, lab in ((1 / 6, "alpha_meta/pi"), (2 * math.pi ** 2 / 3 / (16 * math.pi), "g2/16pi"), (1 / 36, "lambda^2")):
            for k in (1, 2):
                row = []
                for TRH in grid + extra:
                    e, r = eta_run(law, k, kappa, TRH)
                    row.append((TRH, e))
                    res[f"{law}_{lab}_{k}_{TRH}"] = float(e)
                    if abs(r['B']) > 1e-12:
                        c3 = max(c3, abs(r['B'] + r['L']) / abs(r['B']), abs(r['Q']), abs(r['Npsi']))
                    if TRH in (1.38e10, 7e8):
                        e2, _ = eta_run(law, 6 - k, kappa, TRH)
                        if abs(e) > 1e-15:
                            c2 = max(c2, abs(e2 / e + 1))
                row.sort()
                cross = []
                for (T1, e1), (T2, e2) in zip(row, row[1:]):
                    for a, b in ((abs(e1), abs(e2)),):
                        if (a - ETA) * (b - ETA) < 0 and a > 0 and b > 0:
                            lt = math.log(T1) + (math.log(ETA) - math.log(a)) * (math.log(T2) - math.log(T1)) / (math.log(b) - math.log(a))
                            cross.append(math.exp(lt))
                sgn = "+" if row[len(row)//2][1] > 0 else "-"
                print(f"  {law} kappa={lab:13s} k={k}: eta(T_RH) = " + ", ".join(f"{T:.1e}:{e:+.1e}" for T, e in row)
                      + f"\n      sinal de eta: {sgn}; T_RH* (eta = 6,1e-10) = {[f'{c:.2e}' for c in cross]}")
    for law in ("L1",):
        for k in (0, 3):
            e, _ = eta_run(law, k, 1 / 6, 1.38e10)
            c1 = max(c1, abs(e))
    print(f"\n  C1 k = 0, 3: |eta| = {c1:.1e}  {'OK' if c1 < 1e-15 else 'FALHA'}")
    print(f"  C2 k -> 6-k inverte eta: desvio {c2:.1e}  {'OK' if c2 < 1e-6 else 'FALHA'}")
    print(f"  C3 max(|B+L|/|B|, |Q|, |N_psi|) = {c3:.1e}  {'OK' if c3 < 1e-8 else 'FALHA'}")
    print("  (pontos T_RH = 2e10-5e10 sao extras declarados pos-registro, para a raiz do lado da lavagem)")
    json.dump(res, open('P14RUNG_results.json', 'w'), indent=1)
