#!/usr/bin/env python3
"""P14 / WP-TRH: bifurcacao de nivel contra a temperatura de reaquecimento (ver P14TRH_DEFINITION.md)."""
import io, contextlib, json, math, runpy
import numpy as np
from scipy.integrate import solve_ivp
with contextlib.redirect_stdout(io.StringIO()):
    B = runpy.run_path('P14BIF_sim.py', run_name='lib')
mk, summ, goH, Y0, TEND = B['make_rhs'], B['summarize'], B['gamma_over_H'], B['Y0'], B['T_END']
EPS, ETA = 1e-3, 6.1e-10

def run(law, TRH, cG):
    f = mk(law, EPS, 'S1', lambda x: goH(TRH * math.exp(-x), cG))
    y0 = np.concatenate([np.full(12, 1 / 12), np.zeros(80)])
    sol = solve_ivp(f, [0, math.log(TRH / TEND)], y0, method='LSODA', rtol=1e-11, atol=1e-22)
    return summ(sol.y[:, -1])

if __name__ == "__main__":
    Teq = (1.66 * math.sqrt(106.75) * 1.2609e12 ** 4 / ((2 * math.pi ** 2 / 3) ** 2 * 1.2209e19)) ** (1 / 3)
    print(f"T_eq (Gamma/H = 1, c_G = 1) = {Teq:.3e} GeV; Lambda_meta = 1.2609e12 GeV")
    a = run('L1', 1.38e11, 1e-3); b = run('L1', 1.38e10, 1.0)
    print(f"equivalencia (c_G=1e-3, T_RH=1.38e11) vs (c_G=1, T_RH=1.38e10): B = {a['B']:.6e} vs {b['B']:.6e}, "
          f"diferenca relativa {abs(a['B']/b['B']-1):.1e}")
    res, c1, c2, c3 = {}, 0.0, 0.0, 0.0
    for law in ('L1', 'L2'):
        print(f"\n  lei {law}:")
        for TRH in (1.38e11, 7e10, 4e10, 2e10, 1.38e10, 7e9, 4e9, 2e9, 1.4e9, 7e8):
            r = run(law, TRH, 1.0)
            e = 7.04 * Y0 * r['B'] / EPS
            need = ETA / abs(e) if e else float('inf')
            if abs(r['B']) > 1e-12:
                c1 = max(c1, abs(r['B'] + r['L']) / abs(r['B']))
            c2 = max(c2, abs(r['Q'])); c3 = max(c3, abs(r['Npsi']))
            res[f"{law}_{TRH}"] = dict(eta_per_eps=float(e), eps_needed=float(need), GH=float(goH(TRH, 1.0)))
            print(f"    T_RH = {TRH:.2e}  Gamma/H = {goH(TRH,1.0):8.1e}  eta/eps = {e:+.3e}  eps necessario = {need:.2e}")
        sub = {k: v for k, v in res.items() if k.startswith(law)}
        kmax = max(sub, key=lambda k: abs(sub[k]['eta_per_eps']))
        ok = [float(k.split('_')[1]) for k, v in sub.items() if v['eps_needed'] <= 1e-5]
        print(f"    maximo: T_RH = {float(kmax.split('_')[1]):.2e}, eta/eps = {sub[kmax]['eta_per_eps']:+.2e}; "
              f"eps necessario <= 1e-5 para T_RH em [{min(ok):.1e}, {max(ok):.1e}] GeV")
    print(f"\n  C1 |B+L|/|B| = {c1:.1e} {'OK' if c1 < 1e-9 else 'FALHA'};  C2 |Q| = {c2:.1e} {'OK' if c2 < 1e-12 else 'FALHA'};"
          f"  C3 |N_psi| = {c3:.1e} {'OK' if c3 < 1e-12 else 'FALHA'}")
    json.dump(res, open('P14TRH_results.json', 'w'), indent=1)
