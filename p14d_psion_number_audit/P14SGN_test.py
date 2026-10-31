#!/usr/bin/env python3
"""P14 / WP-SGN: fase do vacuo Z6 no termo de enrolamento (ver P14SGN_DEFINITION.md)."""
import io, contextlib, json, math, runpy
import numpy as np
from scipy.optimize import minimize

with contextlib.redirect_stdout(io.StringIO()):
    N_ = runpy.run_path('P14SIM_network.py', run_name='lib')
    L = runpy.run_path('P14LOOP_sim.py', run_name='lib')
one, V3, OJS, J, PL, PR, h1, TR = N_['one'], N_['V3'], N_['OJS'], N_['J'], N_['PL'], N_['PR'], N_['h1'], N_['TR']
NT = len(TR)
Jp = J[0] + 1j * J[1]
Xw = one(Jp @ PL) @ one(Jp @ PR)            # J+ Sigma+ no espaco de Fock
proj = lambda M: V3.conj() @ M @ V3.T
LAM, G, TB = 1 / 6, 1.0, 0.5
PLACE = {"PA": (1, 1), "PB": (1, 2), "PC": (0, 1)}  # (multiplo de theta no bond, multiplo no V)


def Hk(k, place):
    th = 2 * math.pi * k / 6
    mb, mv = PLACE[place]
    hb = h1.copy()
    hb[2, 3] = TB * np.exp(1j * mb * th); hb[3, 2] = np.conj(hb[2, 3])
    V = np.exp(1j * mv * th) * Xw
    V = V + V.conj().T
    return proj(one(hb) + G * OJS + LAM * V)


def tri(H):
    best, arg = 0.0, None
    for a in range(NT):
        for b in range(NT):
            if a == b or abs(H[a, b]) < 1e-14:
                continue
            for c in range(NT):
                if c in (a, b):
                    continue
                v = (H[a, b] * H[b, c] * H[c, a]).imag
                if abs(v) > best:
                    best, arg = abs(v), (a, b, c)
    return best, arg


def removability(H, nphase, seeds=50):
    rng = np.random.default_rng(0)
    if nphase == 6:
        A = np.array([[1.0 if k in t else 0.0 for k in range(6)] for t in TR])
        phase = lambda al: A @ al
    else:
        phase = lambda al: al
    def f(al):
        p = phase(al)
        U = np.exp(1j * p)
        Hr = (U[:, None] * H) * U.conj()[None, :]
        return np.sum(Hr.imag ** 2)
    best = np.inf
    for s in range(seeds):
        x0 = rng.uniform(-math.pi, math.pi, nphase)
        r = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-12, fatol=1e-24))
        r = minimize(f, r.x, method="BFGS", options=dict(gtol=1e-14))
        best = min(best, r.fun)
    return math.sqrt(max(best, 0.0))


def dyn(Hm):
    mk = L['make_rhs']
    mk.__globals__['Hphi'] = lambda phi, g=1.0: Hm
    r, ok = L['run']('L2', 0.0, 10.0, 1.0)
    return r, ok


if __name__ == "__main__":
    IP, IE, IN1, IN2, EXO = L['IP'], L['IE'], L['IN1'], L['IN2'], L['EXO']
    out = {}
    for place in PLACE:
        print(f"\n##### colocacao {place} (bond x{PLACE[place][0]} theta, V x{PLACE[place][1]} theta) #####")
        H1 = Hk(1, place)
        _, arg = tri(H1)
        for k in range(6):
            H = Hk(k, place)
            herm = np.abs(H - H.conj().T).max()
            J3, _ = tri(H)
            J3s = (H[arg[0], arg[1]] * H[arg[1], arg[2]] * H[arg[2], arg[0]]).imag if arg else 0.0
            M2a = removability(H, 6, seeds=20)
            M2b = removability(H, NT, seeds=10)
            r, ok = dyn(H)
            dt = r['dt']
            out[f"{place}_{k}"] = dict(J3=float(J3), J3s=float(J3s), M2_modes=M2a, M2_states=M2b, herm=float(herm),
                                     dp=float(dt[IP]), de=float(dt[IE]), dnu=float(dt[IN1] + dt[IN2]),
                                     exo=float(abs(dt[EXO]).max()), ok=ok)
            print(f"  k={k}: herm {herm:.0e}  J3 = {J3:.3e}  J3s = {J3s:+.3e}  ||Im|| min (6 modos) = {M2a:.2e}  "
                  f"(20 estados, pos-hoc) = {M2b:.2e}  Dp = {dt[IP]:+.3e}  De = {dt[IE]:+.3e}  "
                  f"Dnu = {dt[IN1]+dt[IN2]:+.3e}  max|exo| = {abs(dt[EXO]).max():.2e}  ok={ok}")
    print("\n##### Controles #####")
    c1 = max(max(abs(out[f'{p}_{k}'][q]) for q in ('J3', 'dp', 'de')) for p in PLACE for k in (0, 3))
    print(f"  C1 k = 0, 3: maior |J3|, |Dp|, |De| = {c1:.1e}  {'OK' if c1 < 1e-12 else 'FALHA'}")
    c2 = 0.0
    for p in PLACE:
        for k in (1, 2):
            a, b = out[f'{p}_{k}'], out[f'{p}_{6-k}']
            for q in ('J3s', 'dp', 'de'):
                if abs(a[q]) > 1e-14:
                    c2 = max(c2, abs(b[q] / a[q] + 1))
    print(f"  C2 imparidade k -> 6-k: maior desvio = {c2:.1e}  {'OK' if c2 < 1e-6 else 'FALHA'}")
    c3 = max(v['herm'] for v in out.values())
    print(f"  C3 hermiticidade: {c3:.1e}  {'OK' if c3 < 1e-12 else 'FALHA'}")
    json.dump(out, open('P14SGN_results.json', 'w'), indent=1)
