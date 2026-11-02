#!/usr/bin/env python3
"""P14 / WP-C2 (ver P14C2_DEFINITION.md)."""
import io, contextlib, math, runpy
import numpy as np
from scipy.linalg import expm
with contextlib.redirect_stdout(io.StringIO()):
    N_ = runpy.run_path('P14SIM_network.py', run_name='lib')
    L = runpy.run_path('P14LOOP_sim.py', run_name='lib')
one, V3, OJS, J, PL, PR, C, Cd = N_['one'], N_['V3'], N_['OJS'], N_['J'], N_['PL'], N_['PR'], N_['C'], N_['Cd']
mk = np.array([k - 2.5 for k in range(6)])
Jp = J[0] + 1j * J[1]
X = one(Jp @ PL) @ one(Jp @ PR)
H0d = one(np.diag(np.abs(mk)).astype(complex))
bond = np.zeros((6, 6), complex); bond[2, 3] = bond[3, 2] = 0.5
Bond = one(bond)
Mhat = one(np.diag(mk).astype(complex))
proj = lambda M: V3.conj() @ M @ V3.T
LAM = 1 / 6
def HF(k, delta=0.0):
    th = 2 * math.pi * k / 6
    return H0d + OJS + LAM * (np.exp(1j * th) * X + np.exp(-1j * th) * X.conj().T) + delta * Bond
print("##### Parte 1 #####")
cm = lambda A, B: A @ B - B @ A
print(f"  |[M,O_JS]| = {np.abs(cm(Mhat, OJS)).max():.1e}; |[M,diag|m|]| = {np.abs(cm(Mhat, H0d)).max():.1e}; "
      f"|[M,X] - 2X| = {np.abs(cm(Mhat, X) - 2 * X).max():.1e}; |[M,bond]| = {np.abs(cm(Mhat, Bond)).max():.2f} (quebra explicita)")
H0 = HF(0)
for k in range(6):
    th = 2 * math.pi * k / 6
    U = expm(-1j * th / 2 * Mhat)
    dF = np.abs(U @ HF(k) @ U.conj().T - H0).max()
    dT = np.abs(proj(U @ HF(k) @ U.conj().T) - proj(H0)).max()
    print(f"  k={k}: ||U H_k U+ - H_0|| Fock = {dF:.1e}, triplas = {dT:.1e}, herm = {np.abs(HF(k)-HF(k).conj().T).max():.0e}")
print("##### Parte 2 (rede de WP-LOOP como detector) #####")
mkr = L['make_rhs']; IP = L['IP']
res = []
for k in (1, 0):
    for d in (0.0, 0.01, 0.03, 0.1, 0.3, 0.5):
        Hm = proj(HF(k, d))
        mkr.__globals__['Hphi'] = lambda phi, g=1.0, Hm=Hm: Hm
        r, ok = L['run']('L2', 0.0, 10.0, 1.0)
        res.append((k, d, r['dt'][IP]))
        print(f"  k={k} delta={d:5.2f}: Dp = {r['dt'][IP]:+.3e}  ok={ok}")
pts = [(d, abs(x)) for k, d, x in res if k == 1 and d > 0]
sl = [math.log(pts[i+1][1] / pts[i][1]) / math.log(pts[i+1][0] / pts[i][0]) for i in range(len(pts) - 1)]
print(f"  expoente local d ln|Dp|/d ln delta: {[round(s, 2) for s in sl]}")
print("##### Parte 3 #####")
kcp = (0.511e-3 / 1.2609e12) ** 2
p = sl[0]
print(f"  kappa_CP = {kcp:.2e}; com expoente p = {p:.2f}: eps ~ kappa_CP^p ~ {kcp**p:.1e}; eta ~ 0,16 eps ~ {0.16*kcp**p:.1e}")
