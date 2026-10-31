#!/usr/bin/env python3
"""P14-D adendo POS-HOC (nao pre-registrado): lavagem apos reaquecimento.
T_f: temperatura em que Gamma_rearr = H (radiacao, g* = 106.75)."""
import math
import runpy, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    _m = runpy.run_path('P14D_scales.py')
ROWS, gamma_rearr, hubble, LAM, G2, MPL = (_m[k] for k in ('ROWS','gamma_rearr','hubble','LAM','G2','MPL'))
for c in (1e-5, 1e-4, 1e-3, 1e-2, 1.0):
    Tf = (1.66 * math.sqrt(106.75) * LAM ** 4 / (c * G2 ** 2 * MPL)) ** (1 / 3)
    print(f"c_Gamma = {c:7.0e}:  T_f = {Tf:.2e} GeV")
print()
print(f"{'lambda':>8s} {'c':>4s} {'T_RH':>9s}  Gamma/H(T_RH) para c_Gamma = 1e-3 .. 1   c_Gamma maximo sem lavagem")
for lam, c, Tc, Teq, Tn, xi in ROWS:
    r = [gamma_rearr(Teq, cg) / hubble(Teq) for cg in (1e-3, 1.0)]
    cmax = hubble(Teq) / gamma_rearr(Teq, 1.0)
    print(f"{lam:>8s} {c:4.1f} {Teq:9.2e}  {r[0]:9.1e} .. {r[1]:9.1e}          {cmax:.1e}")
