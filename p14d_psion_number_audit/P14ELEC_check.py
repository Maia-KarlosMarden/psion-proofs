"""Checagem exploratoria (nao pre-registrada): como J.Sigma age em cada tripla nao confinada."""
import io, contextlib, runpy
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    N = runpy.run_path('P14SIM_network.py', run_name='lib')
V3, OJS, TR = N['V3'], N['OJS'], N['TR']
Ot = V3.conj() @ OJS @ V3.T
for t in [t for t in TR if sum(t) % 6 in (0, 3)]:
    i = TR.index(t)
    print(f"{t} n_L={sum(k<=2 for k in t)} M={sum(k-2.5 for k in t):+.1f} <J.Sigma>={Ot[i,i].real:+.3f} "
          f"acoplamento={np.linalg.norm(np.delete(Ot[i], i)):.3f}")
