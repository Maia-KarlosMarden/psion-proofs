import math
import numpy as np
Q9 = np.array([1, -5, -5, 10, 4, 4]) / 9
A = 1.212e-15; LAMlen = 1.5650e-28; LAM = 1.2609e12
def r2s(t, a):
    X = np.array([[math.cos(math.pi*k/3), math.sin(math.pi*k/3)] for k in t]) * a
    q = Q9[list(t)]; c = X.mean(0); s = float(np.sum(q*np.sum((X-c)**2, 1))); Qt = q.sum()
    return s/Qt if abs(Qt) > 1e-12 else s
lim_comp = 1e5
print(f"D1 leptons carregados: escala prevista {LAM:.2e} GeV contra limite ~{lim_comp:.0e} GeV -> {'passa' if LAM > lim_comp else 'FALHA'} (ambas as variantes)")
lim_nu = 1e-36
for var, a in (("SZ-T (escala hadronica)", A), ("SZ-R (escala 1/Lambda)", LAMlen)):
    vals = {t: r2s(t, a) for t in ((1,2,3), (0,2,4), (0,1,5))}
    worst = max(abs(v) for v in vals.values())
    print(f"D2 {var}: <r2>_nu (m^2) = " + ", ".join(f"{t}:{v:+.1e}" for t, v in vals.items())
          + f" | limite {lim_nu:.0e} -> {'passa' if worst < lim_nu else f'FALHA (fator {worst/lim_nu:.0e})'}")
