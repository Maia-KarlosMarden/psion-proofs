import math
import numpy as np
from scipy.optimize import brentq
Q9 = np.array([1, -5, -5, 10, 4, 4]) / 9
j = 2.5
def pos(k, G, h=1.0):
    if G == "G1": return np.array([math.cos(math.pi*k/3), math.sin(math.pi*k/3), 0.0])
    if G == "G5":
        c = k % 3; return np.array([math.cos(2*math.pi*c/3), math.sin(2*math.pi*c/3), -h if k <= 2 else h])
    if G == "G6":
        m = k - 2.5; rho = math.sqrt(j*(j+1) - m*m)
        return np.array([rho*math.cos(math.pi*k/3), rho*math.sin(math.pi*k/3), m])
def r2(t, G, h=1.0):
    X = np.array([pos(k, G, h) for k in t]); q = Q9[list(t)]; c = X.mean(0)
    s = float(np.sum(q*np.sum((X-c)**2, 1))); Qt = q.sum()
    return s/Qt if abs(Qt) > 1e-12 else s
def ratios(G, h=1.0):
    rp = r2((0,4,5), G, h); return r2((0,1,5), G, h)/rp, r2((0,1,2), G, h)/rp, rp
def report(G, h=None):
    Rn, RS, rp = ratios(G, h if h is not None else 1.0)
    okn = abs(Rn/-0.164 - 1) < 0.10; okS = abs(RS - 0.86) < 0.21
    a = math.sqrt(0.7071/rp) if rp > 0 else float('nan')
    tag = f" h/a = {h:.4f}" if h is not None else ""
    print(f"  {G}{tag}: R_n = {Rn:+.3f} {'ok' if okn else 'X'} | R_S = {RS:+.3f} {'ok' if okS else 'X'} | a = {a:.3f} fm -> {'ACERTA' if okn and okS else 'nao acerta'}")
report("G1"); report("G6")
hs = np.linspace(0, 5, 5001); f = [ratios("G5", h)[0] + 0.164 for h in hs]
roots = [brentq(lambda x: ratios("G5", x)[0] + 0.164, hs[i], hs[i+1]) for i in range(len(hs)-1) if f[i]*f[i+1] < 0]
print(f"  G5: R_n(h/a) varia de {min(f)-0.164:+.3f} a {max(f)-0.164:+.3f}; raizes para R_n = -0,164: {[round(r,4) for r in roots]}")
for r in roots: report("G5", r)
for h in (0.0, 0.5, 1.0): report("G5", h)
