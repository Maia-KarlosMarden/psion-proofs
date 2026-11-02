import math
import numpy as np
Q9 = np.array([1, -5, -5, 10, 4, 4]) / 9
def P(var, H, g, s=1.0):
    p = {0: (-s, H), 1: (s, H), 2: (0, g), 3: (0, -g)}
    p.update({4: (s, -H), 5: (-s, -H)} if var == "AN" else {4: (-s, -H), 5: (s, -H)})
    return p
def r2(t, p):
    X = np.array([p[k] for k in t], float); q = Q9[list(t)]; c = X.mean(0)
    v = float(np.sum(q*np.sum((X-c)**2, 1))); Qt = q.sum()
    return v/Qt if abs(Qt) > 1e-12 else v
def rat(var, H, g):
    p = P(var, H, g); rp = r2((0,4,5), p)
    return r2((0,1,5), p)/rp, r2((0,1,2), p)/rp, rp
print("Versao principal (g = 0, H/s = sqrt 3):")
for var in ("AN", "AX"):
    Rn, RS, rp = rat(var, math.sqrt(3), 0.0)
    a = math.sqrt(0.7071/rp) if rp > 0 else float('nan')
    ok = abs(Rn/-0.164-1) < 0.1 and abs(RS-0.86) < 0.21
    print(f"  {var}: R_n = {Rn:+.3f} (med -0,164) | R_S = {RS:+.3f} (med 0,86) | s = {a:.3f} fm -> {'ACERTA' if ok else 'nao acerta'}")
print("Varredura informativa:")
for var in ("AN", "AX"):
    for g in (0.0, 0.5):
        hs = np.linspace(0.2, 5, 241)
        vals = [rat(var, h, g) for h in hs]
        good = [round(h,2) for h, v in zip(hs, vals) if abs(v[0]/-0.164-1) < 0.1 and abs(v[1]-0.86) < 0.21 and v[2] > 0]
        Rns = [v[0] for v in vals if v[2] > 0]
        print(f"  {var} g/s={g}: R_n em [{min(Rns):+.2f}, {max(Rns):+.2f}]; H/s que acertam ambos: {good[:3]}{'...' if len(good)>3 else ''} ({len(good)} de {len(hs)})")
