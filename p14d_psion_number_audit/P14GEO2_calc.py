import math, itertools
import numpy as np
Q9 = np.array([1, -5, -5, 10, 4, 4]) / 9
A = 1.212
def P(t): return np.array([[math.cos(math.pi*k/3), math.sin(math.pi*k/3)] for k in t])
def r2(t):
    X = P(t); q = Q9[list(t)]; c = X.mean(0); s = float(np.sum(q*np.sum((X-c)**2, 1))); Qt = q.sum()
    return (s/Qt if abs(Qt) > 1e-12 else s) * A**2
print("##### Parte 1: raios de carga previstos (fm^2) #####")
for nome, t in (("p (controle)", (0,4,5)), ("n (controle)", (0,1,5)), ("Sigma- (0,1,2), minha identificacao", (0,1,2)),
                ("Sigma+ (1,3,5)", (1,3,5)), ("Sigma+ (2,3,4)", (2,3,4)), ("Lambda (1,2,3)", (1,2,3)),
                ("Sigma0 (0,2,4)", (0,2,4)), ("(3,4,5) Q=+2", (3,4,5))):
    v = r2(t); print(f"  {nome:38s} <r2> = {v:+.3f} fm^2" + (f"  (r = {math.sqrt(v):.3f} fm)" if v > 0 else ""))
def apex(t):
    s = sorted(t)
    for i in range(3):
        a, b, c = s[i], s[(i+1)%3], s[(i+2)%3]
        if (b-a) % 6 == 1 and (c-b) % 6 == 1: return b, (a, c)
    return None, None
def mus(t, rule):
    ap, outer = apex(t)
    choices = [(ap, outer)] if ap is not None else [(c, tuple(x for x in t if x != c)) for c in t]
    out = []
    for c, (a, b) in choices:
        out.append((2/3)*(Q9[a]+Q9[b]) - (1/3)*Q9[c] if rule == "VA" else Q9[c])
    return out
meas = {"n": -0.685, "Lambda": -0.2195, "Sigma+": 0.880, "Sigma-": -0.415}
tup = {"p": [(0,4,5)], "n": [(0,1,5)], "Lambda": [(1,2,3)], "Sigma+": [(1,3,5), (2,3,4)], "Sigma-": [(0,1,2)]}
print("##### Parte 2: momentos magneticos #####")
for rule in ("VA", "VB"):
    mp = mus((0,4,5), rule)[0]
    ok = True; line = f"  {rule}: mu_p = {mp:+.3f}"
    for nm, m in meas.items():
        vals = sorted({round(v/mp, 3) for t in tup[nm] for v in mus(t, rule)})
        hit = any(abs(v/m - 1) < 0.2 for v in vals)
        ok &= hit
        line += f" | {nm}/p: {vals} (med {m:+.3f}) {'ok' if hit else 'X'}"
    print(line + f"  -> {'ACERTA' if ok else 'nao acerta'}")
