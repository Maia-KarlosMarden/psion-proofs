#!/usr/bin/env python3
"""P14 / WP-D0: auditoria exata de cargas conservadas (ver P14D0_DEFINITION.md)."""
from fractions import Fraction as F
from itertools import combinations

K = range(6)
L, R = (0, 1, 2), (3, 4, 5)
Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]
M = [F(2 * k - 5, 2) for k in K]


def jp_nonzero(k):   # <k+1|J+|k> != 0
    return 0 <= k < 5


def jm_nonzero(k):   # <k-1|J-|k> != 0
    return 0 < k <= 5


def shifts():
    """(k, l, d): psion k in L -> k+d, psion l in R -> l-d, d = +-1, elementos nao nulos."""
    out = []
    for k in L:
        for l in R:
            for d in (+1, -1):
                ok_k = jp_nonzero(k) if d == 1 else jm_nonzero(k)
                ok_l = jm_nonzero(l) if d == 1 else jp_nonzero(l)
                if ok_k and ok_l and 0 <= k + d <= 5 and 0 <= l - d <= 5:
                    out.append((k, l, d))
    return out


def dQ(k, l, d):
    return Q9[k + d] - Q9[k] + Q9[l - d] - Q9[l]


# ---------- algebra linear exata ----------
def rref(rows, n):
    A = [list(r) for r in rows]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    return A[:r], piv


def nullspace(rows, n):
    """base de {c : c . row = 0 para toda row}."""
    if not rows:
        return [[F(int(i == j)) for j in range(n)] for i in range(n)]
    A, piv = rref([[F(x) for x in r] for r in rows], n)
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = [F(0)] * n
        v[f] = F(1)
        for i, p in enumerate(piv):
            v[p] = -A[i][f]
        basis.append(v)
    return basis


def in_span(v, basis, n):
    if not basis:
        return all(x == 0 for x in v)
    r1 = len(rref(basis, n)[0])
    r2 = len(rref(basis + [v], n)[0])
    return r1 == r2


def fmt(v):
    return "(" + ", ".join(str(x) for x in v) + ")"


# ---------- nivel de constituintes ----------
def e(k):
    v = [0] * 6
    v[k] = 1
    return v


def delta_shift(k, l, d):
    v = [0] * 6
    v[k] -= 1; v[l] -= 1; v[k + d] += 1; v[l - d] += 1
    return v


SH = shifts()
proc_D1 = [delta_shift(*s) for s in SH]
proc_D1p = [delta_shift(*s) for s in SH if dQ(*s) == 0]
proc_cond = [[0, 0, -1, 1, 0, 0]]

charges6 = {
    "N_psi": [F(1)] * 6,
    "M = sum m_J": M,
    "Q (eletrica)": Q9,
    "n_R (quiralidade R)": [F(int(k >= 3)) for k in K],
}

print("##### 1. Nivel de constituintes #####")
print(f"  processos JSigma D1: {len(proc_D1)} ; D1' (Delta Q = 0): {len(proc_D1p)}")
sets6 = {
    "S_sym (JSigma D1)": proc_D1,
    "S_sym' (JSigma D1')": proc_D1p,
    "fase quebrada livre (D1 + condensado)": proc_D1 + proc_cond,
    "fase quebrada livre (D1' + condensado)": proc_D1p + proc_cond,
}
K6 = {}
for name, P in sets6.items():
    B = nullspace(P, 6)
    K6[name] = B
    print(f"  {name}: dim K = {len(B)}")
    for v in B:
        print(f"      base: {fmt(v)}")
    for cn, cv in charges6.items():
        print(f"      {cn:22s} conservada: {in_span(cv, B, 6)}")

# ---------- nivel de compostos ----------
TR = list(combinations(K, 3))
IDX = {t: i for i, t in enumerate(TR)}
NT = len(TR)
OBS = [t for t in TR if sum(t) % 6 in (0, 3)]


def vec_comp(pairs):
    v = [0] * NT
    for t, s in pairs:
        v[IDX[t]] += s
    return v


def canon(labels):
    s = tuple(sorted(labels))
    return s if len(set(s)) == 3 else None


proc_rearr, proc_js_intra, proc_js_inter, proc_cond_c = [], [], [], []
for a, b in combinations(TR, 2):
    ms = sorted(a + b)
    for c, d in combinations(TR, 2):
        if (c, d) != (a, b) and {c, d} != {a, b} and sorted(c + d) == ms:
            proc_rearr.append(vec_comp([(a, -1), (b, -1), (c, 1), (d, 1)]))
for t in TR:
    for (k, l, d) in SH:
        if k in t and l in t:
            new = canon([x for x in t if x not in (k, l)] + [k + d, l - d])
            if new:
                proc_js_intra.append(((k, l, d), vec_comp([(t, -1), (new, 1)])))
    for k_from, k_to in ((2, 3), (3, 2)):
        if k_from in t:
            new = canon([x for x in t if x != k_from] + [k_to])
            if new:
                proc_cond_c.append(vec_comp([(t, -1), (new, 1)]))
for t1 in TR:
    for t2 in TR:
        if t1 == t2:
            continue
        for (k, l, d) in SH:
            if k in t1 and l in t2:
                n1 = canon([x for x in t1 if x != k] + [k + d])
                n2 = canon([x for x in t2 if x != l] + [l - d])
                if n1 and n2:
                    proc_js_inter.append(((k, l, d), vec_comp([(t1, -1), (t2, -1), (n1, 1), (n2, 1)])))


def js_list(filt):
    return [v for (s, v) in proc_js_intra + proc_js_inter if filt(s)]


JS_D1 = js_list(lambda s: True)
JS_D1p = js_list(lambda s: dQ(*s) == 0)


def additive(w):
    return [sum(w[k] for k in t) for t in TR]


SCHEMES = {
    "S1": ({(0, 1, 2), (1, 2, 3)}, {(0, 4, 5), (0, 1, 5), (0, 2, 4), (1, 3, 5), (2, 3, 4)}),
    "S2": ({(0, 1, 2), (0, 1, 5), (0, 2, 4), (1, 2, 3)}, {(0, 4, 5), (1, 3, 5), (2, 3, 4)}),
    "S3": ({(0, 1, 2), (1, 2, 3)}, {(0, 4, 5), (0, 1, 5), (0, 2, 4), (1, 3, 5), (2, 3, 4)}),
}
chargesC = {cn: additive(cv) for cn, cv in charges6.items()}
for sn, (lep, bar) in SCHEMES.items():
    chargesC[f"B ({sn})"] = [F(int(t in bar)) for t in TR]
    chargesC[f"L ({sn})"] = [F(int(t in lep)) for t in TR]
    chargesC[f"B-L ({sn})"] = [F(int(t in bar)) - F(int(t in lep)) for t in TR]

print(f"\n##### 2. Nivel de compostos (20 triplas; {len(OBS)} observaveis) #####")
print(f"  processos: rearranjo {len(proc_rearr)}, JSigma intra {len(proc_js_intra)}, "
      f"inter {len(proc_js_inter)}, condensado {len(proc_cond_c)}")
setsC = {
    "S_brk' (rearranjo + JSigma D1, sem condensado)": proc_rearr + JS_D1,
    "S_brk  (rearranjo + JSigma D1 + condensado)": proc_rearr + JS_D1 + proc_cond_c,
    "S_brk' (rearranjo + JSigma D1', sem condensado)": proc_rearr + JS_D1p,
    "S_brk  (rearranjo + JSigma D1' + condensado)": proc_rearr + JS_D1p + proc_cond_c,
}
KC = {}
for name, P in setsC.items():
    B = nullspace(P, NT)
    KC[name] = B
    print(f"  {name}: dim K = {len(B)}")
    for cn, cv in chargesC.items():
        print(f"      {cn:22s} conservada: {in_span(cv, B, NT)}")

print("\n##### 3. Cargas candidatas: conservadas na fase quebrada, violadas na simetrica #####")
pairs = [
    ("S_sym (JSigma D1)", "S_brk  (rearranjo + JSigma D1 + condensado)"),
    ("S_sym (JSigma D1)", "S_brk' (rearranjo + JSigma D1, sem condensado)"),
    ("S_sym' (JSigma D1')", "S_brk  (rearranjo + JSigma D1' + condensado)"),
    ("S_sym' (JSigma D1')", "S_brk' (rearranjo + JSigma D1', sem condensado)"),
]
verdict = {}
for s_name, b_name in pairs:
    pull = [additive(v) for v in K6[s_name]]          # imagem das cargas simetricas
    cand = [v for v in KC[b_name] if not in_span(v, pull, NT)]
    # dimensao do quociente K(S_brk) / imagem
    dim_q = len(rref(KC[b_name] + pull, NT)[0]) - len(rref(pull, NT)[0]) if pull else len(KC[b_name])
    verdict[(s_name, b_name)] = dim_q
    print(f"  {s_name}  ->  {b_name}")
    print(f"      dim K_brk = {len(KC[b_name])}, dim imagem de K_sym = {len(rref(pull, NT)[0]) if pull else 0}, "
          f"candidatas (dim do quociente) = {dim_q}")
    for v in cand[:4]:
        nz = {TR[i]: x for i, x in enumerate(v) if x != 0}
        print(f"        exemplo: {nz}")

main = verdict[("S_sym (JSigma D1)", "S_brk  (rearranjo + JSigma D1 + condensado)")]
print(f"\n  VEREDITO (criterio pre-registrado, caso principal D1 + condensado): "
      f"{'VIAVEL' if main > 0 else 'INVIAVEL'}")
