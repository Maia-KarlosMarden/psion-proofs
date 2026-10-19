#!/usr/bin/env python3
"""P14 / WP-B: N_psi, B, L accounting for the Z6 three-body census (exact)."""
from fractions import Fraction as F
from itertools import combinations

Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]
mJ = [F(2 * k - 5, 2) for k in range(6)]

triples = [t for t in combinations(range(6), 3) if sum(t) % 6 in (0, 3)]
assert len(triples) == 8, triples
print("##### Censo: triplas observáveis (I in {0,3}) #####")
for t in triples:
    print(f"  {t}  I={sum(t)%6}  Q={sum(Q9[k] for k in t)}  sum mJ={sum(mJ[k] for k in t)}")

schemes = {
    "S1 (paper: 9 familias, geracoes separadas, exotico fora)": {
        "lep": {(0, 1, 2): 3, (1, 2, 3): 1},
        "bar": {(0, 4, 5): 1, (0, 1, 5): 1, (0, 2, 4): 1, (1, 3, 5): 1, (2, 3, 4): 1},
    },
    "S2 (8 triplas, duais = neutrinos)": {
        "lep": {(0, 1, 2): 1, (0, 1, 5): 1, (0, 2, 4): 1, (1, 2, 3): 1},
        "bar": {(0, 4, 5): 1, (1, 3, 5): 1, (2, 3, 4): 1},
    },
    "S3 (8 triplas, duais = barions)": {
        "lep": {(0, 1, 2): 1, (1, 2, 3): 1},
        "bar": {(0, 4, 5): 1, (0, 1, 5): 1, (0, 2, 4): 1, (1, 3, 5): 1, (2, 3, 4): 1},
    },
}


def BL(t, sch):
    return (1 if t in sch["bar"] else 0), (1 if t in sch["lep"] else 0)


print("\n##### Passo 2: recombinacoes 2->2 que conservam o multiconjunto de k #####")
pairs = list(combinations(triples, 2)) + [(t, t) for t in triples]
rearr = []
for a, b in pairs:
    ms = sorted(a + b)
    for c, d in pairs:
        if (c, d) <= (a, b):
            continue
        if sorted(c + d) == ms:
            rearr.append((a, b, c, d))
print(f"  processos encontrados: {len(rearr)}")
for name, sch in schemes.items():
    nB = nL = nBL = 0
    for a, b, c, d in rearr:
        Bi = BL(a, sch)[0] + BL(b, sch)[0]; Li = BL(a, sch)[1] + BL(b, sch)[1]
        Bf = BL(c, sch)[0] + BL(d, sch)[0]; Lf = BL(c, sch)[1] + BL(d, sch)[1]
        nB += Bi != Bf; nL += Li != Lf; nBL += (Bi - Li) != (Bf - Lf)
    print(f"  {name}: violam B: {nB}, violam L: {nL}, violam B-L: {nBL} de {len(rearr)}")
for a, b, c, d in rearr:
    print(f"    {a}+{b} <-> {c}+{d}")

print("\n##### Passo 3: B, L, B-L em unidades de Delta n_psi/3 #####")
res = {}
for name, sch in schemes.items():
    fam = {**sch["lep"], **sch["bar"]}
    N = sum(fam.values())
    B = F(sum(sch["bar"].values()), N); L = F(sum(sch["lep"].values()), N)
    # peso 2J+1, J=|sum mJ|, por tripla (geracoes com o mesmo peso da tripla)
    w = {t: (2 * abs(sum(mJ[k] for k in t)) + 1) * m for t, m in fam.items()}
    W = sum(w.values())
    Bw = sum(w[t] for t in sch["bar"]) / W; Lw = sum(w[t] for t in sch["lep"]) / W
    res[name] = (B, L)
    print(f"  {name}:\n     democratico: B={B} L={L} B-L={B-L}\n     peso 2J+1 : B={Bw} L={Lw} B-L={Bw-Lw}")

print("\n##### Passo 4: conversao por sphaleron #####")
def cs(nf, nh): return F(8 * nf + 4 * nh, 22 * nf + 13 * nh)
print(f"  c_s(SM, N_f=3, N_H=1) = {cs(3,1)} = {float(cs(3,1)):.4f}")
hits = [(nf, nh) for nf in range(1, 5) for nh in range(0, 5) if cs(nf, nh) == F(8, 23)]
print(f"  8/23 de §13.3 corresponde a (N_f, N_H) = {hits}  -> {float(F(8,23)):.4f}")
print("\n  fator que substitui eps_conf = 4/9 na eq. (33):  c_s * (B-L)  [por psion: /3 adicional]")
for name, (B, L) in res.items():
    f = cs(3, 1) * (B - L)
    print(f"  {name}: c_s(B-L) = {f} = {float(f):+.4f}  ;  /(4/9) = {float(f/F(4,9)):+.3f}"
          f"  ;  por psion /3: {float(f/3):+.4f} (/(4/9) = {float(f/3/F(4,9)):+.3f})")
