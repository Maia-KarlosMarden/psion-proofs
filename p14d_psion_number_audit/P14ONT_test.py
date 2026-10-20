#!/usr/bin/env python3
"""P14 / WP-ONT: contabilidade de N_psi e carga em reacoes observadas, ontologias O+ e O-."""
from fractions import Fraction as F
import itertools

Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]
inv = lambda t: tuple(sorted((6 - k) % 6 for k in t))
qt = lambda t: sum(Q9[k] for k in t)
TUP = {"e-": (0, 1, 2), "mu-": (0, 1, 2), "p": (0, 4, 5), "n": (0, 1, 5),
       "nu_e": (1, 2, 3), "nu_mu": (0, 2, 4), "nu_tau": (0, 1, 5), "Sigma+": (1, 3, 5)}
ANTI = {"e+": "e-", "mu+": "mu-", "pbar": "p", "nbar": "n", "nubar_e": "nu_e", "nubar_mu": "nu_mu"}
BOSONS = {"gamma": 0, "pi0": 0, "pi+": 1, "pi-": -1, "X17": 0}     # cargas observadas


def props(name, ont):
    """(N_psi, Q, assinatura de identidade) de um fermion composto."""
    if name in TUP:
        t = TUP[name]
        return 3, qt(t), ("psi", t)
    base = ANTI[name]
    t = TUP[base]
    if ont == "O+":
        ti = inv(t)
        return 3, qt(ti), ("psi", ti)
    return -3, -qt(t), ("psibar", t)


REACTIONS = {
    "R1 n -> p e- nubar_e": (["n"], ["p", "e-", "nubar_e"]),
    "R2 p -> n e+ nu_e": (["p"], ["n", "e+", "nu_e"]),
    "R3 e- p -> n nu_e": (["e-", "p"], ["n", "nu_e"]),
    "R4 e+ e- -> 2 gamma": (["e+", "e-"], ["gamma", "gamma"]),
    "R5 e+ e- -> 3 gamma": (["e+", "e-"], ["gamma", "gamma", "gamma"]),
    "R6 gamma -> e+ e-": (["gamma"], ["e+", "e-"]),
    "R7 p pbar -> pi+ pi- pi0": (["p", "pbar"], ["pi+", "pi-", "pi0"]),
    "R8 pi0 -> 2 gamma": (["pi0"], ["gamma", "gamma"]),
    "R9 pi+ -> mu+ nu_mu": (["pi+"], ["mu+", "nu_mu"]),
    "R10 mu- -> e- nubar_e nu_mu": (["mu-"], ["e-", "nubar_e", "nu_mu"]),
    "R11 X17 -> e+ e-": (["X17"], ["e+", "e-"]),
}
UNK = list(BOSONS)


def rref(rows, n):
    A = [list(r) for r in rows]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    return A[:r], piv


def solve(eqs, fixed):
    """eqs: lista de (coef dict boson->int, const). sistema sum coef*x + const = 0."""
    unk = [u for u in UNK if u not in fixed]
    n = len(unk)
    rows = []
    for coef, const in eqs:
        c = const + sum(coef.get(u, 0) * fixed[u] for u in fixed)
        rows.append([F(coef.get(u, 0)) for u in unk] + [F(c)])
    A, piv = rref(rows, n + 1)
    inconsistent = n in piv
    sol = {}
    if not inconsistent:
        for i, p in enumerate(piv):
            sol[unk[p]] = -A[i][n]
    return (not inconsistent), sol


for ont in ("O+", "O-"):
    print(f"\n========== Ontologia {ont} ==========")
    eqs, qok = [], {}
    for rname, (lhs, rhs) in REACTIONS.items():
        coef, const, q = {}, 0, F(0)
        for side, sgn in ((lhs, 1), (rhs, -1)):
            for x in side:
                if x in BOSONS:
                    coef[x] = coef.get(x, 0) + sgn
                    q += sgn * BOSONS[x]
                else:
                    N, Qx, _ = props(x, ont)
                    const += sgn * N
                    q += sgn * Qx
        eqs.append((coef, const))
        qok[rname] = (q == 0)
    for mode, fixed in (("bosons livres", {}), ("fiel ao paper (X17 com N_psi = 0)", {"X17": 0})):
        ok, sol = solve(eqs, fixed)
        print(f"  T-N [{mode}]: {'CONSISTENTE' if ok else 'INCONSISTENTE'}"
              + (f"  N_psi dos bosons = { {k: str(v) for k, v in sol.items()} }" if ok else ""))
        if not ok:
            # menor subconjunto conflitante (pares e trios)
            names = list(REACTIONS)
            for k in (2, 3):
                hits = []
                for comb in itertools.combinations(range(len(names)), k):
                    if not solve([eqs[i] for i in comb], fixed)[0]:
                        hits.append(" + ".join(names[i].split()[0] for i in comb))
                if hits:
                    print(f"      conflitos minimos de {k} reacoes: {hits[:8]}")
                    break
    bad = [r for r, v in qok.items() if not v]
    print(f"  T-Q carga conservada em todas: {not bad}" + (f"  falhas: {bad}" if bad else ""))
    print("  T-A antiparticulas:")
    for a, base in (("e+", "e-"), ("pbar", "p"), ("nbar", "n")):
        Na, Qa, ida = props(a, ont); Nb, Qb, idb = props(base, ont)
        print(f"     {a:5s}: N_psi={Na:+d} Q={str(Qa):>3s}  (part.: Q={str(Qb):>3s})  "
              f"carga oposta: {Qa == -Qb}  distinta da particula: {ida != idb}")
    # Sigma+ (paper: (1,3,5) Sigma+-like); antiparticula observada tem Q = -1
    t = TUP["Sigma+"]
    qa = qt(inv(t)) if ont == "O+" else -qt(t)
    print(f"     Sigma+bar: Q = {qa} (observado -1): {qa == -1}")
