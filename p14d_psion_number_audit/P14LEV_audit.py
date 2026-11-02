#!/usr/bin/env python3
"""P14 / WP-LEV1: cargas que a regra de niveis precisa conservar (ver P14LEV_DEFINITION.md)."""
import math
import sympy as sp

ST = ["p", "n", "Lam", "Sig", "e", "mu", "nue", "num", "nut"]
CLASS = {"p": "H", "n": "H", "Lam": "H", "Sig": "H", "e": "E", "mu": "E", "nue": "nu", "num": "nu", "nut": "nu"}
TUP = {"p": (0, 4, 5), "n": (0, 1, 5), "Lam": (1, 2, 3), "Sig": (0, 2, 4), "e": (0, 1, 2), "mu": (0, 1, 2),
       "nue": (1, 2, 3), "num": (0, 2, 4), "nut": (0, 1, 5)}
Q = {"p": 1, "n": 0, "Lam": 0, "Sig": 0, "e": -1, "mu": -1, "nue": 0, "num": 0, "nut": 0}
BQ = {"g": 0, "pi0": 0, "pi+": 1, "pi-": -1}


def vec(lhs, rhs):
    """lhs/rhs: listas de nomes; prefixo '~' = antiparticula; bosons em BQ."""
    v = {s: 0 for s in ST}; q = 0; npsi = 0
    for side, sg in ((lhs, 1), (rhs, -1)):
        for x in side:
            if x in BQ:
                q += sg * BQ[x]; continue
            a = -1 if x.startswith("~") else 1
            s = x.lstrip("~")
            v[s] += sg * a; q += sg * a * Q[s]; npsi += sg * a * 3
    return [v[s] for s in ST], q, npsi


OBS = {
    "n -> p e- nubar_e": (["n"], ["p", "e", "~nue"]),
    "p -> n e+ nu_e (beta+)": (["p"], ["n", "~e", "nue"]),
    "e- p -> n nu_e": (["e", "p"], ["n", "nue"]),
    "pi+ -> mu+ nu_mu": (["pi+"], ["~mu", "num"]),
    "mu- -> e- nubar_e nu_mu": (["mu"], ["e", "~nue", "num"]),
    "nu_mu n -> mu- p": (["num", "n"], ["mu", "p"]),
    "nubar_e p -> e+ n": (["~nue", "p"], ["~e", "n"]),
    "Lambda -> p pi-": (["Lam"], ["p", "pi-"]),
    "Lambda -> n pi0": (["Lam"], ["n", "pi0"]),
    "Sigma0 -> Lambda gamma": (["Sig"], ["Lam", "g"]),
    "oscilacao nu_mu -> nu_e": (["num"], ["nue"]),
    "oscilacao nu_e -> nu_tau": (["nue"], ["nut"]),
}
FORB = {
    "p -> nu pi+": (["p"], ["num", "pi+"]),
    "n -> nu pi0": (["n"], ["num", "pi0"]),
    "n -> nu nu nubar (invisivel)": (["n"], ["num", "nut", "~num"]),
    "H(p e-) -> nu_mu nu_tau": (["p", "e"], ["num", "nut"]),
    "n -> nu_tau gamma (mesma tupla)": (["n"], ["nut", "g"]),
    "Lambda -> nu_e gamma (mesma tupla)": (["Lam"], ["nue", "g"]),
    "Sigma0 -> nu_mu gamma (mesma tupla)": (["Sig"], ["num", "g"]),
}

if __name__ == "__main__":
    print("##### Parte 1: cargas aditivas conservadas exigidas pelas reacoes observadas #####")
    rows = []
    for name, (l, r) in OBS.items():
        v, q, npsi = vec(l, r)
        rows.append(v)
        print(f"  {name:28s} Q conservada: {q == 0}  N_psi conservado: {npsi == 0}")
    M = sp.Matrix(rows)
    ns = M.nullspace()
    print(f"\n  dimensao do espaco de cargas conservadas = {len(ns)}")
    for b in ns:
        b = b / max(abs(x) for x in b)
        print("   ", {s: str(b[i]) for i, s in enumerate(ST)})
    span = sp.Matrix.hstack(*ns)
    def inspan(c):
        c = sp.Matrix(c)
        return sp.Matrix.hstack(span, c).rank() == span.rank()
    B = [1 if CLASS[s] == "H" else 0 for s in ST]
    L = [1 if CLASS[s] != "H" else 0 for s in ST]
    N3 = [1] * len(ST)
    Le = [1 if s in ("e", "nue") else 0 for s in ST]
    print(f"  N_psi/3 conservado: {inspan(N3)}; B (classe H): {inspan(B)}; L (classes E, nu): {inspan(L)}; "
          f"L_e separado: {inspan(Le)}")

    print("\n##### Parte 2: processos permitidos por N_psi e Q, nao observados #####")
    for name, (l, r) in FORB.items():
        v, q, npsi = vec(l, r)
        dB = -sum(v[i] * B[i] for i in range(len(ST)))
        dL = -sum(v[i] * L[i] for i in range(len(ST)))
        same = all(TUP[x.lstrip('~')] == TUP[r[0].lstrip('~')] for x in l) and len(l) == 1 and r[1] == "g"
        print(f"  {name:36s} N_psi ok: {npsi == 0}  Q ok: {q == 0}  dB = {dB:+d}  dL = {dL:+d}  "
              f"mesma tupla: {same}")

    print("\n##### Parte 3: supressao exigida #####")
    HBAR, YR = 6.582e-25, 3.156e7
    mN, alpha, Lam = 0.939, 1 / 137.036, 1.2387e12
    G_em = alpha * mN
    G_6 = mN ** 5 / Lam ** 4
    for name, Gnat, lim in (("p -> nu pi+   (contato dim-6)", G_6, 3.9e32),
                            ("n -> nu pi0   (contato dim-6)", G_6, 1.1e33),
                            ("n -> 3 nu     (contato dim-6)", G_6, 5.8e29),
                            ("H -> nu nu    (contato dim-6, ref. 1e29 anos)", G_6, 1e29),
                            ("n -> nu_tau g (EM, mesma tupla, ref. 1e29 anos)", G_em, 1e29)):
        tau = HBAR / Gnat / YR
        S = math.sqrt(tau / lim)
        print(f"  {name:48s} tau natural = {tau:.1e} anos; limite {lim:.1e} anos -> supressao de amplitude S <= {S:.1e}")
