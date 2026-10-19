#!/usr/bin/env python3
"""P1 / Passo 4c: D1' (ver P1_D1PRIME_DEFINITION.md).  j = 5/2 apenas."""
import math
import itertools
import numpy as np
import P_2pi_step4 as s4
import P_2pi_step4b as s4b

Q9 = [1, -5, -5, 10, 4, 4]


def graded_pairs(J, PL, PR, Q):
    n = len(Q)
    out = []
    for Ja in J:
        A, B = Ja @ PL, Ja @ PR
        for q in range(-30, 31):
            Aq = np.zeros((n, n), complex)
            Bq = np.zeros((n, n), complex)
            for k in range(n):
                for l in range(n):
                    if abs(A[k, l]) > 0 and Q[k] - Q[l] == q:
                        Aq[k, l] = A[k, l]
                    if abs(B[k, l]) > 0 and Q[k] - Q[l] == -q:
                        Bq[k, l] = B[k, l]
            if np.abs(Aq).max() > 0 and np.abs(Bq).max() > 0:
                out.append((Aq, Bq))
    return out


def make_Wq(j2, ND, sigma):
    n, j, J = s4.spin_mats(j2)
    L, R = s4.LR(n, j)
    PL, PR = s4.proj(n, L), s4.proj(n, R)
    pairs = graded_pairs(J, PL, PR, Q9)

    def W(s):
        HS = np.trace(s) ** 2
        FS = np.trace(s @ s)
        HJ = sum(np.trace(a @ s) * np.trace(b @ s) for a, b in pairs)
        FJ = sum(np.trace(a @ s @ b @ s) for a, b in pairs)
        return (HS - FS / ND).real + sigma * (HJ.real - FJ.real / ND)
    return W, n, j


def retained_weight():
    n, j, J = s4.spin_mats(5)
    L, R = s4.LR(n, j)
    PL, PR = s4.proj(n, L), s4.proj(n, R)
    full = sum(np.abs(Ja @ PL).sum() * 0 for Ja in J)
    tot = kept = 0.0
    for Ja in J:
        A, B = Ja @ PL, Ja @ PR
        for k in range(n):
            for l in range(n):
                for k2 in range(n):
                    for l2 in range(n):
                        w = abs(A[k, l] * B[k2, l2]) ** 2
                        if w == 0:
                            continue
                        tot += w
                        if (Q9[k] - Q9[l]) + (Q9[k2] - Q9[l2]) == 0:
                            kept += w
    return kept / tot


ALLOWED = {(k, l) for k, l in itertools.combinations(range(1, 6), 2) if Q9[k] == Q9[l]}


def build_prime(j2, ND, sigma):
    W, n, j = make_Wq(j2, ND, sigma)
    support = list(range(1, n))
    basis, nr = s4.hermitian_basis(n, support, ALLOWED)
    B = np.array([E / math.sqrt(x) for E, x in zip(basis, nr)])
    d = len(B)
    A = np.zeros((d, d))
    W1 = [W(b) for b in B]
    for i in range(d):
        A[i, i] = W1[i]
        for k in range(i + 1, d):
            A[i, k] = A[k, i] = 0.5 * (W(B[i] + B[k]) - W1[i] - W1[k])
    return n, j, W, support, B, A


def main():
    print("##### PASSO 4c — D1' (vértice e condensado neutros em carga), j = 5/2 #####")
    print(f"pares permitidos no condensado: {sorted(ALLOWED)}")
    print(f"peso do vértice D1 que sobrevive ao filtro de carga: {retained_weight():.3f} "
          f"(o violador era 0.598)")
    # verificação: o vértice filtrado não tem componentes violadores (por construção)
    n, j, J = s4.spin_mats(5)
    L, R = s4.LR(n, j)
    pairs = graded_pairs(J, s4.proj(n, L), s4.proj(n, R), Q9)
    print(f"componentes (A^q, B^-q) retidos: {len(pairs)}")
    s4b.build = build_prime
    b0 = (4 * 2.5 + 3) / (2 * 2.5 + 1)
    Nc = 6.0
    res = {}
    for sigma in (+1, -1):
        for ND in (2, 4):
            out = s4b.analyse(5, ND, sigma, [1e-2, 1e-3])
            lam = out["lam"]
            Wp = 5 * lam
            Z = 2 ** 5 / Wp
            be = b0 / (Z * Nc)
            res[(sigma, ND)] = (lam, Z, be)
            print(f"  s={sigma:+d} N_D={ND}: lam_max={lam:.5f} (deg {out['deg']})  "
                  f"W'={Wp:.4f}  Z'={Z:.4f}  b_e'={be:.5f}")
            for name, delta, r in out["rows"]:
                if r is None:
                    print(f"      [{name}] delta={delta:g}: sem solução")
                else:
                    R_, ov, mu, Wv, it, ok = r
                    print(f"      [{name}] delta={delta:g}: R/lam-1={R_/lam-1:+.2e}  "
                          f"overlap={ov:.6f}  conv={ok}")
    print("\n##### COMPARAÇÃO FINAL (só agora) #####")
    print("  alvo: Z_NP = 32/25 = 1.2800 ; b_e = 325/1152 = 0.28212")
    for (sigma, ND), (lam, Z, be) in res.items():
        flag = "  <-- dentro de 1%" if abs(Z / 1.28 - 1) <= 0.01 else ""
        print(f"  s={sigma:+d} N_D={ND}: Z'={Z:.4f} ({100*(Z/1.28-1):+.1f}%)  "
              f"b_e'={be:.5f} ({100*(be/(325/1152)-1):+.1f}%){flag}")


if __name__ == "__main__":
    main()
