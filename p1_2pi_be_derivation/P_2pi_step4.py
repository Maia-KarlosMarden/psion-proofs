#!/usr/bin/env python3
"""
P1 / Passo 4: peso crítico W_crit e Z_ren, conforme P1_STEP3_DEFINITION.md
(SHA-256 4e0ad99b...2a2c).  Convenção VERTEX_CONVENTION = "D1".

Procedimento (pré-registrado):
  W[s] = (H_S - F_S/N_D) + sigma*(H_J - F_J/N_D),  s = matriz de pesos no espaço
  de rótulos (Hermitiana, suportada nos estados livres).
    H_S = (Tr s)^2        F_S = Tr(s s)
    H_J = Re sum_a Tr(A_a s) Tr(B_a s)     F_J = Re sum_a Tr(A_a s B_a s)
    A_a = J_a P_L , B_a = J_a P_R           (D1; D1-i: Re = vértice + conjugado)
  W é forma quadrática em s; a equação de gap linearizada em G_c dá a direção
  crítica s_c = autovetor de maior autovalor de W (no produto de Frobenius).
  Normalização ADICIONADA ANTES DE RODAR (registro D3-N): Tr s_c = n_livres.
  W_crit = W[s_c normalizado].   Z_ren = 2^(n_livres) / W_crit.
  b_e = b0 / (Z_ren * N_c).

Convenção de sinal (registro D4-G, antes de rodar): G > 0 é atrativo no canal
escalar de Hartree; portanto a direção crítica é a de maior autovalor POSITIVO.

Nenhum passo usa 325/1152, 32/25, m_e ou rho_Lambda; a comparação com 32/25 é
feita só na impressão final.

Desvios do pré-registro, declarados:
  - Teste de esquema: em regularização dimensional sem massa não há ponto crítico
    (divergência quadrática ausente), então o teste (1) fica "não aplicável";
    o regulador só reescala G_c e não entra em W_crit.
"""
import itertools
import math
import numpy as np
from fractions import Fraction as Fr

VERTEX_CONVENTION = "D1"


def spin_mats(j2):
    n = j2 + 1
    j = j2 / 2
    Jz = np.zeros((n, n), complex)
    Jp = np.zeros((n, n), complex)
    for k in range(n):
        m = k - j
        Jz[k, k] = m
        if k + 1 < n:
            Jp[k + 1, k] = math.sqrt(j * (j + 1) - m * (m + 1))
    Jm = Jp.conj().T
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / 2j
    return n, j, [Jx, Jy, Jz]


def LR(n, j):
    L = [k for k in range(n) if k - j < 0]
    R = [k for k in range(n) if k - j > 0]
    return L, R


def proj(n, ks):
    P = np.zeros((n, n), complex)
    for k in ks:
        P[k, k] = 1
    return P


def make_W(j2, ND, sigma):
    n, j, J = spin_mats(j2)
    L, R = LR(n, j)
    PL, PR = proj(n, L), proj(n, R)
    A = [Ja @ PL for Ja in J]
    B = [Ja @ PR for Ja in J]

    def W(s):
        HS = np.trace(s) ** 2
        FS = np.trace(s @ s)
        HJ = sum(np.trace(a @ s) * np.trace(b @ s) for a, b in zip(A, B))
        FJ = sum(np.trace(a @ s @ b @ s) for a, b in zip(A, B))
        # D1-i: vértice + conjugado, meia soma -> parte real
        return (HS - FS / ND).real + sigma * (HJ.real - FJ.real / ND)
    return W, n, j


def hermitian_basis(n, support, allowed_pairs=None):
    """base real do espaço de matrizes Hermitianas suportadas em `support`."""
    basis, norms = [], []
    for k in support:
        E = np.zeros((n, n), complex)
        E[k, k] = 1
        basis.append(E)
        norms.append(1.0)
    for k, l in itertools.combinations(support, 2):
        if allowed_pairs is not None and (k, l) not in allowed_pairs:
            continue
        E1 = np.zeros((n, n), complex)
        E1[k, l] = E1[l, k] = 1
        E2 = np.zeros((n, n), complex)
        E2[k, l] = 1j
        E2[l, k] = -1j
        basis += [E1, E2]
        norms += [2.0, 2.0]
    return basis, np.array(norms)


def critical(W, n, support, allowed_pairs=None):
    basis, nr = hermitian_basis(n, support, allowed_pairs)
    d = len(basis)
    Q = np.zeros((d, d))
    W1 = [W(E) for E in basis]
    for i in range(d):
        Q[i, i] = W1[i]
        for k in range(i + 1, d):
            Q[i, k] = Q[k, i] = 0.5 * (W(basis[i] + basis[k]) - W1[i] - W1[k])
    D = np.diag(1 / np.sqrt(nr))
    w, V = np.linalg.eigh(D @ Q @ D)
    order = np.argsort(w)[::-1]
    lam = w[order[0]]
    degenerate = (len(w) > 1 and abs(w[order[0]] - w[order[1]]) < 1e-9)
    y = V[:, order[0]]
    x = D @ y
    s = sum(xi * E for xi, E in zip(x, basis))
    tr = np.trace(s).real
    info = dict(lam=lam, gap=(w[order[0]] - w[order[1]]) if len(w) > 1 else None,
                degenerate=degenerate, trace_unit=tr, dim=d)
    if abs(tr) < 1e-9:
        info["status"] = "traço nulo: normalização impossível"
        return None, info
    s = s * (len(support) / tr)          # normalização D3-N: Tr s = n_livres
    info["status"] = "ok"
    return W(s), info


def fr(x, md=200):
    return str(Fr(x).limit_denominator(md))


def run_family():
    print("\n##### TESTE PRINCIPAL: suporte = estados livres (k >= 1) #####")
    for j2 in (1, 3, 5, 7, 9):
        n, j, J = spin_mats(j2)
        nl = n - 1
        print(f"\n=== j = {fr(j)}  (n_livres = {nl}, 2^n = {2**nl}, (2j)^2 = {(2*j)**2:.0f}) ===")
        for sigma in (+1, -1):
            for ND in (2, 4):
                W, _, _ = make_W(j2, ND, sigma)
                Wc, info = critical(W, n, list(range(1, n)))
                if Wc is None:
                    print(f"  s={sigma:+d} N_D={ND}: {info['status']} (lam={info['lam']:.4f})")
                    continue
                Z = 2 ** nl / Wc
                Nc = 2 * j + 1
                b0 = (4 * j + 3) / (2 * j + 1)
                be = b0 / (Z * Nc)
                tag = "" if info["lam"] > 0 else "  [lam<=0: sem G_c>0]"
                tag += "  [DEGENERADO]" if info["degenerate"] else ""
                print(f"  s={sigma:+d} N_D={ND}: lam_max={info['lam']:.4f}  W_crit={Wc:.4f} "
                      f"(~{fr(Wc)})  Z_ren={Z:.4f}  b_e={be:.5f}{tag}")


def run_anchor_included():
    print("\n##### TESTE 3: âncora k = 0 incluída (j = 5/2) #####")
    j2 = 5
    n, j, J = spin_mats(j2)
    for sigma in (+1, -1):
        for ND in (2, 4):
            W, _, _ = make_W(j2, ND, sigma)
            Wc, info = critical(W, n, list(range(n)))
            if Wc is None:
                print(f"  s={sigma:+d} N_D={ND}: {info['status']}")
                continue
            print(f"  s={sigma:+d} N_D={ND}: lam_max={info['lam']:.4f}  W_crit={Wc:.4f}  "
                  f"Z_ren={2**n/Wc:.4f}  (2^6/W)")


def run_charge_tests():
    print("\n##### TESTE 4: carga elétrica e Z6 (j = 5/2) #####")
    j2 = 5
    n, j, J = spin_mats(j2)
    Q9 = [1, -5, -5, 10, 4, 4]
    L, R = LR(n, j)
    # componentes do vértice sum_a J_a x J_a = JzJz + 1/2 (J+J- + J-J+)
    comps = []
    for k in L:
        for l in R:
            comps.append(("JzJz", (k, k), (l, l), 1.0))
    def jp(k):
        m = k - j
        return math.sqrt(j * (j + 1) - m * (m + 1)) if k + 1 < n else 0.0
    def jm(k):
        m = k - j
        return math.sqrt(j * (j + 1) - m * (m - 1)) if k - 1 >= 0 else 0.0
    for k in L:
        for l in R:
            if k + 1 < n and l - 1 >= 0:
                comps.append(("J+J-", (k, k + 1), (l, l - 1), 0.5 * jp(k) * jm(l)))
            if k - 1 >= 0 and l + 1 < n:
                comps.append(("J-J+", (k, k - 1), (l, l + 1), 0.5 * jm(k) * jp(l)))
    tot = 0.0
    viol = 0.0
    nviol = 0
    ksum_ok = True
    for name, (k, kp), (l, lp), amp in comps:
        w = amp ** 2
        tot += w
        dQ = (Q9[kp] - Q9[k]) + (Q9[lp] - Q9[l])
        if dQ != 0:
            viol += w
            nviol += 1
        if (kp + lp) != (k + l):
            ksum_ok = False
    print(f"  Z6 (soma de k conservada em todos os componentes): {ksum_ok}")
    print(f"  carga elétrica Q_k: componentes com Delta Q != 0: {nviol} de {len(comps)}; "
          f"peso violador / peso total = {viol/tot:.3f}")
    return viol / tot


def run_charge_restricted():
    print("\n##### VARIANTE (tentativa adicional): condensado só entre estados de mesma carga #####")
    j2 = 5
    n, j, J = spin_mats(j2)
    Q9 = [1, -5, -5, 10, 4, 4]
    allowed = {(k, l) for k, l in itertools.combinations(range(1, n), 2) if Q9[k] == Q9[l]}
    print(f"  pares fora da diagonal permitidos (k,l), k,l >= 1: {sorted(allowed)}")
    for sigma in (+1, -1):
        for ND in (2, 4):
            W, _, _ = make_W(j2, ND, sigma)
            Wc, info = critical(W, n, list(range(1, n)), allowed_pairs=allowed)
            if Wc is None:
                print(f"  s={sigma:+d} N_D={ND}: {info['status']}")
                continue
            print(f"  s={sigma:+d} N_D={ND}: lam_max={info['lam']:.4f}  W_crit={Wc:.4f}  "
                  f"Z_ren={2**(n-1)/Wc:.4f}")


def sanity_vs_step2():
    """W[s = projetor livre] deve reproduzir o passo 2 (peso unitário)."""
    print("##### VERIFICAÇÃO: peso unitário reproduz o passo 2 #####")
    j2 = 5
    n, j, J = spin_mats(j2)
    s = proj(n, range(1, n))
    for sigma in (+1, -1):
        for ND in (2, 4):
            W, _, _ = make_W(j2, ND, sigma)
            print(f"  s={sigma:+d} N_D={ND}: W = {W(s):.4f} = {fr(W(s))}")


if __name__ == "__main__":
    sanity_vs_step2()
    run_family()
    run_anchor_included()
    run_charge_tests()
    run_charge_restricted()
    print("\n##### COMPARAÇÃO FINAL (só agora) em j = 5/2 #####")
    print("  alvo: (2j)^2 = 25 ; Z_NP = 32/25 = 1.2800 ; b_e = 325/1152 = 0.28212")
