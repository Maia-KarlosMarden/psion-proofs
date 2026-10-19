#!/usr/bin/env python3
"""
P1 / Passo 2: tensor de vértice D1 e contrações de Hartree (H) e Fock (F)
em aritmética racional exata, para j = 1/2 .. 9/2.

Convenção (ver seção 7 do documento de especificação):
  VERTEX_CONVENTION = "D1"
  D1   : A^a = J^a P_L , B^a = J^a P_R   (geradores completos; projetor só no
         estado de entrada; J± cruza a fronteira L|R)
  D1-i : hermiticidade pelo conjugado do vértice inteiro:
         V = (1/2) sum_a [ A^a x B^a + A^a† x B^a† ]
  D1-ii: coeficiente unitário em sum_a A^a x B^a (produto escalar J.Sigma);
         o fator 1/2 da forma J+S- + J-S+ já está no produto vetorial
  D1-iii: sinal relativo escalar/JΣ: ambos reportados (+1 e -1)
  D1-iv : estrutura de Dirac: N_D = 2 e N_D = 4 reportados

Propagador: diagonal em k, peso unitário, suporte S = estados livres
(k >= 1) ou todos os estados.  Isto é uma contração de grupo, NÃO a solução
da equação de gap (isso é o passo 4).

Sinais fermiônicos: Hartree tem 2 voltas fechadas (+1), Fock 1 volta (-1).
Fatores de Dirac: Hartree ~ N_D^2, Fock ~ N_D.

Só biblioteca padrão; tudo em Fraction.
"""
from fractions import Fraction as F

VERTEX_CONVENTION = "D1"


def labels(j2):
    j = F(j2, 2)
    n = j2 + 1
    m = [F(k) - j for k in range(n)]
    L = [k for k in range(n) if m[k] < 0]
    R = [k for k in range(n) if m[k] > 0]
    return j, n, m, L, R


def contractions_JS(j2, support):
    """H, F do termo J.Sigma com vértice D1 e propagador = projetor em `support`."""
    j, n, m, L, R = labels(j2)
    Ls = [k for k in L if k in support]
    Rs = [k for k in R if k in support]
    # Hartree: so a = z tem traço (J± tem traço nulo)
    H = sum(m[k] for k in Ls) * sum(m[k] for k in Rs)
    # Fock: sum_a (J_a)_{k k'} (J_a)_{k' k} = |J+_{k',k}|^2 / 2 , apenas k' = k+1
    Fk = F(0)
    for k in Ls:
        kp = k + 1
        if kp in Rs:
            Fk += (j * (j + 1) - m[k] * (m[k] + 1)) / 2
    # D1-i: o conjugado do vértice dá as mesmas contrações (verificado por simetria)
    return H, Fk


def scalar(j2, support):
    n = len(support)
    return F(n * n), F(n)


def report(j2):
    j, n, m, L, R = labels(j2)
    allk = set(range(n))
    free = set(range(1, n))
    print(f"\n===== j = {j} =====  (convenção {VERTEX_CONVENTION})")
    for name, sup in (("todos os estados", allk), ("livres (k >= 1)", free)):
        Hs, Fs = scalar(j2, sup)
        Hj, Fj = contractions_JS(j2, sup)
        print(f"  suporte: {name}")
        print(f"    escalar  : H = {Hs}, F = {Fs}, F/H = {Fs/Hs}")
        ratio = (abs(Fj) / abs(Hj)) if Hj != 0 else None
        print(f"    J.Sigma  : H = {Hj}, F = {Fj}, |F|/|H| = "
              f"{ratio if ratio is not None else 'H=0'}")
        for ND in (2, 4):
            # coeficientes com sinais fermiônicos e Dirac
            Hs_t = ND ** 2 * Hs
            Fs_t = -ND * Fs
            Hj_t = ND ** 2 * Hj
            Fj_t = -ND * Fj
            rs = abs(Fs_t) / Hs_t
            rj = (abs(Fj_t) / abs(Hj_t)) if Hj_t != 0 else None
            print(f"      N_D={ND}: escalar H={Hs_t}, F={Fs_t}  (|F|/H={rs}) | "
                  f"JΣ H={Hj_t}, F={Fj_t}  (|F|/|H|={rj})")


def formula_check():
    """Forma fechada observada nos livres:  |F|/|H| = 2/(j-1/2)^2 (a verificar)."""
    print("\n--- verificação da forma fechada |F|/|H|(J.Sigma, livres) = 2/(j-1/2)^2 ---")
    for j2 in (3, 5, 7, 9):
        j, n, m, L, R = labels(j2)
        free = set(range(1, n))
        H, Fk = contractions_JS(j2, free)
        pred = F(2) / (j - F(1, 2)) ** 2
        print(f"  j={j}: |F|/|H| = {abs(Fk)/abs(H)} ; 2/(j-1/2)^2 = {pred} ; "
              f"{'ok' if abs(Fk)/abs(H) == pred else 'DIFERE'}")


if __name__ == "__main__":
    for j2 in (1, 3, 5, 7, 9):
        report(j2)
    formula_check()
