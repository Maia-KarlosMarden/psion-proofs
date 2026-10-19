#!/usr/bin/env python3
"""
P1 / Passo 1: espaço interno explícito para o cálculo de Z_NP (PH v85).

Constrói J^a (spin j) no espaço de rótulos k = 0..2j (m = k - j), com
convenção do repositório (scripts/pmns_models.py, H7): J atua nos
constituintes L (k <= j-1/2), Sigma nos R (k >= j+1/2), mesma representação.

Verifica:  (a) [J_a,J_b] = i eps_abc J_c;  (b) Casimir J^2 = j(j+1);
           (c) Tr J_a = 0, Tr(J_a J_b) = delta_ab j(j+1)(2j+1)/3;
           (d) fecho (ou não) da álgebra dos operadores projetados P_L J P_L.
Calcula as contrações de traço relevantes para o passo 2 (Hartree vs Fock)
com propagador diagonal em k, para o termo escalar e para o termo J.Sigma.

Só biblioteca padrão. Matrizes complexas como listas de listas.
"""
import math
from fractions import Fraction as F

# ---------- álgebra linear mínima ----------
def zeros(n): return [[0j] * n for _ in range(n)]
def eye(n):
    M = zeros(n)
    for i in range(n): M[i][i] = 1 + 0j
    return M
def mm(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def add(A, B, s=1):
    n = len(A)
    return [[A[i][j] + s * B[i][j] for j in range(n)] for i in range(n)]
def scal(A, c): return [[c * x for x in r] for r in A]
def tr(A): return sum(A[i][i] for i in range(len(A)))
def comm(A, B): return add(mm(A, B), mm(B, A), -1)
def nrm(A): return max(abs(x) for r in A for x in r)
def dag(A):
    n = len(A)
    return [[A[j][i].conjugate() for j in range(n)] for i in range(n)]


def spin_ops(j2):
    """j2 = 2j (inteiro). Estados k = 0..j2, m = k - j."""
    n = j2 + 1
    j = F(j2, 2)
    Jz, Jp = zeros(n), zeros(n)
    for k in range(n):
        m = k - j
        Jz[k][k] = float(m) + 0j
        if k + 1 < n:   # J+ |k> = amp |k+1>
            Jp[k + 1][k] = math.sqrt(float(j * (j + 1) - m * (m + 1))) + 0j
    Jm = dag(Jp)
    Jx = scal(add(Jp, Jm), 0.5)
    Jy = scal(add(Jp, Jm, -1), -0.5j)
    return Jx, Jy, Jz


def proj(n, ks):
    P = zeros(n)
    for k in ks: P[k][k] = 1 + 0j
    return P


def run(j2):
    j = F(j2, 2)
    n = j2 + 1
    print(f"\n===== j = {j}  (dimensão do multipleto N = {n}) =====")
    Jx, Jy, Jz = spin_ops(j2)
    J = [Jx, Jy, Jz]
    eps = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}

    # (a) comutadores
    worst = 0.0
    for (a, b, c), s in eps.items():
        lhs = comm(J[a], J[b])
        rhs = scal(J[c], 1j * s)
        worst = max(worst, nrm(add(lhs, rhs, -1)))
    print(f"(a) max |[J_a,J_b] - i eps J_c| = {worst:.2e}")

    # (b) Casimir
    C = add(add(mm(Jx, Jx), mm(Jy, Jy)), mm(Jz, Jz))
    target = float(j * (j + 1))
    print(f"(b) Casimir = {C[0][0].real:.6f} * I ; esperado j(j+1) = {target:.6f} ; "
          f"desvio {nrm(add(C, scal(eye(n), target), -1)):.2e}")

    # (c) traços
    t1 = [abs(tr(Ja)) for Ja in J]
    exp_tr2 = float(j * (j + 1) * (2 * j + 1) / 3)
    t2 = [tr(mm(J[a], J[b])).real for a in range(3) for b in range(3)]
    print(f"(c) Tr J_a = {[round(x,12) for x in t1]} ; "
          f"Tr(J_a J_b) diag = {[round(t2[i*4],6) for i in range(3)]} "
          f"(esperado {exp_tr2:.6f}), fora-da-diagonal máx = "
          f"{max(abs(t2[1]),abs(t2[2]),abs(t2[3]),abs(t2[5]),abs(t2[6]),abs(t2[7])):.1e}")

    # (d) projetores L/R e fecho da álgebra projetada
    L = [k for k in range(n) if k - j < 0]
    R = [k for k in range(n) if k - j > 0]
    PL, PR = proj(n, L), proj(n, R)
    JL = [mm(PL, mm(Ja, PL)) for Ja in J]
    SR = [mm(PR, mm(Ja, PR)) for Ja in J]
    JLx, JLy, JLz = JL
    dev = nrm(add(comm(JLx, JLy), scal(JLz, 1j), -1))
    print(f"(d) álgebra projetada P_L J P_L: |[Jx,Jy] - i Jz| = {dev:.3f}  "
          f"({'fecha' if dev < 1e-9 else 'NÃO fecha: J_L não é gerador de SU(2)'})")

    # contrações de traço (propagador diagonal em k, peso unitário)
    # escalar: Hartree = (Tr 1)^2 ; Fock = Tr(1*1)
    def counts(ks, name):
        m = len(ks)
        print(f"    [{name}, {m} estados] escalar: Hartree = {m*m}, Fock = {m}, Fock/Hartree = {F(1,m)}")
    counts(list(range(n)), "todos")
    counts(list(range(1, n)), "livres (sem k=0)")

    # J.Sigma, J em L, Sigma em R:
    # Hartree = sum_a Tr(J_L^a) Tr(Sigma_R^a) ; Fock = sum_a Tr(J_L^a Sigma_R^a)
    def hartree(ks_L, ks_R):
        PLs, PRs = proj(n, ks_L), proj(n, ks_R)
        tot = 0
        for Ja in J:
            tot += tr(mm(PLs, Ja)) * tr(mm(PRs, Ja))
        return tot.real
    def fock(ks_L, ks_R):
        PLs, PRs = proj(n, ks_L), proj(n, ks_R)
        tot = 0
        for Ja in J:
            tot += tr(mm(mm(PLs, Ja), mm(PRs, Ja)))
        return tot.real
    allL, allR = L, R
    freeL = [k for k in L if k != 0]
    print(f"    J.Sigma (J em L, Sigma em R), propagador diagonal em k:")
    print(f"      todos os estados : Hartree = {hartree(allL, allR):+.4f} , Fock = {fock(allL, allR):+.4f}")
    print(f"      sem k=0          : Hartree = {hartree(freeL, allR):+.4f} , Fock = {fock(freeL, allR):+.4f}")
    # Fock com mistura k=2 <-> 3 (ligação do condensado de enrolamento)
    if len(L) >= 1 and len(R) >= 1:
        kb = L[-1]
        Sbond = zeros(n)
        Sbond[kb][kb + 1] = 1 + 0j
        Sbond[kb + 1][kb] = 1 + 0j
        tot = 0
        for Ja in J:
            tot += tr(mm(mm(PL, Ja), mm(Sbond, mm(PR, mm(Ja, Sbond)))))
        print(f"      Fock com ligação S_(k={kb},{kb+1}) = 1: {tot.real:+.4f} "
              f"(não nulo reativa o termo de Fock)")


if __name__ == "__main__":
    for j2 in (1, 3, 5, 7, 9):
        run(j2)
