#!/usr/bin/env python3
"""P1 / Passo 5: Fock-space trace and the numerator 2^(2j) (ver P1_STEP5_DEFINITION.md)."""
import math
from fractions import Fraction as Fr
import numpy as np


def ops(n):
    sm = np.array([[0, 1], [0, 0]], float)      # annihilation on one mode
    Z = np.diag([1.0, -1.0])
    I = np.eye(2)
    cs = []
    for k in range(n):
        m = np.array([[1.0]])
        for l in range(n):
            m = np.kron(m, Z if l < k else (sm if l == k else I))
        cs.append(m)
    return cs


def t1():
    print("##### T1: espaço de Fock, estado maximamente misto #####")
    print(" n   Tr1   <N>   <N(N-1)>   n(n-1)/4   Hartree n^2/4  Fock -n/4  4<N>^2")
    for n in range(1, 10):
        cs = ops(n)
        d = 2 ** n
        N = sum(c.T @ c for c in cs)
        V = N @ (N - np.eye(d))
        # V como soma de c†_k c†_l c_l c_k, conferência direta
        V2 = sum(cs[k].T @ cs[l].T @ cs[l] @ cs[k] for k in range(n) for l in range(n))
        assert np.allclose(V, V2)
        trN, trV = np.trace(N) / d, np.trace(V) / d
        print(f" {n}  {int(round(np.trace(np.eye(d))))}  {trN:5.2f}  {trV:8.3f}  {n*(n-1)/4:8.3f}"
              f"  {n*n/4:10.3f}  {-n/4:10.3f}  {4*trN**2:7.2f}")


def t2():
    print("\n##### T2: normalizações do quociente  (numerador) / (4 <N>^2 ou (Tr N)^2) #####")
    print(" j    n  (a) ambos norm.   (b) ambos nao norm.   (c) hibrido   2^n/n^2 paper")
    for j2 in (1, 3, 5, 7, 9):
        n = j2
        a = 1 / n ** 2                                  # Tr rho =1 ; 4<N>^2 = n^2
        trN = n * 2 ** (n - 1)
        b = 2 ** n / trN ** 2                           # rho = 1
        c = 2 ** n / n ** 2
        print(f" {j2}/2  {n}   {a:12.5f}   {b:16.6e}   {c:10.4f}   {2**n/n**2:10.4f}")


def t3():
    n = 5
    print("\n##### T3: Gibbs, energia uniforme eps nos n=5 modos: Z/(2<N>)^2 #####")
    print("  beta*eps   Z=(1+e^-x)^n   <N>=n/(e^x+1)   Z/(2<N>)^2")
    for x in (0.0, 0.01, 0.1, 0.5, 1.0, 3.0, 10.0):
        Z = (1 + math.exp(-x)) ** n
        N = n / (math.exp(x) + 1)
        print(f"  {x:7.2f}   {Z:12.5f}   {N:12.5f}   {Z/(2*N)**2:12.4f}")
    print("  (valor do paper: 32/25 = 1.2800 só em beta*eps -> 0)")


def t4():
    print("\n##### T4: be(j) = b0/(Z_NP Nc) vs fórmula fechada #####")
    for j2 in (1, 3, 5, 7, 9):
        j = Fr(j2, 2)
        b0 = (4 * j + 3) / (2 * j + 1)
        Nc = 2 * j + 1
        Z = Fr(2 ** j2, j2 ** 2)
        be = b0 / (Z * Nc)
        closed = (4 * j + 3) * (2 * j) ** 2 / (2 ** j2 * (2 * j + 1) ** 2)
        print(f"  j={j}: be={be} = {float(be):.5f}  (fechada igual: {be == closed})")


if __name__ == "__main__":
    t1()
    t2()
    t3()
    t4()
