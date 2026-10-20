#!/usr/bin/env python3
"""P14 / WP-TOP: B topologico ou de simetria na eq. (1)? (ver P14TOP_DEFINITION.md)."""
import io, contextlib, runpy
from fractions import Fraction as F
from itertools import combinations
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    S = runpy.run_path('P14SIM_network.py', run_name='lib')
A, B, J = S['A'], S['B'], S['J']
PL, PR = S['PL'], S['PR']
N = 6
mk = np.array([k - 2.5 for k in range(N)])


def herm_basis():
    out = []
    for i in range(N):
        for j in range(i, N):
            M = np.zeros((N, N), complex)
            if i == j:
                M[i, i] = 1
                out.append(M)
            else:
                M1 = M.copy(); M1[i, j] = M1[j, i] = 1; out.append(M1)
                M2 = M.copy(); M2[i, j] = -1j; M2[j, i] = 1j; out.append(M2)
    return out


HB = herm_basis()                       # 36 hermitianas; gerador anti-hermitiano X = iH


def flat(*Ms):
    return np.concatenate([np.r_[M.real.ravel(), M.imag.ravel()] for M in Ms])


def null(cols, tol=1e-10):
    C = np.array(cols).T
    u, s, vt = np.linalg.svd(C)
    r = int((s > tol * s.max()).sum())
    return vt[r:]


def sym_var(pairs, X):
    ad = lambda M: M @ X - X @ M
    return sum(np.kron(ad(a), b) + np.kron(a, ad(b)) + np.kron(ad(b), a) + np.kron(b, ad(a)) for a, b in pairs) / 2


if __name__ == "__main__":
    pairs = list(zip(A, B))
    print("##### Parte 1a: simetria continua dos rotulos (J.Sigma; G_UV e invariante sob qualquer U) #####")
    ns = null([flat(sym_var(pairs, 1j * H)) for H in HB])
    gens = [sum(c * H for c, H in zip(v, HB)) for v in ns]
    print(f"  dimensao = {len(gens)}")
    for G in gens:
        d = G.diagonal().real
        off = np.abs(G - np.diag(G.diagonal())).max()
        fit = np.linalg.lstsq(np.c_[np.ones(N), mk], d, rcond=None)
        print(f"   gerador diagonal? {off < 1e-10}; = {fit[0][0]:+.3f}*1 {fit[0][1]:+.3f}*m_k (residuo {np.abs(np.c_[np.ones(N), mk] @ fit[0] - d).max():.1e})")
    comm = max(np.abs(a @ b - b @ a).max() for a in gens for b in gens)
    print(f"  abeliana: {comm < 1e-10}  -> grupo U(1)_N x U(1)_M (M = soma de m_k)")

    print("\n##### Parte 1b: extensao quiral (X = iH x gamma5) #####")
    # variacao: sum {H,A_a} (x) B_a  e  sum A_a (x) {H,B_a}  (estruturas P x S independentes); G_UV exige {H,1} = 0
    cols_noUV = [flat(sum(np.kron(H @ a + a @ H, b) for a, b in pairs), sum(np.kron(a, H @ b + b @ H) for a, b in pairs)) for H in HB]
    n1 = null(cols_noUV)
    cols_UV = [np.r_[c, flat(2 * H)] for c, H in zip(cols_noUV, HB)]
    n2 = null(cols_UV)
    print(f"  geradores quirais continuos: sem G_UV = {len(n1)}, com G_UV = {len(n2)}")

    print("\n##### Parte 2: quebra pelo condensado de enrolamento J+Sigma+ + J-Sigma- #####")
    Jp = J[0] + 1j * J[1]; Jm = J[0] - 1j * J[1]
    cond = [(Jp @ PL, Jp @ PR), (Jm @ PL, Jm @ PR)]
    for lab, G in (("U(1)_N (identidade)", np.eye(N)), ("U(1)_M (diag m_k)", np.diag(mk))):
        v = np.abs(sym_var(cond, 1j * G)).max()
        print(f"  {lab:22s}: variacao do condensado = {v:.2e} -> {'preservado' if v < 1e-10 else 'QUEBRADO'}")
    print("  variedade de vacuo: U(1)_M quebrado -> S^1 (ou S^1/Z_n): pi_0 discreto, pi_1 = Z, pi_2 = pi_3 = 0")
    print("  carga topologica pontual (pi_2 ou pi_3): NENHUMA -> sem skyrmion/monopolo; so cordas e paredes")

    print("\n##### Parte 3: cargas de Noether como B #####")
    Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]
    st = {"p": (0, 4, 5), "n": (0, 1, 5), "Lam": (1, 2, 3), "Sig": (0, 2, 4), "e": (0, 1, 2), "mu": (0, 1, 2),
          "nue": (1, 2, 3), "num": (0, 2, 4), "nut": (0, 1, 5)}
    Bv = {s: 1 if s in ("p", "n", "Lam", "Sig") else 0 for s in st}
    rows = []
    for s, t in st.items():
        M = sum(F(2 * k - 5, 2) for k in t); Qs = sum(Q9[k] for k in t)
        rows.append([float(M), float(Qs), 3.0, float(Bv[s])])
        print(f"  {s:4s} {t}: M = {str(M):>5s}  Q = {str(Qs):>3s}  3Q - N/2 = {str(3*Qs - F(3, 2)):>5s}  B = {Bv[s]}")
    R = np.array(rows)
    rk = np.linalg.matrix_rank(R[:, :3]); rk2 = np.linalg.matrix_rank(R)
    print(f"  posto de (M, Q, N_psi) = {rk}; com B = {rk2} -> B no espaco das cargas de Noether: {rk2 == rk}")

    print("\n##### Parte 4: o nivel como grau de liberdade #####")
    trip = list(combinations(range(N), 3))
    neutral = [t for t in trip if sum(Q9[k] for k in t) == 0]
    print(f"  estados internos de N = 3 (rotulos distintos): {len(trip)}; neutros: {neutral}")
    print("  rotulos internos da eq. (1) por estado: (tupla) apenas; nenhum numero quantico interno de nivel.")
    print("  n e nu_tau = mesma tupla (0,1,5): diferem so por excitacao espacial -> transicao radiativa sem regra de selecao interna")
