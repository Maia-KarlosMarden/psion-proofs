#!/usr/bin/env python3
"""P14-D: pre-checagem de escalas (ver P14D_DEFINITION.md)."""
import math

LAM = 1.2609e12
G2 = 2 * math.pi ** 2 / 3
MPL = 1.2209e19
MZ6 = 2.996e-3
EPS_PSI = 25 / 169
T_SPH_MAX, T_SPH_MIN = 1.2e11, 130.0
YEAR_GEV = 3.156e7 * 1.519e24          # 1 ano em GeV^-1

# linhas resolvidas de WP-A passo 2: (lambda, c, T_c, T_eq, T_nuc, xi)
ROWS = [
    ("1/25", 1.0, 1.389e11, 1.38e11, 3.57e9, 38.8),
    ("1/25", 0.3, 2.535e11, 1.38e11, 2.32e10, 10.9),
    ("eps_psi", 0.3, 4.876e11, 1.91e11, 1.24e11, 3.9),
    ("eps_psi", 3.0, 1.542e11, 1.91e11, 4.32e10, 3.6),
    ("1/Nc", 1.0, 2.835e11, 1.97e11, 7.18e10, 3.9),
    ("1/Nc", 3.0, 1.637e11, 1.97e11, 3.54e10, 4.6),
    ("b_e", 1.0, 3.688e11, 2.24e11, 4.35e10, 8.5),
    ("b_e", 3.0, 2.129e11, 2.24e11, 2.70e10, 7.9),
]


def gstar(T):
    return 106.75 if T > 100 else (10.75 if T > 5e-4 else 3.36)


def hubble(T):
    return 1.66 * math.sqrt(gstar(T)) * T * T / MPL


def gamma_rearr(T, c):
    return c * G2 ** 2 * T ** 5 / LAM ** 4


def T_of_t(t):
    # radiacao: t = 0.301 g*^{-1/2} M_Pl / T^2  (itera g*)
    T = 1e-3
    for _ in range(20):
        T = math.sqrt(0.301 * MPL / (math.sqrt(gstar(T)) * t))
    return T


def sph(T):
    return "ativo" if T_SPH_MIN < T < T_SPH_MAX else "inativo"


print("##### 1. Versao 'paredes Z6' (K.10) #####")
T_ann = T_of_t(4 * YEAR_GEV)
for name, T in (("formacao (T = m_Z6)", MZ6), ("aniquilacao (t = 4 anos)", T_ann)):
    r = [gamma_rearr(T, c) / hubble(T) for c in (1e-3, 1.0)]
    print(f"  {name:28s} T = {T:.3e} GeV  Gamma_rearr/H = {r[0]:.1e} .. {r[1]:.1e}  sphaleron {sph(T)}")
print("  -> nenhuma violacao de B ativa enquanto as paredes Z6 existem.")

print("\n##### 2. Versao 'bolhas da transicao P14-c' #####")
print(f"{'lambda':>8s} {'c':>4s} {'T_nuc':>9s} {'T_RH':>9s} {'alpha':>9s} {'D_TI':>9s} "
      f"{'Gam/H(Tn)':>18s} {'eta_cota (c_G=1e-3..1)':>24s}")
best = 0.0
for lam, c, Tc, Teq, Tn, xi in ROWS:
    alpha = (Teq / Tn) ** 4
    D = (Tn / Teq) ** 3
    gh = [gamma_rearr(Tn, cg) / hubble(Tn) for cg in (1e-3, 1.0)]
    eta = [45 / (2 * math.pi ** 2 * 106.75) * EPS_PSI * math.pi / 3 * gamma_rearr(Tn, cg) / Tn * D
           for cg in (1e-3, 1.0)]
    best = max(best, eta[1])
    print(f"{lam:>8s} {c:4.1f} {Tn:9.2e} {Teq:9.2e} {alpha:9.2e} {D:9.2e} "
          f"{gh[0]:8.1e}..{gh[1]:8.1e} {eta[0]:11.1e}..{eta[1]:9.1e}")
print(f"\n  maior eta_cota (c_Gamma = 1): {best:.1e} ; criterio para seguir: >= 1e-10")
print("  escala: eta_cota ∝ (T_nuc/Lambda)^4 (T_nuc/T_RH)^3  ~  T_nuc^7 a T_RH fixo")
print(f"\n  comparacao com a estimativa de cabeca registrada (1e-17 para lambda=1/25, c=1)")
