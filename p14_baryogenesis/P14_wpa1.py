#!/usr/bin/env python3
"""P14 / WP-A passo 1: ordem da transicao e teto de xi* = T_c/T_0 em M1, M2, M3."""
import math
import numpy as np
from scipy import integrate, optimize

LAM = 1.2609e12
SIG0 = 2.9961e-3
NC = 6

# ---------------- M1: NJL campo medio -----------------
def gap_integral(s, L=LAM):
    """I(s) - I(0) = int_0^L p^2 dp (1/E - 1/p) exato a O(s^4/L^2): s^2/4 - (s^2/2) ln(2L/s)."""
    return s * s / 4 - s * s / 2 * math.log(2 * L / s)


def V_T(s, T, invG, nd):
    """(2 pi^2)[V(s,T)-V(0,T)] com o termo Lambda^2 s^2 cancelado analiticamente:
       s^2 g(sigma0)/2 - [s^4/32 - (s^4/8) ln(2L/s)] - th(s,T),  g = gap_integral.
       invG e ignorado (fixado pelo gap em sigma0)."""
    if s == 0:
        return 0.0
    s = abs(s)
    vac_res = s ** 4 / 32 - s ** 4 / 8 * math.log(2 * LAM / s)
    f = lambda p: p * p * 2 * T * (math.log1p(math.exp(-math.sqrt(p * p + s * s) / T)) - math.log1p(math.exp(-p / T)))
    th = integrate.quad(f, 0, 40 * T + 10 * s, limit=300)[0]
    return s * s * gap_integral(SIG0) / 2 - vac_res - th


def m1():
    print("##### M1: NJL campo medio, corte Lambda_meta, sigma_0 = 2.9961 MeV #####")
    # gap: 1/(4G') = I(sigma0)/(2 pi^2), com G' absorvendo degenerescencia
    # dV/ds = s*invG/2 - s*I(s)/(2pi^2) = 0 em s = sigma0  ->  invG = 2 I(sigma0)/(2 pi^2)
    # (primeira execucao usava 4*I: erro de normalizacao, corrigido; ver REGISTRY)
    invG = 2 * (LAM ** 2 / 2 + gap_integral(SIG0)) / (2 * math.pi ** 2)
    g0 = optimize.minimize_scalar(lambda s: V_T(s, 1e-9 * SIG0, invG, 4), bounds=(0, 3 * SIG0), method="bounded").x
    print(f"  verificacao do gap a T=0: sigma_min/sigma_0 = {g0/SIG0:.4f} (deve ser 1)")
    # coeficientes de Landau em s -> 0 a temperatura T: a(T) = d2V/ds2|0 ; b(T) = d4V/ds4|0 /24
    def a(T):
        h = SIG0 * 0.05
        return (V_T(h, T, invG, 4) - 2 * V_T(0, T, invG, 4) + V_T(-h, T, invG, 4)) / h ** 2
    def b(T):
        h = SIG0 * 0.1
        f = [V_T(k * h, T, invG, 4) for k in (-2, -1, 0, 1, 2)]
        d4 = (f[0] - 4 * f[1] + 6 * f[2] - 4 * f[3] + f[4]) / h ** 4
        return d4 / 24
    Tc = optimize.brentq(a, 0.1 * SIG0, 20 * SIG0, xtol=1e-9)
    Tc_formula = SIG0 * math.sqrt(3 / math.pi ** 2 * (math.log(2 * LAM / SIG0) - 0.5))
    print(f"  T_c = {Tc*1e3:.4f} MeV ; T_c/sigma_0 = {Tc/SIG0:.4f} ; formula (3/pi^2)(ln(2L/s0)-1/2): {Tc_formula/SIG0:.4f}")
    print(f"  coeficiente quartico b(T_c) = {b(Tc):+.3e}  -> {'POSITIVO: 2a ordem' if b(Tc) > 0 else 'negativo: 1a ordem'}")
    # espinodal: a(T_0) = 0 -> para 2a ordem T_0 = T_c
    print("  T_0 (espinodal da fase simetrica) = T_c  ->  teto xi* = T_c/T_0 = 1.00")
    # checar que o minimo abaixo de T_c cresce continuamente (sem salto)
    for frac in (0.999, 0.99, 0.9):
        T = frac * Tc
        res = optimize.minimize_scalar(lambda s: V_T(s, T, invG, 4), bounds=(0, 1.5 * SIG0), method="bounded")
        print(f"    T = {frac} T_c : sigma_min/sigma_0 = {res.x/SIG0:.4f}")
    print(f"  Escala: T_c = {Tc*1e3:.2f} MeV vs T_nuc(paper) = {LAM/434:.2e} GeV : razao {LAM/434/Tc:.1e}")
    print("  -> o condensado de 3 MeV nao e a transicao da bariogenese.")


# ---------------- M2: Potts q=3 campo medio -----------------
def m2():
    print("\n##### M2: Potts q=3 em campo medio (fator Z3 de V.5) #####")
    q = 3
    Tc = (q - 2) / (2 * (q - 1) * math.log(q - 1))   # em unidades de zJ
    T0 = 1 / q
    print(f"  T_c/zJ = {Tc:.4f} ; T_0/zJ = {T0:.4f} ; teto xi* = T_c/T_0 = {Tc/T0:.4f}")
    # verificacao numerica: f(m) = -(1/2) z J [ (1+(q-1)m)^2 + (q-1)(1-m)^2 ]/q^2 ... usar forma padrao
    def free(m, T):
        p1 = (1 + (q - 1) * m) / q
        p2 = (1 - m) / q
        e = -0.5 * (p1 ** 2 + (q - 1) * p2 ** 2)
        s = -(p1 * math.log(p1) + (q - 1) * p2 * math.log(p2)) if 0 < m < 1 else -math.log(1 / q) * (-1) if m == 0 else 0
        return e - T * s
    for T in (Tc * 1.02, Tc, Tc * 0.98, T0 * 1.01, T0 * 0.99):
        ms = np.linspace(1e-6, 0.999, 4000)
        f = np.array([free(m, T) for m in ms])
        f0 = free(1e-9, T)
        i = f.argmin()
        bar = f[: i].max() - f0 if i > 5 else 0.0
        print(f"    T/zJ={T:.4f}: m_min={ms[i]:.3f}  f_min-f_0={f[i]-f0:+.4e}  barreira={bar:+.2e}")


# ---------------- M3: CW conforme (exploratorio, P14-c) -----------------
def m3():
    print("\n##### M3 (EXPLORATORIO, nao esta no paper): CW conforme + massa termica #####")
    print("  V = (lam/4) s^4 [ln(s^2/L^2) - 1/2] + (c/2) T^2 s^2 ; T_0 = 0, sem teto.")
    print("  xi* pela acao de bounce S3/T = 4 ln(M_Pl/T) + ... (criterio Gamma = H^4), grade em lam:")
    MPL = 1.2209e19
    c = 1.0
    for lam in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01):
        # T_c: V(min) = 0 em T_c; minimo em s*: aproximacao numerica
        def Vmin(T):
            def V(s): return lam / 4 * s ** 4 * (math.log(s * s / LAM ** 2) - 0.5) + 0.5 * c * T * T * s * s
            r = optimize.minimize_scalar(V, bounds=(1e-3 * LAM, 3 * LAM), method="bounded")
            return r.fun
        Tc = optimize.brentq(Vmin, 1e-6 * LAM, 5 * LAM)
        # bounce em aproximacao de parede fina nao vale; usar escala: S3/T ~ A / (lam * ln(Tc/T)) com A ~ 50 (ordem de grandeza)
        # criterio: S3/T = 4 ln(MPl/T) - ... ~ 140 em T ~ 1e9
        def nuc(lnT):
            T = math.exp(lnT)
            S3T = 50.0 / (lam * max(math.log(Tc / T), 1e-9))
            return S3T - (4 * math.log(MPL / T) - 1.5 * math.log(S3T) if S3T > 1 else 1)
        try:
            lnTn = optimize.brentq(nuc, math.log(Tc) - 40, math.log(Tc) - 1e-6)
            xi = Tc / math.exp(lnTn)
            print(f"    lam={lam:5.2f}: T_c={Tc:.2e} GeV  xi*~{xi:8.1f}  (escala de ordem de grandeza)")
        except ValueError:
            print(f"    lam={lam:5.2f}: T_c={Tc:.2e} GeV  sem nucleacao (inflacao termica)")
    print("  A faixa-alvo de WP-C (44 a 302) corresponde a lam de ordem 0.02 a 0.1 neste modelo.")
    print("  Nada no paper fixa lam; isto e um postulado novo (P14-c).")


if __name__ == "__main__":
    m1(); m2(); m3()
