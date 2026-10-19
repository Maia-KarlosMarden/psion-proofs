#!/usr/bin/env python3
"""P14 / WP-C: sphaleron activity at T_nuc, Ccore reproducibility, entropy dilution."""
import math
from fractions import Fraction as F
from itertools import product

LAM, XI, GSTAR, MPL, MZ = 1.2609e12, 434.09, 106.75, 1.2209e19, 91.19
ALPHA2_MZ, B2 = 0.0338, 19 / 6
EPS_PSI, N_NG = 25 / 169, 0.749
T_NUC = LAM / XI


def alpha_w(T):
    return 1 / (1 / ALPHA2_MZ + B2 / (2 * math.pi) * math.log(T / MZ))


def hubble(T, gstar=GSTAR):
    return 1.66 * math.sqrt(gstar) * T ** 2 / MPL


def sph_ratio(T, kappa):
    return kappa * alpha_w(T) ** 5 * T ** 4 / T ** 3 / hubble(T)


print("##### Passo 2: atividade do sphaleron, Gamma/(T^3 H) #####")
print(f"  T_nuc = Lambda/xi* = {T_NUC:.3e} GeV")
for T, name in ((T_NUC, "T_nuc"), (LAM, "Lambda_meta"), (1e12, "1e12"), (1e2, "1e2")):
    r = [sph_ratio(T, k) for k in (15, 18, 21)]
    print(f"  {name:12s} T={T:.2e}  alpha_W={alpha_w(T):.4f}  "
          f"Gamma/(T^3H) = {r[1]:.3g} (kappa 15..21: {r[0]:.3g}..{r[2]:.3g})  "
          f"{'ATIVO' if r[1] > 1 else 'FORA DE EQUILIBRIO'}")
# temperatura onde a razao cruza 1
lo, hi = 1e6, 1e15
for _ in range(100):
    mid = math.sqrt(lo * hi)
    if sph_ratio(mid, 18) > 1:
        lo = mid
    else:
        hi = mid
print(f"  T_max(sphaleron ativo, kappa=18) = {lo:.2e} GeV")

print("\n##### Passo 3: reprodutibilidade de Ccore = 8.4434e-4 #####")
target = 8.4434e-4
cands = {}
gS_opts = {"g*(4/11)": GSTAR * 4 / 11, "3.91": 3.91, "g*": GSTAR, "43/11": 43 / 11, "2": 2.0,
           "g*(4/11)^(4/3)": GSTAR * (4 / 11) ** (4 / 3)}
pref_opts = {"8/23": 8 / 23, "28/79": 28 / 79, "1": 1.0}
for (gn, gS), (pn, pf) in product(gS_opts.items(), pref_opts.items()):
    for en, e in (("4/3", 4 / 3), ("1", 1.0), ("0", 0.0)):
        v = pf * GSTAR / gS ** e * 1e-4
        cands[f"{pn} * g*/({gn})^{en} * 1e-4"] = v
        v2 = pf * (GSTAR / gS) ** e * 1e-4
        cands[f"{pn} * (g*/{gn})^{en} * 1e-4"] = v2
best = sorted(cands.items(), key=lambda kv: abs(kv[1] / target - 1))[:5]
for k, v in best:
    print(f"  {k:45s} = {v:.4e}  ({100*(v/target-1):+.1f}%)")
hit = [k for k, v in cands.items() if abs(v / target - 1) < 0.01]
print(f"  reproduzido a 1%: {hit if hit else 'NENHUMA combinacao'}")
print(f"  eq.(14): Ccore * eps_psi = {target*EPS_PSI:.3e}; para dar 6.02e-10 o fator "
      f"(g_1/2/g_S) teria de ser {6.02e-10/(target*EPS_PSI):.2e}, nao O(1)")

print("\n##### Passo 4: diluicao de entropia #####")
gS0 = 3.91
dil = GSTAR / gS0
print(f"  g_S(T_nuc)/g_S(hoje) = {GSTAR}/{gS0} = {dil:.1f}")
print(f"  n_f/n_b por grau de liberdade = 3/4 = 0.75 ; paper n/n_gamma = {N_NG}")
print("  A eq.(33) nao contem fator de diluicao explicito.")

print("\n##### Passo 5: eta_B pela cadeia WP-B com xi* = 434 (provisorio) #####")
cs = F(28, 79)
BL = {"S1 (+1/9)": F(1, 9), "S3 (+3/7)": F(3, 7)}
print("  eta_B = c_s (B-L) * eps_psi * (n/n_gamma) * xi*^-3  [x 1/3 por psion] [/ diluicao]")
for name, bl in BL.items():
    base = float(cs * bl) * EPS_PSI * N_NG / XI ** 3
    print(f"  {name}: por composto {base:.2e} | por psion {base/3:.2e} | "
          f"por composto/diluicao {base/dil:.2e} | por psion/diluicao {base/3/dil:.2e}")
print(f"  eq.(33) original (4/9): {4/9*EPS_PSI*N_NG/XI**3:.2e} ; observado 6.13e-10")
print("\n  xi* que a cadeia WP-B exigiria para 6.13e-10 (so para registro, nao e previsao):")
for name, bl in BL.items():
    for lab, extra in (("por composto", 1.0), ("por psion", 1 / 3), ("por psion/diluicao", 1 / 3 / dil)):
        xi_req = (float(cs * bl) * EPS_PSI * N_NG * extra / 6.13e-10) ** (1 / 3)
        print(f"    {name} {lab:20s}: xi* = {xi_req:.0f}")
