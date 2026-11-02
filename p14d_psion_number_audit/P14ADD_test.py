#!/usr/bin/env python3
"""P14 / WP-ADD1: o menor ingrediente que conserva B (ver P14ADD_DEFINITION.md)."""
import io, contextlib, math, runpy
from fractions import Fraction as F
from itertools import combinations

with contextlib.redirect_stdout(io.StringIO()):
    L = runpy.run_path('P14LEV_audit.py', run_name='lib')
ST, vec, OBS, FORB = L['ST'], L['vec'], L['OBS'], L['FORB']
SPEC = {s: ("B" if s in ("p", "n", "Lam", "Sig") else "L") for s in ST}
Q9 = [F(1, 9), F(-5, 9), F(-5, 9), F(10, 9), F(4, 9), F(4, 9)]

if __name__ == "__main__":
    print("##### Parte 1: catalogo com duas especies #####")
    ok_obs = True
    for tag, cat in (("observadas", OBS), ("nao observadas", FORB)):
        print(f"  {tag}:")
        for name, (l, r) in cat.items():
            v, q, npsi = vec(l, r)
            dNB = sum(3 * v[i] for i, s in enumerate(ST) if SPEC[s] == "B")
            dNL = sum(3 * v[i] for i, s in enumerate(ST) if SPEC[s] == "L")
            cons = dNB == 0 and dNL == 0 and q == 0
            if tag == "observadas":
                ok_obs &= cons
            print(f"    {name:38s} dN_B = {dNB:+d}  dN_L = {dNL:+d}  dQ = {q:+d}  -> {'permitido' if cons else 'PROIBIDO'}")
    print(f"  todas as observadas permitidas: {ok_obs}")

    print("\n##### Parte 2: espectro de tres corpos por especie #####")
    trip = list(combinations(range(6), 3))
    qs = {t: sum(Q9[k] for k in t) for t in trip}
    integ = [t for t in trip if qs[t].denominator == 1]
    frac = [t for t in trip if qs[t].denominator != 1]
    print(f"  por especie: {len(integ)} tuplas de carga inteira, {len(frac)} de carga fracionaria")
    for t in integ:
        print(f"    {t}: Q = {qs[t]}")
    used_L = {(0, 1, 2), (1, 2, 3), (0, 2, 4), (0, 1, 5)}
    newL = [(t, qs[t]) for t in integ if t not in used_L]
    print(f"  setor L, carga inteira alem de e e dos tres nu: {[(t, str(q)) for t, q in newL]}")
    print(f"  -> leptons de MATERIA com Q > 0 previstos: {[(t, str(q)) for t, q in newL if q > 0]} (nunca vistos: precisam ser pesados)")
    used_B = {(0, 4, 5), (0, 1, 5), (1, 2, 3), (0, 2, 4)}
    print(f"  setor B, carga inteira alem de p, n, Lambda, Sigma0: {[(t, str(qs[t])) for t in integ if t not in used_B]}")
    print(f"  em cada setor, {len(frac)} tuplas de carga fracionaria (ja presentes na PH original; dobram)")

    print("\n##### Parte 3: custo nas previsoes numericas #####")
    K, GN, Nc, ZNP, Ne, g2 = 1.16335, 6.70883e-39, 6, F(32, 25), 10, 2 * math.pi ** 2 / 3
    def preds(ND, Nf):
        b0 = (F(ND) + F(2, Nf)) / 2
        be = b0 / (ZNP * Nc)
        S = 6 * Ne / float(be)
        Lam = math.sqrt(K / GN) * math.exp(-S / (2 * g2))
        me = Lam * math.exp(-Ne / float(be))
        rho = Nc / (27 * math.pi ** 4) * Lam ** 4 * math.exp(-S)
        return be, Lam, me, rho
    be0, Lam0, me0, rho0 = preds(4, 6)
    print(f"  base (N_Dirac = 4, N_f = 6): b_e = {be0} = {float(be0):.4f}, Lambda = {Lam0:.4e} GeV, "
          f"m_e = {me0*1e3:.4f} MeV, rho_Lambda = {rho0:.3e} GeV^4")
    for lab, ND, Nf in (("R0 indice fora da corrida", 4, 6), ("R1 N_f 6 -> 12", 4, 12), ("R2 N_Dirac 4 -> 8", 8, 6)):
        be, Lam, me, rho = preds(ND, Nf)
        print(f"  {lab:26s}: b_e = {float(be):.4f}  m_e = {me*1e3:.4e} MeV (x{me/me0:.3e})  "
              f"rho_Lambda = {rho:.3e} (x{rho/rho0:.2e})")

    print("\n##### Parte 4: bariogenese com mistura de especies (A1') #####")
    MPL, gs, Lam, mN = 1.2209e19, 106.75, 1.2387e12, 0.939
    HBAR, YR = 6.582e-25, 3.156e7
    tau_c1 = HBAR / (mN ** 5 / Lam ** 4) / YR
    cmax = math.sqrt(tau_c1 / 3.9e32)
    Tstar = lambda c: (1.66 * math.sqrt(gs) * Lam ** 4 / (c ** 2 * MPL)) ** (1 / 3)
    print(f"  c = 1: tau(p -> nu pi+) = {tau_c1:.1e} anos; limite 3.9e32 -> c <= {cmax:.1e}")
    print(f"  T* (mistura ativa acima de T*): c = 1 -> {Tstar(1):.2e} GeV; c = c_max -> {Tstar(cmax):.2e} GeV; "
          f"Lambda_meta = {Lam:.2e} GeV")
    print(f"  janela de violacao de B abaixo de Lambda_meta com c permitido: {'SIM' if Tstar(cmax) < Lam else 'NAO'}")
