#!/usr/bin/env python3
"""P14 / WP-TPROG: programa de temperatura livre (ver P14TPROG_DEFINITION.md)."""
import json, math, sys
from scipy.special import kve

MPL, ZETA3 = 1.2209e19, 1.2020569
INV_B = math.log(610 / 0.511) / 2
LAM = 0.511e-3 * math.exp(10 * INV_B)
m_lev = lambda l: LAM * math.exp(-l * INV_B)
ALPHA = 1 / 137.036


def gstar(T):
    for cut, g in ((80, 106.75), (4, 86.25), (1, 75.75), (0.15, 61.75), (0.1, 17.25), (5e-4, 10.75)):
        if T > cut:
            return g
    return 3.91


def sectors(m8=0.938):
    mt = m_lev(0.8)
    return [dict(name="topo", m=mt, g=2, sv=1 / mt ** 2),
            dict(name="N8", m=m8, g=2, sv=1 / 0.135 ** 2),
            dict(name="e10", m=0.511e-3, g=2, sv=math.pi * ALPHA ** 2 / 0.511e-3 ** 2)]


def feq(x):
    if x < 1e-6:
        return 1.0
    if x > 700:
        return 0.0
    return 0.5 * x * x * kve(2, x) * math.exp(-x)


def Yeq(sec, T):
    s = 2 * math.pi ** 2 / 45 * gstar(T) * T ** 3
    return sec["g"] * T ** 3 / math.pi ** 2 * feq(sec["m"] / T) / s


def simulate(T0, prog, eta_in, dl=0.004, m8=0.938):
    secs = sectors(m8)
    T = T0
    st = []
    for sec in secs:
        A = eta_in.get(sec["name"], 0.0) / 7.04
        e = Yeq(sec, T)
        st.append([A, (-A + math.sqrt(A * A + 4 * e * e)) / 2])
    dil = 1.0
    for kind, Tt in prog:
        if kind == "reheat":
            if Tt > T:
                f = gstar(T) * T ** 3 / (gstar(Tt) * Tt ** 3)
                dil *= f
                for o in st:
                    o[0] *= f; o[1] *= f
                T = Tt
            continue
        while T > Tt:
            T = max(Tt, T * math.exp(-dl))
            H = 1.66 * math.sqrt(gstar(T)) * T * T / MPL
            s = 2 * math.pi ** 2 / 45 * gstar(T) * T ** 3
            for sec, o in zip(secs, st):
                k = s * sec["sv"] / H * dl
                e = Yeq(sec, T)
                b = 1 + k * o[0]
                c = o[1] + k * e * e
                o[1] = 2 * c / (b + math.sqrt(b * b + 4 * k * c))   # raiz positiva, forma estavel
    out = {}
    for sec, (A, yb) in zip(secs, st):
        out[sec["name"]] = dict(X=7.04 * (yb + A), Xbar=7.04 * yb, A=7.04 * A,
                                Oh2=2.74e8 * sec["m"] * (2 * yb + A))
    out["dil"] = dil
    return out


HIST = {
    "H1": [("cool", 1e-5)],
    "H2": [("cool", 0.01), ("reheat", 2.0), ("cool", 1e-5)],
    "H3": [("cool", 1e6), ("reheat", 1e11), ("cool", 1e-5)],
    "H4": [("cool", 0.05), ("reheat", 5.0), ("cool", 0.05), ("reheat", 5.0), ("cool", 1e-5)],
    "H5": [("cool", 1e-4), ("reheat", 0.05), ("cool", 1e-5)],
}

if __name__ == "__main__":
    print(f"massas: topo {m_lev(0.8):.3e} GeV, N8 0.938 GeV, e10 5.11e-4 GeV")
    res, c1, c2, c3 = {}, 0.0, 0.0, 0.0
    for h, prog in HIST.items():
        for lab, eta in (("eta0", {}), ("eta1e-9", {"N8": 1e-9, "e10": 1e-9})):
            r = simulate(1e12, prog, eta)
            rf = simulate(1e12, prog, eta, dl=0.001)
            res[f"{h}_{lab}"] = r
            if lab == "eta0":
                c1 = max(c1, max(abs(r[k]["A"]) for k in ("topo", "N8", "e10")))
            else:
                c2 = max(c2, abs(r["N8"]["A"] / (1e-9 * r["dil"]) - 1), abs(r["e10"]["A"] / (1e-9 * r["dil"]) - 1))
            for k in ("topo", "N8", "e10"):
                for q in ("X", "Xbar"):
                    if r[k][q] > 1e-300:
                        c3 = max(c3, abs(rf[k][q] / r[k][q] - 1))
            print(f"\n{h} {lab}  (diluicao total {r['dil']:.2e})")
            for k in ("topo", "N8", "e10"):
                x = r[k]
                extra = f"  Omega h2 = {x['Oh2']:.2e}" if k == "topo" else ""
                print(f"   {k:5s} n_X/n_g = {x['X']:.3e}  n_Xbar/n_g = {x['Xbar']:.3e}  assim = {x['A']:+.3e}{extra}")
    print(f"\nC1 eta_in = 0 -> maior |assimetria| = {c1:.1e}  {'OK' if c1 == 0 else 'FALHA'}")
    print(f"C2 assimetria = entrada x diluicao: maior desvio relativo = {c2:.1e}  {'OK' if c2 < 1e-12 else 'FALHA'}")
    print(f"C3 convergencia no passo (0.004 vs 0.001): maior diferenca = {c3:.1e}  {'OK' if c3 < 0.01 else 'FALHA'}")
    json.dump(res, open("P14TPROG_results.json", "w"), indent=1)
