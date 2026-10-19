#!/usr/bin/env python3
"""P14 / WP-A passo 2: bounce O(3) no potencial CW + massa termica; xi* = T_c/T_nuc."""
import math
import numpy as np
from scipy import integrate, optimize

V0 = 1.2609e12          # v = Lambda_meta (GeV)
MP = 2.435e18           # Planck reduzida
GSTAR = 106.75
LAMS = {"1/(2j)^2": 0.04, "eps_psi": 25 / 169, "1/Nc": 1 / 6, "N_c g^2/(4pi)^2": 0.25, "b_e": 325 / 1152}


def make_V(lam, c, T, v=V0):
    def V(p):
        if p <= 0:
            return 0.0
        return lam / 4 * p ** 4 * (math.log(p * p / (v * v)) - 0.5) + 0.5 * c * T * T * p * p
    def dV(p):
        if p <= 0:
            return 0.0
        return lam * p ** 3 * math.log(p * p / (v * v)) + c * T * T * p
    return V, dV


def phi_min(V, lo, hi):
    r = optimize.minimize_scalar(V, bounds=(lo, hi), method="bounded", options={"xatol": lo * 1e-6})
    return r.x, r.fun


def Tc_of(lam, c):
    def f(lnT):
        T = math.exp(lnT)
        V, _ = make_V(lam, c, T)
        return phi_min(V, 1e-4 * V0, 3 * V0)[1]
    return math.exp(optimize.brentq(f, math.log(1e-6 * V0), math.log(3 * V0)))


def bounce_S3(lam, c, T):
    """S_3 via shooting em unidades v=1, V/v^4."""
    V, dV = make_V(lam, c, T, v=1.0)
    pmin, _ = phi_min(V, 1e-4, 3.0)
    # barreira: V(p)>0 para p pequeno (massa termica), V(pmin)<0
    def rhs(r, y):
        return [y[1], dV(y[0]) - 2 * y[1] / r if r > 0 else dV(y[0]) / 3]
    def shoot(p0):
        sol = integrate.solve_ivp(rhs, [1e-6, 1e4], [p0, 0.0], rtol=1e-7, atol=1e-12,
                                  events=[lambda r, y: y[0], lambda r, y: y[1]], dense_output=True, max_step=2.0)
        # undershoot: phi' hits 0 while phi>0 (returns)-> event 1 ; overshoot: phi crosses 0 -> event 0
        t0 = sol.t_events[0]; t1 = sol.t_events[1]
        t1 = t1[t1 > 1e-3]
        if len(t0) and (not len(t1) or t0[0] < t1[0]):
            return +1, sol
        return -1, sol
    lo, hi = 1e-3 * pmin, pmin * 0.999999
    # garantir lo undershoot, hi overshoot
    s_lo, _ = shoot(lo); s_hi, _ = shoot(hi)
    if not (s_lo == -1 and s_hi == +1):
        return float("nan")
    for _ in range(45):
        mid = math.sqrt(lo * hi)
        s, _ = shoot(mid)
        if s == +1:
            hi = mid
        else:
            lo = mid
    s, sol = shoot(math.sqrt(lo * hi))
    r_end = sol.t_events[0][0] if len(sol.t_events[0]) else sol.t[-1]
    rr = np.linspace(1e-6, r_end, 4000)
    y = sol.sol(rr)
    integrand = 4 * math.pi * rr ** 2 * (0.5 * y[1] ** 2 + np.array([V(p) for p in y[0]]))
    return float(np.trapezoid(integrand, rr))  # em unidades v^4 * (1/v)^3 = v


def run(lam, c):
    Tc = Tc_of(lam, c)
    dV = lam * V0 ** 4 / 8
    Teq = (30 * dV / (math.pi ** 2 * GSTAR)) ** 0.25
    TRH = Teq
    def lnGamma_over_H4(lnT):
        T = math.exp(lnT)
        S3 = bounce_S3(lam, c, T / V0) * V0   # GeV
        if math.isnan(S3):
            return float("nan")
        S3T = S3 / T
        H2 = (math.pi ** 2 / 30 * GSTAR * T ** 4 + dV) / (3 * MP ** 2)
        return 4 * math.log(T) + 1.5 * math.log(S3T / (2 * math.pi)) - S3T - 2 * math.log(H2)
    # varrer T abaixo de T_c
    grid = np.exp(np.linspace(math.log(0.98 * Tc), math.log(1e-8 * Tc), 28))
    vals = []
    for T in grid:
        g = lnGamma_over_H4(math.log(T))
        vals.append(g)
        if not math.isnan(g) and g > 0:
            break
    vals = np.array(vals)
    ok = np.where(~np.isnan(vals) & (vals > 0))[0]
    if len(ok) == 0:
        return dict(Tc=Tc, Teq=Teq, Tnuc=None, xi=None, S3T=None, note="sem nucleacao ate 1e-8 T_c (inflacao eterna no modelo)")
    i = ok[0]
    a, b = math.log(grid[i - 1]), math.log(grid[i])
    lnTn = optimize.brentq(lnGamma_over_H4, a, b, xtol=2e-2)
    Tn = math.exp(lnTn)
    S3T = bounce_S3(lam, c, Tn / V0) * V0 / Tn
    return dict(Tc=Tc, Teq=Teq, Tnuc=Tn, xi=Tc / Tn, S3T=S3T,
                note=("INFLACAO TERMICA: T_nuc < T_eq, reaquece a T_RH~%.1e" % TRH) if Tn < Teq else "sem dominacao de vacuo")


if __name__ == "__main__":
    print("##### WP-A passo 2: CW + massa termica, bounce O(3), Gamma = H^4 #####")
    print(f"v = Lambda_meta = {V0:.4e} GeV ; g* = {GSTAR} ; M_P(red) = {MP:.3e}")
    print(f"{'lambda':>18s} {'c':>4s} {'T_c [GeV]':>10s} {'T_eq':>9s} {'T_nuc':>9s} {'S3/T':>6s} {'xi*':>9s}  nota")
    rows = []
    for name, lam in LAMS.items():
        for c in (1.0, 0.3, 3.0):
            try:
                r = run(lam, c)
            except Exception as e:
                r = dict(Tc=float("nan"), Teq=float("nan"), Tnuc=None, xi=None, S3T=None, note=f"erro: {e}")
            xi = f"{r['xi']:9.1f}" if r["xi"] else "      ---"
            tn = f"{r['Tnuc']:9.2e}" if r["Tnuc"] else "      ---"
            s3 = f"{r['S3T']:6.1f}" if r["S3T"] else "   ---"
            print(f"{name:>18s} {c:4.1f} {r['Tc']:10.3e} {r['Teq']:9.2e} {tn} {s3} {xi}  {r['note']}")
            rows.append((name, lam, c, r))
    print("\nfaixa-alvo de WP-C: 44 <= xi* <= 302 e T_nuc > T_eq")
    for name, lam, c, r in rows:
        if r["xi"] and 44 <= r["xi"] <= 302 and r["Tnuc"] > r["Teq"]:
            print(f"  COMPATIVEL: lambda={name} ({lam:.4f}), c={c}")
