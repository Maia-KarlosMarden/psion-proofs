#!/usr/bin/env python3
"""P14 / WP-CAL: calibracao do universo PH (ver P14CAL_DEFINITION.md)."""
import io, contextlib, json, math, runpy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.optimize import linprog

with contextlib.redirect_stdout(io.StringIO()):
    D_ = runpy.run_path('P14DOM_sim.py', run_name='lib')
    T_ = runpy.run_path('P14TPROG_sim.py', run_name='lib')
TR, NT, EDGES, W0 = D_['TR'], D_['NT'], D_['EDGES'], D_['W0']
BE, BP, make_rhs = D_['BE'], D_['BP'], D_['make_rhs']
gamma_over_H, T_RH, T_END = D_['gamma_over_H'], D_['T_RH'], D_['T_END']
TI = {t: i for i, t in enumerate(TR)}
IP, IE, IN1, IN2 = TI[(0, 4, 5)], TI[(0, 1, 2)], TI[(1, 2, 3)], TI[(0, 2, 4)]
NBE, NBP = null_space(BE), null_space(BP)
nE, nP = NBE.shape[1], NBP.shape[1]
NV = nE + nP
Wv = np.array([W0[i, j] for i, j in EDGES])
Y0 = 0.208 * 12 / 106.75
ETA = 6.1e-10
DT = ETA / (7.04 * Y0)


def pattern(v):
    y = NBE @ v[:nE]
    c = NBP @ v[nE:]
    X = np.zeros((NT, NT))
    for e, (i, j) in enumerate(EDGES):
        X[i, j], X[j, i] = y[e], -y[e]
    return X, c, max(np.abs(y / Wv).max(), np.abs(c).max())


def run(v, cG, law, cGmul=1.0):
    X, c, _ = pattern(v)
    f = make_rhs(law, 1.0, X, c, lambda x: gamma_over_H(T_RH * math.exp(-x), cG * cGmul))
    Y0v = np.concatenate([np.full(12, 1 / 12), np.zeros(2 * NT)])
    sol = solve_ivp(f, [0, math.log(T_RH / T_END)], Y0v, method="LSODA", rtol=1e-11, atol=1e-22)
    Y = sol.y[:, -1]
    return np.concatenate([Y[12:12 + NT] - Y[12 + NT:], Y[:6] - Y[6:12]])   # 20 tuplas + 6 livres


def response(cG, law):
    G = np.zeros((26, NV))
    for k in range(NV):
        e = np.zeros(NV); e[k] = 1.0
        s = 1e-3 / pattern(e)[2]
        G[:, k] = run(s * e, cG, law) / s
    return G


def constraints(G):
    rows, b = [G[IP], G[IE], G[IN1] + G[IN2]], [1.0, 1.0, -2.0]
    for t in range(NT):
        if t not in (IP, IE, IN1, IN2):
            rows.append(G[t]); b.append(0.0)
    for k in range(5):
        rows.append(G[20 + k]); b.append(0.0)
    return np.array(rows), np.array(b)


def solve_lp(G):
    Aeq, beq = constraints(G)
    Mx = (NBE / Wv[:, None])
    nx = Mx.shape[0]
    Aub = []
    for M, off in ((Mx, 0), (NBP, nE)):
        for r in M:
            row = np.zeros(NV + 1); row[off:off + M.shape[1]] = r; row[-1] = -1
            Aub.append(row.copy()); row[:NV] *= -1; Aub.append(row)
    Aub = np.array(Aub)
    sc = 1 / np.abs(Aeq).max()
    res = linprog(np.r_[np.zeros(NV), 1.0], A_ub=Aub, b_ub=np.zeros(len(Aub)),
                  A_eq=np.c_[Aeq * sc, np.zeros(len(Aeq))], b_eq=beq * sc,
                  bounds=[(None, None)] * NV + [(0, None)], method="highs")
    return res


def eta_of(d):
    return 7.04 * Y0 * d


if __name__ == "__main__":
    print(f"alvo: eta = {ETA}, Delta_p = Delta_e = D = {DT:.4e}; dimensao dos botoes de direcao = {NV}")
    best, out = None, {}
    for cG in (1e-5, 1e-3, 1.0):
        for law in ("L1", "L2"):
            G = response(cG, law)
            Aeq, _ = constraints(G)
            r = solve_lp(G)
            rank = np.linalg.matrix_rank(Aeq, tol=1e-6 * np.abs(Aeq).max())
            if r.status == 0:
                eps_ef = DT * r.x[-1]
                print(f"  cG={cG:6.0e} {law}: posto das restricoes {rank}/24, viavel, pitada minima eps_ef = {eps_ef:.3e}")
                out[f"{cG}_{law}"] = dict(eps_ef=eps_ef, rank=int(rank))
                if best is None or eps_ef < best[0]:
                    best = (eps_ef, cG, law, DT * r.x[:NV], G)
            else:
                print(f"  cG={cG:6.0e} {law}: posto {rank}/24, INVIAVEL ({r.message})")
                out[f"{cG}_{law}"] = dict(eps_ef=None, rank=int(rank))
    eps_ef, cG, law, v, G = best
    print(f"\nRECEITA: c_G = {cG:.0e}, lei {law}, eps_ef = {eps_ef:.3e} "
          f"(V4 positividade: {'OK' if eps_ef < 1 else 'FALHA'})")

    print("\n##### V1 corrida nao linear com a receita #####")
    d = run(v, cG, law)
    others = [abs(d[t]) for t in range(NT) if t not in (IP, IE, IN1, IN2)] + list(abs(d[20:]))
    print(f"  eta(p) = {eta_of(d[IP]):.4e}  Delta_e/Delta_p = {d[IE]/d[IP]:.6f}  "
          f"(Delta_nu)/(Delta_p) = {(d[IN1]+d[IN2])/d[IP]:.6f}  maior |exotico|/D = {max(others)/DT:.2e}")
    v1 = abs(eta_of(d[IP]) / ETA - 1) < 0.01 and abs(d[IE] / d[IP] - 1) < 0.01 and max(others) / DT < 0.01
    print(f"  V1 {'OK' if v1 else 'FALHA'}")
    print(f"  N_psi liquido = {d[20:].sum() + 3 * d[:20].sum():.1e}")

    print("\n##### V2 sintonia fina #####")
    for m in (0.5, 2.0):
        dm = run(v, cG, law, cGmul=m)
        print(f"  c_G x {m}: eta = {eta_of(dm[IP]):.3e} (x{eta_of(dm[IP])/ETA:.3f}), "
              f"Delta_e/Delta_p = {dm[IE]/dm[IP]:.3f}, maior |exotico|/D = "
              f"{max(abs(dm[t]) for t in range(NT) if t not in (IP, IE, IN1, IN2))/DT:.2e}")
    rng = np.random.default_rng(0)
    for lvl in (0.01, 0.10):
        etas, ratios, exo = [], [], []
        for _ in range(10):
            vn = v * (1 + lvl * rng.standard_normal(NV))
            dn = G @ vn                       # linear (verificado em V1)
            etas.append(eta_of(dn[IP]) / ETA); ratios.append(dn[IE] / dn[IP])
            exo.append(max(abs(dn[t]) for t in range(NT) if t not in (IP, IE, IN1, IN2)) / DT)
        print(f"  ruido {lvl:.0%} na direcao: eta/alvo em [{min(etas):.2f}, {max(etas):.2f}], "
              f"Delta_e/Delta_p em [{min(ratios):.2f}, {max(ratios):.2f}], |exotico|/D ate {max(exo):.2f}")

    print("\n##### V3 neutrinos #####")
    dnu = 2 * ETA / 2                          # por sabor, em n/n_gamma
    xi = dnu / 0.2488
    print(f"  excesso de nu-bar por sabor = {dnu:.1e} n_gamma -> xi = {xi:.1e}; limite |xi| < 0.03: OK; "
          f"mu_nu/T = 0.123 de 13.5 seria {0.123/xi:.0e} vezes maior")

    print("\n##### A5 / V5 materia escura do topo (historia H1) #####")
    base = T_['sectors']
    for kap in (1.0, 4 * math.pi):
        def omega(l, kap=kap):
            m = T_['m_lev'](l)
            T_['sectors'] = lambda m8=0.938: [dict(name="topo", m=m, g=2, sv=kap / m ** 2)] + base(m8)[1:]
            r = T_['simulate'](max(1e3 * m, 1e3), [("cool", 1e-5)], {}, dl=0.004)
            return r['topo']['Oh2'], m
        lo, hi = 0.0, 12.0
        for _ in range(40):
            mid = (lo + hi) / 2
            if omega(mid)[0] > 0.12:
                lo = mid
            else:
                hi = mid
        oh, m = omega(hi)
        print(f"  kappa = {kap:.2f}: l_DM = {hi:.3f}, m_DM = {m:.3e} GeV, Omega h2 = {oh:.3f}; "
              f"nivel 0,8 do v85 daria {omega(0.8)[0]:.1e}")
    T_['sectors'] = base
    json.dump(dict(out=out, recipe=dict(cG=cG, law=law, eps_ef=eps_ef, v=list(v))), open("P14CAL_results.json", "w"), indent=1)
