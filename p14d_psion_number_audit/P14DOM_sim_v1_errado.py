#!/usr/bin/env python3
"""P14 / WP-DOM: montagem num dominio Z6 que viola CP (ver P14DOM_DEFINITION.md)."""
import io, contextlib, json, math, runpy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.sparse.csgraph import connected_components

with contextlib.redirect_stdout(io.StringIO()):
    S = runpy.run_path('P14NOW_sim.py', run_name='lib')
TR, NT, REARR, weights = S['TR'], S['NT'], S['REARR'], S['weights']
gamma_over_H, T_RH, T_END, K_EQ = S['gamma_over_H'], S['T_RH'], S['T_END'], S['K_EQ']
IP, IE = S['IP'], S['IE']
idx = np.array(TR)
ra, rb, rc, rd = REARR.T
NR = len(REARR)
RTOL = 1e-11
YEQ, YEQT = 1 / 12, K_EQ * (1 / 12) ** 3

_, W0 = weights(1.0, "L2", 0.0, 1)          # W nao depende da lei
EDGES = [(i, j) for i in range(NT) for j in range(i + 1, NT) if W0[i, j] > 1e-14]
NE = len(EDGES)
BE = np.zeros((NT, NE))
for e, (i, j) in enumerate(EDGES):
    BE[i, e], BE[j, e] = -1, 1
BR = np.zeros((NT, NR))
for r, (a, b, c, d) in enumerate(REARR):
    BR[a, r] -= 1; BR[b, r] -= 1; BR[c, r] += 1; BR[d, r] += 1


def cp_pattern(variant, seed):
    rng = np.random.default_rng(seed)
    if variant == "V1":
        N = null_space(BE)
        y = N @ rng.standard_normal(N.shape[1]) if N.shape[1] else np.zeros(NE)
        cr = np.zeros(NR)
    else:
        N = null_space(np.hstack([BE, BR]))
        v = N @ rng.standard_normal(N.shape[1])
        y, cr = v[:NE], v[NE:]
    x = np.array([y[e] / W0[i, j] for e, (i, j) in enumerate(EDGES)])
    s = max(np.abs(x).max(initial=0), np.abs(cr).max(initial=0))
    if s > 0:
        x, cr = x / s, cr / s
    X = np.zeros((NT, NT))
    for e, (i, j) in enumerate(EDGES):
        X[i, j], X[j, i] = W0[i, j] * x[e], -W0[i, j] * x[e]
    return X, cr


def make_rhs(law, eps, X, cr, rate):
    kb, W = weights(1.0, law, 0.0, 1)
    sec = []
    for s in (+1, -1):
        R = W + s * eps * X
        assert (R >= -1e-15).all()
        sec.append((R, 1 + s * eps * cr, 1 - s * eps * cr))

    def f(x, Y):
        r = rate(x)
        psi, apsi = Y[:6], Y[6:12]
        t, at = Y[12:12 + NT], Y[12 + NT:]
        out = []
        dd = []
        for (P, Tt), (R, fw, bw) in zip(((psi, t), (apsi, at)), sec):
            dP, dT = np.zeros(6), np.zeros(NT)
            J = kb * (P[idx[:, 0]] * P[idx[:, 1]] * P[idx[:, 2]] - Tt / K_EQ)
            dT += J
            for col in range(3):
                np.subtract.at(dP, idx[:, col], J)
            dT += R.T @ Tt - R.sum(axis=1) * Tt
            Jr = fw * Tt[ra] * Tt[rb] - bw * Tt[rc] * Tt[rd]
            np.subtract.at(dT, ra, Jr); np.subtract.at(dT, rb, Jr)
            np.add.at(dT, rc, Jr); np.add.at(dT, rd, Jr)
            dd.append((dP, dT))
        (dpsi, dt), (dapsi, dat) = dd
        Ja = psi * apsi - YEQ ** 2
        dpsi -= Ja; dapsi -= Ja
        Jt = t * at - YEQT ** 2
        dt -= Jt; dat -= Jt
        return r * np.concatenate([dpsi, dapsi, dt, dat])
    return f


def summarize(Y):
    psi, apsi = Y[:6], Y[6:12]
    t, at = Y[12:12 + NT], Y[12 + NT:]
    return dict(dp=t[IP] - at[IP], de=t[IE] - at[IE], dfree=(psi - apsi).sum(),
                dtup=(t - at).sum(), Nnet=(psi - apsi).sum() + 3 * (t - at).sum())


def runA(variant, cG, law, eps, seed):
    X, cr = cp_pattern(variant, seed)
    f = make_rhs(law, eps, X, cr, lambda x: gamma_over_H(T_RH * math.exp(-x), cG))
    Y0 = np.concatenate([np.full(12, 1 / 12), np.zeros(2 * NT)])
    sol = solve_ivp(f, [0, math.log(T_RH / T_END)], Y0, method="LSODA", rtol=RTOL, atol=1e-22)
    return summarize(sol.y[:, -1])


if __name__ == "__main__":
    nc, lab = connected_components((W0 > 1e-14).astype(int))
    print(f"grafo de transicoes: {NE} arestas, {nc} componentes; "
          f"componente do eletron: {[TR[i] for i in range(NT) if lab[i] == lab[IE]]}")
    print(f"dim nucleo V1 (circulacoes de transicao) = {null_space(BE).shape[1]}, "
          f"V2 (conjunto) = {null_space(np.hstack([BE, BR])).shape[1]}")

    print("\n##### Controles #####")
    c1 = max(max(abs(v) for v in runA(V, cG, law, 0.0, 0).values())
             for V in ("V1", "V2") for cG in (1e-5, 1e-3, 1.0) for law in ("L1", "L2"))
    print(f"  C1 eps = 0: maior |assimetria| = {c1:.1e}  {'OK' if c1 < 1e-12 else 'FALHA'}")
    c4 = 0.0
    for V in ("V1", "V2"):
        for law in ("L1", "L2"):
            X, cr = cp_pattern(V, 0)
            f = make_rhs(law, 1e-3, X, cr, lambda x: 1e3)
            Yeq = np.concatenate([np.full(12, YEQ), np.full(2 * NT, YEQT)])
            sol = solve_ivp(f, [0, 20], Yeq, method="LSODA", rtol=RTOL, atol=1e-22)
            c4 = max(c4, max(abs(v) for v in summarize(sol.y[:, -1]).values()))
    print(f"  C4 equilibrio preservado (eps = 1e-3): maior |assimetria| = {c4:.1e}  {'OK' if c4 < 1e-12 else 'FALHA'}")

    print("\n##### Parte A: 120 corridas #####")
    res, c2, c3 = {}, 0.0, 0.0
    for V in ("V1", "V2"):
        for cG in (1e-5, 1e-3, 1.0):
            for law in ("L1", "L2"):
                for eps in (1e-6, 1e-3):
                    rows = []
                    for seed in range(5):
                        r = runA(V, cG, law, eps, seed)
                        c2 = max(c2, abs(r['Nnet']))
                        if eps == 1e-3 and seed == 0:
                            rm = runA(V, cG, law, -eps, seed)
                            for k in ("dp", "de", "dfree"):
                                if abs(r[k]) > 1e-12:
                                    c3 = max(c3, abs(rm[k] / r[k] + 1))
                        res[f"{V}_{cG}_{law}_{eps}_{seed}"] = {k: float(v) for k, v in r.items()}
                        rows.append(r)
                    dp = np.array([x['dp'] for x in rows]); de = np.array([x['de'] for x in rows])
                    df = np.array([x['dfree'] for x in rows])
                    hyd = sum(1 for x in rows if x['dp'] > 0 and x['de'] > 0 and 0.1 < x['de'] / x['dp'] < 10)
                    print(f"  {V} cG={cG:6.0e} {law} eps={eps:6.0e}: |Dp|/eps med={np.median(abs(dp))/eps:9.2e} "
                          f"[{abs(dp).min()/eps:8.1e},{abs(dp).max()/eps:8.1e}] sinais p {''.join('+' if v > 0 else '-' for v in dp)} | "
                          f"|De|/eps max={abs(de).max()/eps:8.1e} | Dlivre/eps med={np.median(df)/eps:+9.2e} | hidrogenio {hyd}/5")
    print(f"\n  C2 N_psi liquido: maior |N_psi| = {c2:.1e}  {'OK' if c2 < 1e-12 else 'FALHA'}")
    print(f"  C3 dominio conjugado: maior desvio de D(-eps) = -D(eps): {c3:.1e}  {'OK' if c3 < 1e-6 else 'FALHA'}")
    lin = max(abs(res[k.replace('_0.001_', '_1e-06_')]['dp'] / 1e-6 / (res[k]['dp'] / 1e-3) - 1)
              for k in res if '_0.001_' in k and abs(res[k]['dp']) > 1e-9)
    print(f"  linearidade em eps (Dp/eps entre 1e-6 e 1e-3): maior desvio {lin:.1e}")

    print("\n##### Parte B: epoca e tamanho dos dominios #####")
    MPL, T0, GPC = 1.2209e19, 2.348e-13, 3.086e27
    for Tdom, gs in ((3e-3, 10.75), (1.38e11, 106.75)):
        gh = gamma_over_H(Tdom, 1.0)
        t_s = 0.301 * gs ** -0.5 * MPL / Tdom ** 2 * 6.582e-25
        d0 = 2 * 2.998e10 * t_s * (Tdom / T0) * (gs / 3.91) ** (1 / 3)
        print(f"  T_dom = {Tdom:.2e} GeV: Gamma/H da montagem (c_G = 1) = {gh:.1e}; "
              f"tamanho comovel maximo hoje = {d0:.2e} cm = {d0/GPC:.1e} Gpc")
    print("  limite observacional para dominios de antimateria: >~ 1 Gpc (Cohen, De Rujula, Glashow 1998)")
    print("  media global (volumes iguais de eps e -eps, vacuos degenerados): Delta = 0 por C3")
    json.dump(res, open("P14DOM_results.json", "w"), indent=1)
