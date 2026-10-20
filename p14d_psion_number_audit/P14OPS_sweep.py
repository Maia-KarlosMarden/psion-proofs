#!/usr/bin/env python3
"""P14 / WP-OPS: termos de CP simples feitos de operadores da PH (ver P14OPS_DEFINITION.md)."""
import io, contextlib, itertools, json, runpy
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    C = runpy.run_path('P14CAL_fit.py', run_name='lib')
D_ = C['D_']
S = D_['S']
H0t, Ot = S['H0t'], S['Ot']
TR, NT, EDGES, REARR = C['TR'], C['NT'], C['EDGES'], D_['REARR']
NBE, NBP, nE, NV, Wv = C['NBE'], C['NBP'], C['nE'], C['NV'], C['Wv']
response, constraints, pattern, DT = C['response'], C['constraints'], C['pattern'], C['DT']
Q9 = [1, -5, -5, 10, 4, 4]

diag = {
    "Q": np.array([sum(Q9[k] for k in t) / 9 for t in TR]),
    "M": np.array([sum(k - 2.5 for k in t) for t in TR]),
    "NR": np.array([sum(1 for k in t if k >= 3) for t in TR], float),
    "Z": np.array([sum(t) % 6 for t in TR], float),
    "E": (H0t + Ot).diagonal().copy(),
}
offd = {"h": H0t, "O": Ot}

cands = {}
for a, A in diag.items():
    for b, B in offd.items():
        cands[f"[{a},{b}]"] = ("x", np.diag(A) @ B - B @ np.diag(A))
cands["[h,O]"] = ("x", H0t @ Ot - Ot @ H0t)
for a, f in diag.items():
    cands[f"rearr {a}{a}"] = ("c", np.array([f[i] * f[j] - f[k] * f[l] for i, j, k, l in REARR]))


def to_v(kind, M):
    v = np.zeros(NV)
    if kind == "x":
        y = np.array([M[i, j] for i, j in EDGES])
        pv = NBE.T @ y
        frac = np.linalg.norm(NBE @ pv) / max(np.linalg.norm(y), 1e-300)
        v[:nE] = pv
    else:
        pv = NBP.T @ M
        frac = np.linalg.norm(NBP @ pv) / max(np.linalg.norm(M), 1e-300)
        v[nE:] = pv
    n = np.linalg.norm(v)
    return (v / n if n > 0 else v), frac


def fit(A, U, b):
    """min ||A U e - b||; devolve rho, e, v."""
    M = A @ U
    e, *_ = np.linalg.lstsq(M, b, rcond=None)
    return np.linalg.norm(M @ e - b) / np.linalg.norm(b), e, U @ e


if __name__ == "__main__":
    out = {}
    for law in ("L2", "L1"):
        G = response(1e-3, law)
        A, b0 = constraints(G)
        b = b0 * DT
        print(f"\n##### c_G = 1e-3, lei {law} #####")
        V = {}
        for name, (kind, M) in cands.items():
            u, frac = to_v(kind, M)
            if np.linalg.norm(u) == 0:
                print(f"  {name:12s}: anulado pela projecao (fracao permitida 0)")
                continue
            V[name] = u
            rho, e, v = fit(A, u[:, None], b)
            out[f"{law}_{name}"] = dict(rho=float(rho), frac=float(frac), eps_ef=float(pattern(v)[2]))
            print(f"  {name:12s}: fracao permitida {frac:5.2f}  rho = {rho:.4f}  eps_ef = {pattern(v)[2]:.2e}")
        names = list(V)
        best2 = min(((fit(A, np.stack([V[p], V[q]], 1), b)[0], p, q) for p, q in itertools.combinations(names, 2)))
        r2, p, q = best2
        v2 = fit(A, np.stack([V[p], V[q]], 1), b)[2]
        best1 = min(out[f"{law}_{n}"]["rho"] for n in names)
        print(f"  melhor 1 termo: rho = {best1:.4f}; melhor par: {p} + {q}, rho = {r2:.4f}, eps_ef = {pattern(v2)[2]:.2e}")

        rng = np.random.default_rng(0)
        R1 = [fit(A, (lambda u: u / np.linalg.norm(u))(rng.standard_normal(NV))[:, None], b)[0] for _ in range(2000)]
        R2 = [fit(A, rng.standard_normal((NV, 2)), b)[0] for _ in range(2000)]
        p1, p2 = np.percentile(R1, [1, 50]), np.percentile(R2, [1, 50])
        print(f"  linha de base 1 termo: 1% = {p1[0]:.4f}, mediana = {p1[1]:.4f}, minimo = {min(R1):.4f}")
        print(f"  linha de base 2 termos: 1% = {p2[0]:.4f}, mediana = {p2[1]:.4f}, minimo = {min(R2):.4f}")
        succ = min(best1, r2) < 0.1
        struct = best1 < p1[0] or r2 < p2[0]
        print(f"  CRITERIO: sucesso {'SIM' if succ else 'nao'}; estrutura acima do acaso {'SIM' if struct else 'nao'}")
        out[f"{law}_summary"] = dict(best1=float(best1), best2=float(r2), pair=[p, q],
                                     base1=[float(x) for x in p1], base2=[float(x) for x in p2],
                                     success=bool(succ), structure=bool(struct))
    json.dump(out, open("P14OPS_results.json", "w"), indent=1)
