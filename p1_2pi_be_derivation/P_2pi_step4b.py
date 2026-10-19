#!/usr/bin/env python3
"""
P1 / Passo 4b: D3' (ver P1_STEP3_D3PRIME_DEFINITION.md).
Resolve a equação de gap NÃO LINEAR em g = g_c (1+delta) e mede
R = W[s*]/Tr(s*^2) -> lambda_max.  Usa make_W e hermitian_basis de P_2pi_step4.py.
W_crit' = n_livres * lambda_max ; Z' = 2^n / W_crit' ; be' = b0/(Z' Nc).
Nada usa 325/1152, 32/25, m_e ou rho_Lambda antes da impressão final.
"""
import math
import numpy as np
from P_2pi_step4 import make_W, hermitian_basis, spin_mats, critical

EULER = 0.5772156649015329


def E1(x):
    x = max(x, 1e-300)
    if x <= 5:
        s = 0.0
        term = 1.0
        for k in range(1, 80):
            term *= -x / k
            s += -term / k          # -sum (-x)^k/(k k!)
        return -EULER - math.log(x) + s
    t = 1.0
    s = 1.0
    for k in range(1, 12):
        t *= -k / x
        s += t
    return math.exp(-x) / x * s


def r_cut(x):
    return 1.0 if x < 1e-300 else 1.0 - x * math.log1p(1.0 / x)


def r_pt(x):
    return 1.0 if x < 1e-300 else math.exp(-x) - x * E1(x)


REG = {"corte": r_cut, "tempo-próprio": r_pt}


def build(j2, ND, sigma):
    n, j, _ = spin_mats(j2)
    W, _, _ = make_W(j2, ND, sigma)
    support = list(range(1, n))
    basis, nr = hermitian_basis(n, support)
    B = np.array([E / math.sqrt(nr_i) for E, nr_i in zip(basis, nr)])  # ortonormal
    d = len(B)
    A = np.zeros((d, d))
    W1 = [W(b) for b in B]
    for i in range(d):
        A[i, i] = W1[i]
        for k in range(i + 1, d):
            A[i, k] = A[k, i] = 0.5 * (W(B[i] + B[k]) - W1[i] - W1[k])
    return n, j, W, support, B, A


def s_of(x, B, rf):
    M = np.tensordot(x, B, axes=1)
    M = (M + M.conj().T) / 2
    w, U = np.linalg.eigh(M)
    fw = np.array([wi * rf(wi * wi) for wi in w])
    s = (U * fw) @ U.conj().T
    return np.array([np.trace(b @ s).real for b in B])


def solve(A, B, g, x0, rf, maxit=300):
    x = x0.copy()

    def res(x):
        return x - g * (A @ s_of(x, B, rf))
    for it in range(maxit):
        F = res(x)
        nx = np.linalg.norm(x)
        if np.linalg.norm(F) < 1e-13 * nx:
            return x, it, True
        h = 1e-6 * nx
        d = len(x)
        J = np.zeros((d, d))
        for i in range(d):
            e = np.zeros(d)
            e[i] = h
            J[:, i] = (res(x + e) - res(x - e)) / (2 * h)
        dx = np.linalg.lstsq(J, -F, rcond=1e-10)[0]
        t = 1.0
        while t > 1e-4:
            xn = x + t * dx
            if np.linalg.norm(xn) > 0.3 * nx and np.linalg.norm(res(xn)) < np.linalg.norm(F):
                break
            t /= 2
        x = x + t * dx
    return x, maxit, np.linalg.norm(res(x)) < 1e-11 * np.linalg.norm(x)


def start(A, B, v, lam, delta, rf):
    lo, hi = 1e-9, 1.0
    def rho(mu):
        return (v @ s_of(mu * v, B, rf)) / mu
    if (1 + delta) * rho(hi) > 1:
        return None
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if (1 + delta) * rho(mid) > 1:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi) * v


def observe(A, B, x, rf, lam, top):
    xs = s_of(x, B, rf)
    R = (xs @ A @ xs) / (xs @ xs)
    ov = np.linalg.norm(top.T @ x) / np.linalg.norm(x)
    return R, ov, np.linalg.norm(x), xs @ A @ xs


def analyse(j2, ND, sigma, deltas, verbose=True):
    n, j, W, support, B, A = build(j2, ND, sigma)
    w, V = np.linalg.eigh(A)
    lam = w[-1]
    top = V[:, np.abs(w - lam) < 1e-9]
    v = top[:, 0]
    out = dict(n=n - 1, lam=lam, deg=top.shape[1], rows=[])
    if lam <= 0:
        return out
    gc = 1.0 / lam
    for name, rf in REG.items():
        for delta in deltas:
            x0 = start(A, B, v, lam, delta, rf)
            if x0 is None:
                out["rows"].append((name, delta, None))
                continue
            x, it, ok = solve(A, B, gc * (1 + delta), x0, rf)
            R, ov, mu, Wv = observe(A, B, x, rf, lam, top)
            out["rows"].append((name, delta, (R, ov, mu, Wv, it, ok)))
    out["A"], out["B"], out["v"], out["top"], out["gc"] = A, B, v, top, gc
    return out


def phi2_second_derivative(out, delta, rf):
    A, B, v, top, gc, lam = out["A"], out["B"], out["v"], out["top"], out["gc"], out["lam"]
    h = delta / 2
    vals = []
    for k in (-1, 0, 1):
        lng = math.log(gc * (1 + delta)) + k * h
        g = math.exp(lng)
        x0 = start(A, B, v, lam, g / gc - 1, rf)
        x, it, ok = solve(A, B, g, x0, rf)
        xs = s_of(x, B, rf)
        vals.append(0.5 * g * (xs @ A @ xs))
    return (vals[0] - 2 * vals[1] + vals[2]) / h ** 2, vals[1]


def main():
    print("##### PASSO 4b — D3' (equação de gap não linear) #####")
    print("W_crit' = n_livres * lambda_max ;  Z' = 2^n/W_crit' ;  b_e' = b0/(Z' Nc)\n")
    deltas = [1e-1, 1e-2, 1e-3, 1e-4]
    summary = {}
    for j2 in (1, 3, 5, 7, 9):
        n, j, _ = spin_mats(j2)
        nl = n - 1
        b0 = (4 * j + 3) / (2 * j + 1)
        Nc = 2 * j + 1
        print(f"=== j = {j2}/2  (n_livres = {nl}) ===")
        for sigma in (+1, -1):
            for ND in (2, 4):
                out = analyse(j2, ND, sigma, deltas if j2 in (3, 5) else [1e-2])
                lam = out["lam"]
                Wp = nl * lam
                Z = 2 ** nl / Wp
                be = b0 / (Z * Nc)
                W, nn, _ = make_W(j2, ND, sigma)
                Wc3, info = critical(W, nn, list(range(1, nn)))
                Z3 = (2 ** nl / Wc3) if Wc3 else float("nan")
                print(f"  s={sigma:+d} N_D={ND}: lam_max={lam:.5f} (deg {out['deg']})  "
                      f"W'={Wp:.4f}  Z'={Z:.4f}  b_e'={be:.5f}   | D3: Z={Z3:.4f}")
                summary[(j2, sigma, ND)] = (lam, Z, be)
                for name, delta, res in out["rows"]:
                    if res is None:
                        print(f"      [{name}] delta={delta:g}: sem solução")
                        continue
                    R, ov, mu, Wv, it, ok = res
                    print(f"      [{name}] delta={delta:g}: |M|={mu:.3e}  R={R:.6f}  "
                          f"R/lam-1={R/lam-1:+.2e}  overlap={ov:.6f}  it={it}  conv={ok}")
        print()

    print("##### d2 Phi2 / d(ln g)^2 em j=5/2 (deve -> 0) #####")
    for sigma in (+1, -1):
        for ND in (2, 4):
            out = analyse(5, ND, sigma, [1e-2])
            for delta in (1e-1, 1e-2, 1e-3):
                d2, p0 = phi2_second_derivative(out, delta, r_cut)
                print(f"  s={sigma:+d} N_D={ND} delta={delta:g}: Phi2={p0:.3e}  d2Phi2={d2:.3e}")

    print("\n##### COMPARAÇÃO FINAL (só agora), j = 5/2 #####")
    print("  alvo: Z_NP = 32/25 = 1.2800 ; b_e = 325/1152 = 0.28212")
    for sigma in (+1, -1):
        for ND in (2, 4):
            lam, Z, be = summary[(5, sigma, ND)]
            print(f"  s={sigma:+d} N_D={ND}: Z'={Z:.4f} ({100*(Z/1.28-1):+.1f}%)  "
                  f"b_e'={be:.5f} ({100*(be/(325/1152)-1):+.1f}%)")


if __name__ == "__main__":
    main()
