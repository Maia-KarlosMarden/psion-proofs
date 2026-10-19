#!/usr/bin/env python3
"""
Teste cego ESTENDIDO de look-elsewhere para a PH v85.

(i)  Cadeia de eta_B com xi* livre em +-0.6%.
(ii) Gramática de fórmulas estruturais para be = X / (Y * n):
       X = a/b  (a<=26, b<=12)      Y = c/d (c<=32, d<=32)      n in 2..12
     (a fórmula do paper: (13/6) / ((32/25) * 6) = 325/1152)
     Conta quantas fórmulas reproduzem, COM A MESMA TOLERÂNCIA DA PH:
       - o invariante mestre rho_L G_N^2   (erro paper: 0.93%)
       - m_e do cascade                     (erro paper: 0.42%)
       - os dois juntos (teste conjunto)
     e quantas têm complexidade <= a da fórmula do paper
     (complexidade = a+b+c+d+n; paper: 13+6+32+25+6 = 82).

Só biblioteca padrão. Valores e fórmulas: PH v85 (eq. 6, eq. 9, §J.4, §12).
Limitação: é um teste de ordem de grandeza; a gramática é escolha minha, e o
resultado depende dela. Reportar sempre junto com a definição da gramática.
"""
import math
from fractions import Fraction

Nc0, Ne = 6, 10
K = 1.16335
GN = 6.70883e-39                                   # GeV^-2 (CODATA 2018)
OBS_INV = 2.515e-47 * GN ** 2
ME_OBS = 0.511e-3                                  # GeV
G2 = 2 * math.pi ** 2 / 3                          # g^2 (ponto fixo UV)
BE_PAPER = 325 / 1152
ETA_OBS, ETA_N, XI0 = 6.13e-10, 0.749, 434.09


def inv(be, nc=Nc0):
    S = nc * Ne / be
    return nc / (27 * math.pi ** 4) * K ** 2 * math.exp(-S * (1 + 3 / math.pi ** 2))


def me(be, nc=Nc0):
    S = nc * Ne / be
    return math.sqrt(K / GN) * math.exp(-S / (2 * G2) - Ne / be)


def err(p, o):
    return abs(p - o) / o


def window(fn, obs, tol, be0):
    """janela de be onde |fn(be)-obs|/obs <= tol (fn monotônica)."""
    lo, hi = be0 * 0.8, be0 * 1.2
    for _ in range(200):
        m = (lo + hi) / 2
        if fn(m) < obs:
            lo = m
        else:
            hi = m
    star = (lo + hi) / 2
    def edge(sign):
        a, b = star, star * (1 + 0.05 * sign)
        for _ in range(200):
            m = (a + b) / 2
            if err(fn(m), obs) <= tol:
                a = m
            else:
                b = m
        return a
    l, h = sorted([edge(-1), edge(1)])
    return star, l, h


def eta(eps, xi):
    return eps * ETA_N / xi ** 3


def eps_psi(j):
    d = 1 / (2 * j) ** 2
    return 4 * d / (1 + d) ** 2


def test_i():
    print("=== (i) eta_B com xi* livre em +-0.6% ===")
    p = eps_psi(2.5) * (4 / 9) * ETA_N / XI0 ** 3
    tol = err(p, ETA_OBS)
    print(f"previsão paper {p:.4e}, erro {tol*100:.2f}% (tolerância usada)")
    js = [0.5, 1.5, 2.5, 3.5, 4.5]
    fr = sorted({Fraction(a, b) for b in range(1, 14) for a in range(1, b + 1)})
    grid = [XI0 * (1 + k * 0.006 / 50) for k in range(-50, 51)]
    for label, xs in (("xi* fixo", [XI0]), ("xi* livre +-0.6%", grid)):
        hits = []
        for j in js:
            for q in fr:
                if any(err(eps_psi(j) * float(q) * ETA_N / x ** 3, ETA_OBS) <= tol
                       for x in xs):
                    hits.append((j, q))
        print(f"  {label}: {len(hits)} de {len(js)*len(fr)} combinações "
              f"({100*len(hits)/(len(js)*len(fr)):.1f}%)")
    print()


def test_ii():
    tol_inv = err(inv(BE_PAPER), OBS_INV)
    tol_me = err(me(BE_PAPER), ME_OBS)
    print("=== (ii) gramática be = X / (Y * n) ===")
    print(f"paper: inv {inv(BE_PAPER):.4e} (erro {tol_inv*100:.3f}%), "
          f"m_e {me(BE_PAPER)*1e3:.4f} MeV (erro {tol_me*100:.3f}%)")
    s_i, li, hi = window(inv, OBS_INV, tol_inv, BE_PAPER)
    s_m, lm, hm = window(me, ME_OBS, tol_me, BE_PAPER)
    print(f"janela inv: [{li:.6f}, {hi:.6f}]  be*={s_i:.6f}")
    print(f"janela m_e: [{lm:.6f}, {hm:.6f}]  be*={s_m:.6f}")
    ov_lo, ov_hi = max(li, lm), min(hi, hm)
    print(f"interseção: [{ov_lo:.6f}, {ov_hi:.6f}]  "
          f"({'vazia' if ov_lo > ov_hi else 'não vazia'})\n")

    Xs = sorted({Fraction(a, b) for b in range(1, 13) for a in range(1, 27)})
    Ys = sorted({Fraction(c, d) for d in range(1, 33) for c in range(1, 33)})
    ns = range(2, 13)
    total = len(Xs) * len(Ys) * len(ns)
    hit_inv = hit_me = hit_both = 0
    cplx_paper = 13 + 6 + 32 + 25 + 6
    simple_inv = simple_both = 0
    best = []
    for X in Xs:
        xa, xb = X.numerator, X.denominator
        for Y in Ys:
            yc, yd = Y.numerator, Y.denominator
            base = float(X / Y)
            for n in ns:
                be = base / n
                a = li <= be <= hi
                b = lm <= be <= hm
                if a:
                    hit_inv += 1
                    c = xa + xb + yc + yd + n
                    if c <= cplx_paper:
                        simple_inv += 1
                if b:
                    hit_me += 1
                if a and b:
                    hit_both += 1
                    c = xa + xb + yc + yd + n
                    if c <= cplx_paper:
                        simple_both += 1
                    best.append((c, X, Y, n))
    best.sort(key=lambda t: t[0])
    print(f"fórmulas na gramática: {total:,}")
    print(f"  reproduzem o invariante (0.93%): {hit_inv}")
    print(f"  reproduzem m_e (0.42%):          {hit_me}")
    print(f"  reproduzem os DOIS:               {hit_both}")
    print(f"  com complexidade <= {cplx_paper} (a do paper): "
          f"invariante {simple_inv}, conjunto {simple_both}")
    print("  5 fórmulas conjuntas mais simples (complexidade, X, Y, n):")
    for c, X, Y, n in best[:5]:
        print(f"     {c}: be = ({X}) / (({Y}) * {n}) = {float(X/Y/n):.6f}")
    print()
    paper_in = li <= BE_PAPER <= hi and lm <= BE_PAPER <= hm
    print(f"be do paper dentro das duas janelas: {paper_in}")
    return hit_both


if __name__ == "__main__":
    test_i()
    test_ii()
