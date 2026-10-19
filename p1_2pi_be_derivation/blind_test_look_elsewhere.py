#!/usr/bin/env python3
"""
Teste cego de look-elsewhere para a Hipótese do Psion (PH v85).

Pergunta: dado o espaço de escolhas estruturais que alguém poderia ter feito,
quantas reproduzem os observáveis com a MESMA tolerância que a PH reporta?

Teste A: invariante mestre rho_L * G_N^2, varrendo be racional a/b.
Teste B: cadeia de eta_B, varrendo j (spin) e eps_conf racional a/b.
Teste C: sensibilidade (quanta precisão cada entrada precisa ter).

Só biblioteca padrão. Uso: python3 blind_test_look_elsewhere.py
Fórmulas e valores: PH v85 (eq. 9, §5.1, §J.4). Reconferir contra o repositório.
"""
import math
from fractions import Fraction

# ---------- constantes do paper ----------
Nc, Ne = 6, 10
K = 1.16335
OBS_INV = 2.515e-47 * (6.70883e-39) ** 2      # rho_L,obs * G_N^2  (~1.132e-123)
BE_PAPER = 325 / 1152
ETA_OBS = 6.13e-10                             # cadeia base Planck 2018
ETA_N = 0.749
XI_STAR = 434.09


def master_invariant(be, ne=Ne, nc=Nc):
    """rho_L G_N^2 = C K^2 exp(-S (1+3/pi^2)),  S = nc*ne/be,  C = nc/(27 pi^4)."""
    S = nc * ne / be
    C = nc / (27 * math.pi ** 4)
    return C * K ** 2 * math.exp(-S * (1 + 3 / math.pi ** 2))


def rel_err(pred, obs):
    return abs(pred - obs) / obs


# =====================================================================
# Teste A: be racional
# =====================================================================
def test_A(bmax=2000):
    pred_paper = master_invariant(BE_PAPER)
    tol = rel_err(pred_paper, OBS_INV)           # tolerância = erro da própria PH
    print("=== TESTE A: invariante mestre vs. be racional ===")
    print(f"previsão com be=325/1152: {pred_paper:.4e}   observado: {OBS_INV:.4e}")
    print(f"erro da PH: {tol*100:.3f}%  (usado como tolerância)\n")

    # janela em be: acha por bisseção onde |erro| <= tol (função monotônica em be)
    def f(be):
        return rel_err(master_invariant(be), OBS_INV)

    # be_alvo: onde previsão = observado (busca em torno de be_paper)
    lo, hi = BE_PAPER * 0.9, BE_PAPER * 1.1
    for _ in range(200):
        mid = (lo + hi) / 2
        if master_invariant(mid) < OBS_INV:
            lo = mid
        else:
            hi = mid
    be_star = (lo + hi) / 2
    # bordas da janela
    def edge(sign):
        a, b = be_star, be_star * (1 + sign * 0.05)
        for _ in range(200):
            m = (a + b) / 2
            if f(m) <= tol:
                a = m
            else:
                b = m
        return a
    be_lo, be_hi = sorted([edge(-1), edge(+1)])
    print(f"be exato que casa: {be_star:.6f}  (paper: {BE_PAPER:.6f})")
    print(f"janela de be dentro da tolerância: [{be_lo:.6f}, {be_hi:.6f}]  "
          f"(largura {be_hi-be_lo:.2e}, {100*(be_hi-be_lo)/be_star:.3f}% de be)\n")

    hits = set()
    for b in range(1, bmax + 1):
        a_lo = math.ceil(be_lo * b)
        a_hi = math.floor(be_hi * b)
        for a in range(a_lo, a_hi + 1):
            hits.add(Fraction(a, b))
    hits = sorted(hits, key=lambda q: (q.denominator, q.numerator))
    print(f"racionais a/b (b <= {bmax}) dentro da janela: {len(hits)}")
    for B in (50, 100, 300, 500, 1152, bmax):
        n = sum(1 for q in hits if q.denominator <= B)
        print(f"   denominador <= {B:5d}: {n}")
    print("   5 mais simples (menor denominador):",
          ", ".join(f"{q.numerator}/{q.denominator}" for q in hits[:5]))
    paper = Fraction(325, 1152)
    print(f"   325/1152 está entre os hits? {paper in hits}\n")
    return len(hits)


# =====================================================================
# Teste B: eta_B
# =====================================================================
def eps_psi(j):
    delta = 1 / (2 * j) ** 2
    return 4 * delta / (1 + delta) ** 2


def eta_chain(j, eps_conf, n_ngamma=ETA_N, xi=XI_STAR):
    return eps_psi(j) * eps_conf * n_ngamma / xi ** 3


def test_B(max_den=13):
    paper = eta_chain(2.5, 4 / 9)
    tol = rel_err(paper, ETA_OBS)
    print("=== TESTE B: cadeia de eta_B ===")
    print(f"previsão (j=5/2, eps_conf=4/9, xi*={XI_STAR}): {paper:.4e}  "
          f"vs {ETA_OBS:.2e}  erro {tol*100:.2f}%")
    js = [0.5, 1.5, 2.5, 3.5, 4.5]
    fr = sorted({Fraction(a, b) for b in range(1, max_den + 1)
                 for a in range(1, b + 1)})
    total = len(js) * len(fr)
    hits = []
    for j in js:
        for q in fr:
            if rel_err(eta_chain(j, float(q)), ETA_OBS) <= tol:
                hits.append((j, q))
    print(f"combinações testadas (j x eps_conf=a/b, b<={max_den}): {total}")
    print(f"combinações dentro da tolerância: {len(hits)}")
    for j, q in hits:
        print(f"   j={j}, eps_conf={q}  -> eta_B={eta_chain(j, float(q)):.3e}")
    print("   (xi* e n/n_gamma mantidos nos valores do paper; ver Teste C)\n")


# =====================================================================
# Teste C: sensibilidade
# =====================================================================
def test_C():
    print("=== TESTE C: sensibilidade ===")
    # eta_B ~ xi^-3
    print("eta_B ~ xi*^-3:  1% em xi* -> ~3% em eta_B. "
          "Acerto de 1.7% exige xi* conhecido a ~0.6%.")
    # invariante vs be
    h = 1e-6
    d = (math.log(master_invariant(BE_PAPER * (1 + h))) -
         math.log(master_invariant(BE_PAPER))) / h
    print(f"d ln(rho_L G_N^2) / d ln(be) = {d:.1f}  "
          f"(1% em be muda o invariante em ~{abs(d):.0f}%)")
    print("Logo be precisa ser derivado a precisão de ~0.01% para sustentar o 0.93%;")
    print("o status atual de be (P1, T2, 28% de lacuna teórica) não garante isso.\n")


if __name__ == "__main__":
    test_A()
    test_B()
    test_C()
