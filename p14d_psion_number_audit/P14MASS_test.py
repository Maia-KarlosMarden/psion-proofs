#!/usr/bin/env python3
"""P14 / WP-M: massa de tupla e conjugada (C_A) no mesmo nivel de cascata."""
import math

LAM, BE, A, P = 1.2609e12, 325 / 1152, 13 / 12, 5 / 6
mk = [k - 2.5 for k in range(6)]
conjA = lambda t: tuple(sorted((6 - k) % 6 for k in t))


def m_level(l):
    return LAM * math.exp(-l / BE)


def mJR_content(t):
    return abs(sum(mk[k] for k in t if k >= 3))


def mass(t, l, reading):
    base = m_level(l)
    if reading == "M1":
        return base
    if reading == "M2":
        return base * math.exp(A * mJR_content(t) ** P / BE)
    if reading == "M3":
        return base * math.exp(A * 0.0 ** P / BE)   # estado fundamental: rotulo de excitacao 0
    raise ValueError


print("##### T1/T2: massas no mesmo nivel #####")
pairs = [("eletron/positron", (0, 1, 2), 10, 0.000511), ("proton/antiproton", (0, 4, 5), 8, 0.938272)]
for reading in ("M1", "M2", "M3"):
    print(f"  leitura {reading}:")
    for name, t, l, mobs in pairs:
        ct = conjA(t)
        m1, m2 = mass(t, l, reading), mass(ct, l, reading)
        print(f"    {name:18s} l={l}: m{t} = {m1:.4e} GeV ; m{ct} = {m2:.4e} GeV ; "
              f"razao = {m1/m2:.4e} ; observado {mobs} GeV -> m(t)/obs = {m1/mobs:.3g}")
print(f"\n  |m_J,R| por conteudo: (0,1,2) -> {mJR_content((0,1,2))}, (0,4,5) -> {mJR_content((0,4,5))}")
print(f"  fator M2 para |m_J,R| = 4: exp(A 4^p / b_e) = {math.exp(A * 4 ** P / BE):.4e}")

print("\n##### T3: Hamiltoniano interno de J.6.3 #####")
habs = [abs(x) for x in mk]
for t in ((0, 1, 2), (0, 4, 5)):
    print(f"  <{t}| diag|m| |{t}> = {sum(habs[k] for k in t)}  (o salto 2<->3 nao contribui na diagonal)")
print(f"  diferenca (0,1,2) - (0,4,5) = {sum(habs[k] for k in (0,1,2)) - sum(habs[k] for k in (0,4,5))} unidades de h")
print(f"  espacamento entre niveis: m(l)/m(l+1) = exp(1/b_e) = {math.exp(1/BE):.2f}")

ok = {r: all(abs(mass(t, l, r) / mass(conjA(t), l, r) - 1) < 1e-12 for _, t, l, _ in pairs) for r in ("M1", "M2", "M3")}
print(f"\n  VEREDITO (degenerescencia exata no mesmo nivel): {ok}")
