#!/usr/bin/env python3
"""P14 / WP-UNW (ver P14UNW_DEFINITION.md)."""
import math
be = 325 / 1152; S = 6 * 10 / be; LAM = 1.2609e12
opts = {"A (1/b_e por degrau)": lambda dl: math.exp(-dl / be),
        "B (6/b_e por degrau)": lambda dl: math.exp(-dl * 6 / be),
        "C (barreira S por evento)": lambda dl: math.exp(-S)}
req = [("n -> nu_tau gamma (mesma tupla)", 7, 5.5e-30), ("p -> nu pi+ (contato)", 7, 7.8e-9),
       ("n -> 3 nu (contato)", 7, 3.4e-7), ("V1: p -> (0,4,5)_L gamma", 2, 5.5e-30)]
print(f"b_e = {be:.6f}, S = {S:.3f}")
print("##### Parte 1 #####")
for name, A in opts.items():
    line = f"  {name:28s}:"
    for lab, dl, lim in req:
        a = A(dl)
        line += f"  [{lab}: A = {a:.1e} {'passa' if a <= lim else 'FALHA'}]"
    print(line)
print("  V1 independentemente excluida: preve leptons de materia Q = +1, +2 a m(10) = 0,51 MeV, nunca vistos")
print("##### Parte 2 #####")
for lab, T in (("WP-SIM/NOW T_RH", 1.38e11), ("WP-TRH pico", 1.4e10), ("WP-TRH borda alta", 2e10),
               ("WP-RUNG lado frio (max)", 6.5e8), ("WP-RUNG lado quente (max)", 4.8e10)):
    print(f"  {lab:28s}: T/Lambda_meta = {T/LAM:.1e}  -> {'abaixo do confinamento' if T < LAM else 'acima'}")
