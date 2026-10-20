#!/usr/bin/env python3
"""P14 / WP-NOW: a PH como esta (ontologia O-), etapas A (contato) e B (congelamento p pbar)."""
import io, contextlib, json, math, runpy
import numpy as np
from scipy.integrate import solve_ivp

with contextlib.redirect_stdout(io.StringIO()):
    S = runpy.run_path('P14SIM_network.py', run_name='lib')
TR, NT, REARR = S['TR'], S['NT'], S['REARR']
V3, H0, OJS = S['V3'], S['H0'], S['OJS']
H0t = (V3.conj() @ H0 @ V3.T).real
Ot = (V3.conj() @ OJS @ V3.T).real
gamma_over_H, T_RH, T_END, Y0 = S['gamma_over_H'], S['T_RH'], S['T_END'], S['Y0']
K_EQ = 10.0
RTOL = 1e-11
N = 6
idx = np.array(TR)
IP, IE = TR.index((0, 4, 5)), TR.index((0, 1, 2))


def weights(g, law, phi, sign):
    A = H0t + g * np.exp(1j * sign * phi) * Ot
    diag = (H0t + g * Ot).diagonal()          # diagonal real (O real, fase so fora da diagonal)
    W = np.abs(A) ** 2
    np.fill_diagonal(W, 0.0)
    kb = diag ** 2 if law == "L1" else np.exp(-(diag - diag.min()))
    return kb / kb.mean(), W


def stageA(cG, law, delta, phi, g=1.0):
    kb, Wm = weights(g, law, phi, +1)
    _, Wa = weights(g, law, phi, -1)
    Yeq = 1 / 12
    Yeq_t = K_EQ * Yeq ** 3
    a, b, c, d = REARR.T

    def f(x, Y):
        r = gamma_over_H(T_RH * math.exp(-x), cG)
        psi, apsi = Y[:6], Y[6:12]
        t, at = Y[12:12 + NT], Y[12 + NT:]
        dpsi, dapsi = np.zeros(6), np.zeros(6)
        dt, dat = np.zeros(NT), np.zeros(NT)
        for P, dP, Tt, dT, W in ((psi, dpsi, t, dt, Wm), (apsi, dapsi, at, dat, Wa)):
            J = kb * (P[idx[:, 0]] * P[idx[:, 1]] * P[idx[:, 2]] - Tt / K_EQ)
            dT += J
            for col in range(3):
                np.subtract.at(dP, idx[:, col], J)
            dT += W @ Tt - W.sum(axis=1) * Tt
            Jr = Tt[a] * Tt[b] - Tt[c] * Tt[d]
            np.subtract.at(dT, a, Jr); np.subtract.at(dT, b, Jr)
            np.add.at(dT, c, Jr); np.add.at(dT, d, Jr)
        Ja = psi * apsi - Yeq ** 2               # aniquilacao de pares de psions
        dpsi -= Ja; dapsi -= Ja
        Jt = t * at - Yeq_t ** 2                 # aniquilacao tupla-antitupla
        dt -= Jt; dat -= Jt
        return r * np.concatenate([dpsi, dapsi, dt, dat])

    Y = np.concatenate([np.full(6, (1 + delta) / 12), np.full(6, (1 - delta) / 12), np.zeros(2 * NT)])
    sol = solve_ivp(f, [0, math.log(T_RH / T_END)], Y, method="LSODA", rtol=RTOL, atol=1e-20)
    Yf = sol.y[:, -1]
    t, at = Yf[12:12 + NT], Yf[12 + NT:]
    Nnet0 = 6 * (2 * delta / 12)
    Nnet = (Yf[:6] - Yf[6:12]).sum() + 3 * (t - at).sum()
    return dict(dp=t[IP] - at[IP], de=t[IE] - at[IE], Nnet=Nnet, Nerr=abs(Nnet - Nnet0), ok=sol.success)


def stageB(YB, m=0.938, sv=1 / 0.135 ** 2, gstar=17.25, g=2):
    """Correcao registrada: integra Y_pbar em variavel linear com Radau e jacobiano (a versao
    em log estourava). x de 1 a 1000."""
    MPL = 1.2209e19
    lam = 0.264 * (gstar / math.sqrt(gstar)) * MPL * m * sv
    Yeq = lambda x: 0.145 * (g / gstar) * x ** 1.5 * math.exp(-x)
    f = lambda x, y: [-(lam / x ** 2) * (y[0] * (y[0] + YB) - Yeq(x) ** 2)]
    jac = lambda x, y: [[-(lam / x ** 2) * (2 * y[0] + YB)]]
    y0 = (-YB + math.sqrt(YB ** 2 + 4 * Yeq(1.0) ** 2)) / 2
    sol = solve_ivp(f, [1.0, 1000.0], [y0], method="Radau", jac=jac, rtol=1e-8, atol=1e-35)
    ypbar = max(sol.y[0, -1], 0.0)
    return ypbar + YB, ypbar


if __name__ == "__main__":
    out = {}
    print("##### Etapa A: rede de alta temperatura com materia e antimateria #####")
    print(f"  {'c_G':>6s} {'lei':>3s} {'delta':>6s} {'phi':>5s} {'Delta_p':>11s} {'Delta_e':>11s} {'N_psi liq':>11s} {'erro N':>8s}")
    for cG in (1e-5, 1e-3, 1.0):
        for law in ("L1", "L2"):
            for delta in (0.0, 1e-8):
                for phi in (0.0, math.pi / 4):
                    r = stageA(cG, law, delta, phi)
                    out[f"{cG}_{law}_{delta}_{phi:.3f}"] = {k: float(v) for k, v in r.items() if k != "ok"}
                    print(f"  {cG:6.0e} {law:>3s} {delta:6.0e} {phi:5.2f} {r['dp']:+11.3e} {r['de']:+11.3e} "
                          f"{r['Nnet']:+11.3e} {r['Nerr']:8.1e}")
    c1 = max(abs(v['dp']) + abs(v['de']) + abs(v['Nnet']) for k, v in out.items() if k.split('_')[2] == '0.0')
    c2 = max(v['Nerr'] for v in out.values())
    print(f"\n  C1/C3 (delta = 0, com e sem phi): maior |assimetria| = {c1:.1e}  {'OK' if c1 < 1e-15 else 'FALHA'}")
    print(f"  C2 conservacao de N_psi: maior erro = {c2:.1e}  {'OK' if c2 < 1e-15 else 'FALHA'}")
    phi_eff = max(abs(out[k]['dp'] - out[k.rsplit('_', 1)[0] + '_0.000']['dp'])
                  for k in out if not k.endswith('_0.000'))
    print(f"  efeito de phi sobre Delta_p (maior diferenca): {phi_eff:.1e}")

    print("\n##### Etapa B: congelamento da aniquilacao p pbar #####")
    for m in (0.938, 0.610):
        for label, YB in (("sem assimetria (PH como esta)", 0.0),
                          ("assimetria posta a mao, Y_B = 8.7e-11", 8.7e-11)):
            Yp, Ypb = stageB(YB, m=m)
            print(f"  m = {m} GeV, {label:40s}: n_p/n_gamma = {7.04*Yp:.2e}, n_pbar/n_gamma = {7.04*Ypb:.2e}")
    print("\n  comparacao (so agora): observado n_B/n_gamma = 6.1e-10, sem antiprotons primordiais")
    json.dump(out, open("P14NOW_results.json", "w"), indent=1)
