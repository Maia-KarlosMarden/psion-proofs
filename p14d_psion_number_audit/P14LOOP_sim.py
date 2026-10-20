#!/usr/bin/env python3
"""P14 / WP-LOOP: fase em J.Sigma com dinamica quantica da montagem (ver P14LOOP_DEFINITION.md)."""
import io, contextlib, json, math, runpy, sys, time
import numpy as np
from scipy.integrate import solve_ivp

with contextlib.redirect_stdout(io.StringIO()):
    S = runpy.run_path('P14NOW_sim.py', run_name='lib')
    N_ = runpy.run_path('P14SIM_network.py', run_name='lib')
TR, NT, REARR, weights = S['TR'], S['NT'], S['REARR'], S['weights']
gamma_over_H, T_RH, T_END, K_EQ = S['gamma_over_H'], S['T_RH'], S['T_END'], S['K_EQ']
V3, H0, O = N_['V3'], N_['H0'], N_['O']
TI = {t: i for i, t in enumerate(TR)}
IP, IE, IN1, IN2 = TI[(0, 4, 5)], TI[(0, 1, 2)], TI[(1, 2, 3)], TI[(0, 2, 4)]
EXO = [t for t in range(NT) if t not in (IP, IE, IN1, IN2)]
idx = np.array(TR)
ra, rb, rc, rd = REARR.T
YEQ, YEQT = 1 / 12, K_EQ * (1 / 12) ** 3
OH = (O + O.conj().T) / 2
OA = (O - O.conj().T) / 2j
proj = lambda M: V3.conj() @ M @ V3.T
Ht0, OHt, OAt = proj(H0), proj(OH), proj(OA)
Y0 = 0.208 * 12 / 106.75
ETA = 6.1e-10
DT = ETA / (7.04 * Y0)
B = np.r_[1.0, 1.0, -2.0, np.zeros(len(EXO) + 5)]


def Hphi(phi, g=1.0):
    return Ht0 + g * (math.cos(phi) * OHt + math.sin(phi) * OAt)


def make_rhs(law, phi, omega, gam, rate):
    kb, _ = weights(1.0, law, 0.0, 1)
    Hm = Hphi(phi)
    Ha = Hm.conj()
    sink0 = kb / K_EQ

    def rear(p):
        Jr = p[ra] * p[rb] - p[rc] * p[rd]
        d = np.zeros(NT)
        np.subtract.at(d, ra, Jr); np.subtract.at(d, rb, Jr)
        np.add.at(d, rc, Jr); np.add.at(d, rd, Jr)
        return d

    def f(x, y):
        r = rate(x)
        psi, apsi = y[:6], y[6:12]
        rm = (y[12:412] + 1j * y[412:812]).reshape(NT, NT)
        ram = (y[812:1212] + 1j * y[1212:1612]).reshape(NT, NT)
        pm, pa = rm.diagonal().real, ram.diagonal().real
        outs = []
        dps = []
        for P, rho, H, pbar in ((psi, rm, Hm, pa), (apsi, ram, Ha, pm)):
            P3 = P[idx[:, 0]] * P[idx[:, 1]] * P[idx[:, 2]]
            p = rho.diagonal().real
            J = kb * (P3 - p / K_EQ)
            dP = np.zeros(6)
            for col in range(3):
                np.subtract.at(dP, idx[:, col], J)
            s = sink0 + pbar
            d = -1j * omega * (H @ rho - rho @ H) - gam * (rho - np.diag(rho.diagonal()))
            d -= 0.5 * (s[:, None] * rho + rho * s[None, :])
            d += np.diag(kb * P3 + YEQT ** 2 + rear(p))
            outs.append(d); dps.append(dP)
        Ja = psi * apsi - YEQ ** 2
        dpsi, dapsi = dps[0] - Ja, dps[1] - Ja
        dm, da = outs
        return r * np.concatenate([dpsi, dapsi, dm.real.ravel(), dm.imag.ravel(), da.real.ravel(), da.imag.ravel()])
    return f


def initial(eq=False):
    y = np.zeros(1612)
    y[:12] = 1 / 12
    if eq:
        for off in (12, 812):
            y[off:off + 400] = (YEQT * np.eye(NT)).ravel()
    return y


def summarize(y):
    rm = y[12:412].reshape(NT, NT).diagonal()
    ra_ = y[812:1212].reshape(NT, NT).diagonal()
    dt = rm - ra_
    dfree = y[:6] - y[6:12]
    d = np.r_[dt[IP], dt[IE], dt[IN1] + dt[IN2], dt[EXO], dfree[:5]]
    return dict(dt=dt, dfree=dfree, d=d, Nnet=dfree.sum() + 3 * dt.sum(),
                trace_m=rm.sum(), trace_a=ra_.sum())


def run(law, phi, omega, gam, cG=1e-3):
    f = make_rhs(law, phi, omega, gam, lambda x: gamma_over_H(T_RH * math.exp(-x), cG))
    sol = solve_ivp(f, [0, math.log(T_RH / T_END)], initial(), method="DOP853", rtol=1e-10, atol=1e-22)
    return summarize(sol.y[:, -1]), sol.success


def score(d):
    s = (d @ B) / (d @ d) if d @ d > 0 else 0.0
    return np.linalg.norm(s * d - B) / np.linalg.norm(B), s


if __name__ == "__main__":
    t0 = time.time()
    J = max(abs((Hphi(0.3)[a, b] * Hphi(0.3)[b, c] * Hphi(0.3)[c, a]).imag)
            for a in range(NT) for b in range(NT) for c in range(NT))
    print(f"invariante de CP max|Im H_ab H_bc H_ca| (phi = 0,3) = {J:.3f}; entradas de O_A nas tuplas = "
          f"{int((abs(OAt) > 1e-12).sum())}")
    print("\n##### Controles #####")
    c1 = max(abs(run(l, 0.0, 1.0, 1.0)[0]['d']).max() for l in ("L2", "L1"))
    print(f"  C1 phi = 0: maior |assimetria| = {c1:.1e}  {'OK' if c1 < 1e-12 else 'FALHA'}")
    c4 = 0.0
    for l in ("L2", "L1"):
        f = make_rhs(l, 0.3, 1.0, 1.0, lambda x: 1e3)
        sol = solve_ivp(f, [0, 20], initial(eq=True), method="DOP853", rtol=1e-10, atol=1e-22)
        c4 = max(c4, abs(summarize(sol.y[:, -1])['d']).max())
    print(f"  C4 equilibrio com phi = 0,3: maior |assimetria| = {c4:.1e}  {'OK' if c4 < 1e-12 else 'FALHA'}")
    print(f"  ({time.time()-t0:.0f} s)")

    print("\n##### Grade: c_G = 1e-3, phi = 1e-3 #####")
    phi = 1e-3
    res, c2, c3 = {}, 0.0, 0.0
    best = {}
    for law in ("L2", "L1"):
        for om in (0.1, 1.0, 10.0):
            for gm in (0.1, 1.0, 10.0):
                r, ok = run(law, phi, om, gm)
                rho, s = score(r['d'])
                c2 = max(c2, abs(r['Nnet']))
                need = abs(DT / r['dt'][IP]) * phi if r['dt'][IP] != 0 else float('inf')
                res[f"{law}_{om}_{gm}"] = dict(rho=float(rho), dp=float(r['dt'][IP]), de=float(r['dt'][IE]),
                                               dnu=float(r['dt'][IN1] + r['dt'][IN2]),
                                               exo_max=float(abs(r['dt'][EXO]).max()),
                                               free_max=float(abs(r['dfree']).max()), phi_needed=float(need), ok=ok)
                print(f"  {law} omega={om:5.1f} gamma={gm:5.1f}: rho = {rho:.4f}  Dp/phi = {r['dt'][IP]/phi:+.2e}  "
                      f"De/Dp = {r['dt'][IE]/r['dt'][IP]:+.2f}  Dnu/Dp = {(r['dt'][IN1]+r['dt'][IN2])/r['dt'][IP]:+.2f}  "
                      f"max|exo|/|Dp| = {abs(r['dt'][EXO]).max()/abs(r['dt'][IP]):.2e}  phi p/ eta = {need:.2e}  ok={ok}")
                if law not in best or rho < best[law][0]:
                    best[law] = (rho, om, gm)
        print(f"  ({time.time()-t0:.0f} s)")
    print(f"\n  C2 N_psi: maior |N_psi| = {c2:.1e}  {'OK' if c2 < 1e-12 else 'FALHA'}")

    print("\n##### Melhor de cada lei: C3, linearidade, sensibilidade a c_G #####")
    for law, (rho, om, gm) in best.items():
        rp = run(law, phi, om, gm)[0]; rn = run(law, -phi, om, gm)[0]; r10 = run(law, 10 * phi, om, gm)[0]
        big = abs(rp['d']) > 1e-14
        c3 = max(c3, abs(rn['d'][big] / rp['d'][big] + 1).max())
        lin = abs(r10['dt'][IP] / (10 * rp['dt'][IP]) - 1)
        line = f"  {law} omega={om} gamma={gm}: rho = {rho:.4f}; linearidade Dp (phi 1e-2 vs 1e-3) {lin:.1e}"
        for cG in (1e-5, 1.0):
            rc_ = run(law, phi, om, gm, cG=cG)[0]
            line += f"; c_G={cG:.0e}: rho = {score(rc_['d'])[0]:.4f}, Dp/phi = {rc_['dt'][IP]/phi:+.1e}"
        print(line)
    print(f"  C3 Delta(-phi) = -Delta(phi): maior desvio {c3:.1e}  {'OK' if c3 < 1e-6 else 'FALHA'}")
    bl = {"L2": 0.926, "L1": 0.787}
    for law, (rho, om, gm) in best.items():
        print(f"\n  CRITERIO {law}: melhor rho = {rho:.4f}; sucesso {'SIM' if rho < 0.1 else 'nao'}; "
              f"estrutura acima do acaso (1o percentil {bl[law]}) {'SIM' if rho < bl[law] else 'nao'}")
    print(f"  tempo total {time.time()-t0:.0f} s")
    json.dump(res, open("P14LOOP_results.json", "w"), indent=1)
