import math
import numpy as np
Q9 = np.array([1, -5, -5, 10, 4, 4]) / 9
def pos(t, G):
    if G in ("G1", "G3"):
        return np.array([[math.cos(math.pi * k / 3), math.sin(math.pi * k / 3)] for k in t])
    if G == "G2":
        return np.array([[k - 2.5, 0.0] for k in t])
    if G == "G4":
        return np.array([[math.cos(2 * math.pi * i / 3), math.sin(2 * math.pi * i / 3)] for i in range(3)]) / math.sqrt(3)
def r2(t, G):
    P = pos(t, G); q = Q9[list(t)]
    c = np.zeros(2) if G == "G3" else P.mean(axis=0)
    s = float(np.sum(q * np.sum((P - c) ** 2, axis=1)))
    Qt = q.sum()
    return s / Qt if abs(Qt) > 1e-12 else s
meas = dict(Rn=(-0.1161 / 0.7071, 0.003), RS=(0.61 / 0.7071, 0.15 / 0.7071))
for G in ("G1", "G2", "G3", "G4"):
    rp, rn, rS = r2((0, 4, 5), G), r2((0, 1, 5), G), r2((0, 1, 2), G)
    Rn, RS = rn / rp, rS / rp
    okn = abs(Rn / meas['Rn'][0] - 1) < 0.10; okS = abs(RS - meas['RS'][0]) < meas['RS'][1]
    a = math.sqrt(0.7071 / rp) if rp > 0 else float('nan')
    print(f"{G}: <r2>_p = {rp:+.4f} a^2, <r2>_n = {rn:+.4f} a^2, <r2>_(0,1,2) = {rS:+.4f} a^2 | R_n = {Rn:+.3f} (med -0,164) "
          f"{'ACERTA' if okn else 'erra'} | R_S = {RS:+.3f} (med 0,86+-0,21) {'ACERTA' if okS else 'erra'} | a_p = {a:.3f} fm")
    if rS > 0:
        ae = 2e-20 / math.sqrt(rS) * 1e15  # fm
        print(f"     eletron: para raio < 2e-20 m exige a_e < {ae:.1e} fm -> a_e/a_p < {ae/a:.1e}; Compton daria a_e/a_p = {610.135/0.509:.0f}")
