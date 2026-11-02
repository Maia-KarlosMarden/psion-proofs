#!/usr/bin/env python3
"""WP-CAL, adendo: A5/V5 refeito (correcao registrada do patch do setor do topo)."""
import io, contextlib, math, runpy
with contextlib.redirect_stdout(io.StringIO()):
    T_ = runpy.run_path('P14TPROG_sim.py', run_name='lib')
sim = T_['simulate']; G = sim.__globals__; base = G['sectors']; m_lev = G['m_lev']

def omega(l, kap):
    m = m_lev(l)
    G['sectors'] = lambda m8=0.938: [dict(name="topo", m=m, g=2, sv=kap / m ** 2)] + base(m8)[1:]
    r = sim(max(1e3 * m, 1e3), [("cool", 1e-5)], {}, dl=0.004)
    G['sectors'] = base
    return r['topo']['Oh2'], m

for kap in (1.0, 4 * math.pi):
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if omega(mid, kap)[0] > 0.12 else (lo, mid)
    oh, m = omega(hi, kap)
    print(f"kappa = {kap:5.2f}: l_DM = {hi:.3f}, m_DM = {m:.3e} GeV, Omega h2 = {oh:.3f}; "
          f"nivel 0,8 do v85 daria Omega h2 = {omega(0.8, kap)[0]:.1e}")
