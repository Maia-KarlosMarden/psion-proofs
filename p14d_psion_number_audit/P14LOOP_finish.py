#!/usr/bin/env python3
"""WP-LOOP, final: C3, linearidade e c_G = 1e-5 nos melhores pontos (c_G = 1 inviavel com integrador explicito)."""
import io, contextlib, runpy
with contextlib.redirect_stdout(io.StringIO()):
    L = runpy.run_path('P14LOOP_sim.py', run_name='lib')
run, score, IP = L['run'], L['score'], L['IP']
phi = 1e-3
c3 = 0.0
for law, om, gm in (("L2", 10.0, 1.0), ("L1", 10.0, 10.0)):
    rp = run(law, phi, om, gm)[0]; rn = run(law, -phi, om, gm)[0]; r10 = run(law, 10 * phi, om, gm)[0]
    big = abs(rp['d']) > 1e-14
    c3 = max(c3, abs(rn['d'][big] / rp['d'][big] + 1).max())
    lin = abs(r10['dt'][IP] / (10 * rp['dt'][IP]) - 1)
    rc = run(law, phi, om, gm, cG=1e-5)[0]
    print(f"{law} omega={om} gamma={gm}: rho = {score(rp['d'])[0]:.4f}; linearidade {lin:.1e}; "
          f"c_G=1e-5: rho = {score(rc['d'])[0]:.4f}, Dp/phi = {rc['dt'][IP]/phi:+.1e}")
print(f"C3 Delta(-phi) = -Delta(phi): maior desvio {c3:.1e}  {'OK' if c3 < 1e-6 else 'FALHA'}")
print("CRITERIO L2: melhor rho = 0.9821; sucesso nao; estrutura acima do acaso (1o percentil 0.926) nao")
print("CRITERIO L1: melhor rho = 0.9616; sucesso nao; estrutura acima do acaso (1o percentil 0.787) nao")
