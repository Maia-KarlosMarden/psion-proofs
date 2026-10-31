# ADENDO POS-EXECUCAO (registrado depois de ver L<1): massas, gap, estabilidade dos otimos de WP-EXH. Nao cego.
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import runpy, io, contextlib
with contextlib.redirect_stdout(io.StringIO()): X=runpy.run_path('P14EXH_calc.py',run_name='lib')
TR,AR,T1,S1,gsr,Rs,loss=X['TR'],X['AR'],X['T1'],X['S1'],X['gsr'],X['Rs'],X['loss']
def en(h,t,s,k,l):
    A,B,C,R2=AR[h]; dg=A+k*B+l*C; w=eigsh((sp.diags(dg)-t*T1-s*S1).tocsc(),k=3,which='SA')[0]; w.sort(); return w
pts=[(0.001,0.069,0.26,0.11),(0.000,0.065,0.19,0.0)]
for x in pts:
    t,s,k,l=x; E={h:en(h,*x) for h in TR}; e0={h:E[h][0] for h in TR}
    dn=e0['n']-e0['p']; dL=e0['L']-e0['p']
    print(f"ponto {x}: L={loss(np.array(x)):.3f}")
    print("  gaps (E1-E0) por hadron: "+", ".join(f"{h}:{E[h][1]-E[h][0]:.2e}" for h in TR))
    print(f"  E_n-E_p={dn:+.4f}  E_L-E_p={dL:+.4f}  E_S0-E_L={e0['S0']-e0['L']:+.4f}  E_Sp1-E_p={e0['Sp1']-e0['p']:+.4f}  E_Sp2-E_p={e0['Sp2']-e0['p']:+.4f}")
    if dL>0:
        u=177.4/dL; print(f"  escala u={u:.1f} MeV/un: Sp2-p={u*(e0['Sp2']-e0['p']):.0f} (251), Sp1-p={u*(e0['Sp1']-e0['p']):.0f}, S0-L={u*(e0['S0']-e0['L']):.0f} (77), n-p={u*dn:.2f} (1.29)")
    rng=np.random.default_rng(22); ok=0;N=40
    for _ in range(N):
        y=np.array(x)*(1+rng.uniform(-0.1,0.1,4))
        if loss(y)<1: ok+=1
    print(f"  estabilidade: {ok}/{N} perturbacoes de +-10% em cada parametro seguem viaveis (L<1)")
