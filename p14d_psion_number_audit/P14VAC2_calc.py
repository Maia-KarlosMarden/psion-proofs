import io, contextlib, math, runpy, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    N_=runpy.run_path('P14SIM_network.py',run_name='lib')
one,OJS,J,PL,PR,V3=N_['one'],N_['OJS'],N_['J'],N_['PL'],N_['PR'],N_['V3']
mk=np.array([k-2.5 for k in range(6)]); Jp=J[0]+1j*J[1]
X=one(Jp@PL)@one(Jp@PR); H0d=one(np.diag(np.abs(mk)).astype(complex))
bond=np.zeros((6,6),complex); bond[2,3]=bond[3,2]=0.5; Bond=one(bond)
H=lambda th,d=0.0: H0d+OJS+(1/6)*(np.exp(1j*th)*X+np.exp(-1j*th)*X.conj().T)+d*Bond
import sys
SECT=sys.argv[1] if len(sys.argv)>1 else 'fock'
proj=lambda M: V3.conj()@M@V3.T
E0=lambda th,d=0.0: np.linalg.eigvalsh(proj(H(th,d)) if SECT=='triples' else H(th,d))[0]
th=2*math.pi*np.arange(96)/96; EA=np.array([E0(t) for t in th])
print(f"A: E0(theta) delta=0: min {EA.min():.12f} max {EA.max():.12f} spread {EA.max()-EA.min():.2e}")
EB=np.array([E0(2*math.pi*k/6) for k in range(6)])
C=lambda s: sum(mk[k]*mk[(k+s)%6] for k in range(6)); dC=[C(0)-C(s) for s in range(6)]
print("B: E0(theta_k):",np.round(EB,12),f"spread {EB.max()-EB.min():.2e}; DeltaC(s)=",[round(x,2) for x in dC])
for d in (0.02,0.05,0.1):
    E=np.array([E0(t,d) for t in th]); F=np.abs(np.fft.rfft(E))/96*2
    print(f"C: delta={d}: E0 range {E.max()-E.min():.3e}; harmonicas n=1..6:",[f"{F[n]:.2e}" for n in range(1,7)],f"a6/a1={F[6]/F[1]:.1e}")
