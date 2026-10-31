import io, contextlib, runpy, math, numpy as np
with contextlib.redirect_stdout(io.StringIO()): C=runpy.run_path('P14GEO10_calc.py',run_name='x')
TR,D0,Ot,R2,SEC=C['TR'],C['D0'],C['Ot'],C['R2'],C['SEC']
Q9=np.array([1,-5,-5,10,4,4])/9; mk=np.arange(6)-2.5
MU=np.array([sum(Q9[k]*mk[k] for k in t) for t in TR])
meas={"n":-0.685,"L":-0.2195,"Sp":0.880,"Sm":-0.415}
def states(g):
    S={}
    for s,ts in SEC.items():
        idx=[TR.index(t) for t in ts]; w,v=np.linalg.eigh((D0+g*Ot)[np.ix_(idx,idx)]); P=np.abs(v)**2
        S[s]=[(w[e],P[:,e]@R2[idx],P[:,e]@MU[idx],P[:,e]) for e in range(3)]
    i=TR.index((0,1,2)); S["Sm"]=(D0[i,i].real+g*Ot[i,i].real,R2[i],MU[i])
    return S
def evalg(g):
    S=states(g); p,Sp=S[9][0],S[9][1]; n,L,S0=S[6][0],S[6][1],S[6][2]; Sm=S["Sm"]
    mu={"n":n[2]/p[2],"L":L[2]/p[2],"Sp":Sp[2]/p[2],"Sm":Sm[2]/p[2]}
    o1=all(np.sign(mu[k])==np.sign(meas[k]) and abs(mu[k]/meas[k]-1)<0.3 for k in meas)
    o2=n[0]>p[0]; o3=(S0[0]-n[0])/(L[0]-n[0]); o4=(Sp[0]-p[0])/(L[0]-n[0])
    return dict(mu=mu,mu_p=p[2],o1=o1,o2=o2,En_Ep=n[0]-p[0],o3=o3,o3ok=abs(o3/1.437-1)<0.2,o4=o4,o4ok=abs(o4/1.426-1)<0.2,
                Rn=n[1]/p[1],RSm=Sm[1]/p[1],RL=L[1]/p[1],RS0=S0[1]/p[1],RSp=Sp[1]/p[1])
for g in (1.0,math.sqrt(2*math.pi**2/3)):
    r=evalg(g)
    print(f"g={g:.3f}")
    print("  O1 mu/mu_p: "+", ".join(f"{k}={r['mu'][k]:+.3f} (med {meas[k]:+.3f})" for k in meas)+f" -> {'PASSA' if r['o1'] else 'falha'}   [mu_p={r['mu_p']:+.3f}]")
    print(f"  O2 E_n - E_p = {r['En_Ep']:+.3f} -> {'PASSA' if r['o2'] else 'falha'}")
    print(f"  O3 gap Sigma0/Lambda = {r['o3']:.3f} (med 1,437) -> {'PASSA' if r['o3ok'] else 'falha'}")
    print(f"  O4 gap Sigma+/Lambda = {r['o4']:.3f} (med 1,426) -> {'PASSA' if r['o4ok'] else 'falha'}")
    print(f"  ja vistos: R_n={r['Rn']:+.3f} R_Sm={r['RSm']:+.3f} R_L={r['RL']:+.3f} R_S0={r['RS0']:+.3f} R_Sp={r['RSp']:+.3f}")
print("varredura informativa:")
for g in (0.5,0.75,1,1.5,2,3,4,5):
    r=evalg(g); print(f"  g={g:4.2f} mu n/L/Sp/Sm = {r['mu']['n']:+.2f}/{r['mu']['L']:+.2f}/{r['mu']['Sp']:+.2f}/{r['mu']['Sm']:+.2f}  En-Ep={r['En_Ep']:+.2f}  O3={r['o3']:.2f} O4={r['o4']:.2f}")
