import io, contextlib, runpy, math, numpy as np
with contextlib.redirect_stdout(io.StringIO()): N=runpy.run_path('P14SIM_network.py',run_name='lib')
V3,OJS,TR=N['V3'],N['OJS'],N['TR']
Ot=(V3.conj()@OJS@V3.T); Ot=0.5*(Ot+Ot.conj().T)
mk=[k-2.5 for k in range(6)]; D0=np.diag([sum(abs(mk[k]) for k in t) for t in TR]).astype(complex)
Q9=np.array([1,-5,-5,10,4,4])/9
hexa={k:np.array([math.cos(math.pi*k/3),math.sin(math.pi*k/3)]) for k in range(6)}
def r2(t):
    Y=np.array([hexa[k] for k in t]); q=Q9[list(t)]; c=Y.mean(0); v=float(np.sum(q*np.sum((Y-c)**2,1))); Q=q.sum(); return v/Q if abs(Q)>1e-12 else v
R2=np.array([r2(t) for t in TR])
# verificacao estrutural
leak=max(abs(Ot[i,j]) for i in range(20) for j in range(20) if sum(TR[i])!=sum(TR[j]))
print(f"verificacao: max |O_JS| entre setores de sum k diferentes = {leak:.1e}")
SEC={6:[(0,1,5),(0,2,4),(1,2,3)],9:[(0,4,5),(1,3,5),(2,3,4)]}
def dressed(g, steps=400):
    # continuidade adiabatica (regra registrada): seguir autovetores de 0 ate g por sobreposicao
    out={}
    for s,ts in SEC.items():
        idx=[TR.index(t) for t in ts]
        prev=np.eye(3); lab=list(ts)
        for gg in np.linspace(0,g,steps+1)[1:] if g>0 else []:
            w,v=np.linalg.eigh((D0+gg*Ot)[np.ix_(idx,idx)])
            ov=np.abs(prev.conj().T@v)  # linhas: rotulo anterior, colunas: novo
            perm=[int(np.argmax(ov[i])) for i in range(3)]
            assert len(set(perm))==3, "ambiguidade adiabatica"
            prev=v[:,perm]
        H=(D0+g*Ot)[np.ix_(idx,idx)]
        for i,t in enumerate(ts):
            P=np.abs(prev[:,i])**2; E=float(np.real(prev[:,i].conj()@H@prev[:,i]))
            out[t]=(float(P@R2[idx]),E,P)
        out[f"uniq{s}"]=True
    out[(0,1,2)]=(R2[TR.index((0,1,2))],None,None)
    return out
def ratios(g):
    d=dressed(g); rp=d[(0,4,5)][0]
    return {"n":d[(0,1,5)][0]/rp,"Sm":d[(0,1,2)][0]/rp,"S0":d[(0,2,4)][0]/rp,"L":d[(1,2,3)][0]/rp,"Sp1":d[(1,3,5)][0]/rp,"Sp2":d[(2,3,4)][0]/rp,"u":d["uniq6"] and d["uniq9"],"d":d}
def crit(R):
    core=abs(R["n"]/-0.164-1)<0.1 and abs(R["Sm"]-0.86)<0.21
    sp=(0.8<=R["Sp1"]<=1.15) or (0.8<=R["Sp2"]<=1.15)
    return core, core and abs(R["S0"])<=0.15 and 0<R["L"]<=0.25 and sp
for g in (0.0,1.0,math.sqrt(2*math.pi**2/3)):
    R=ratios(g); c,e=crit(R)
    print(f"g={g:.3f}: n={R['n']:+.3f} Sm={R['Sm']:+.3f} S0={R['S0']:+.3f} L={R['L']:+.3f} Sp(1,3,5)={R['Sp1']:+.3f} Sp(2,3,4)={R['Sp2']:+.3f} rotulos unicos={R['u']} nucleo={c} estendido={e}")
    for lab in [(0,1,5),(0,2,4),(1,2,3),(0,4,5),(1,3,5),(2,3,4)]:
        _,w,P=R["d"][lab]; print(f"     {lab}: E={w:+.3f} pesos={np.round(P,3).tolist()}")
print("varredura g:")
gs=np.linspace(0,4,401); rows=[(g,ratios(g)) for g in gs]
for g,R in rows[::40]: print(f"  g={g:.1f} n={R['n']:+.3f} L={R['L']:+.3f} S0={R['S0']:+.3f} Sp1={R['Sp1']:+.3f} Sp2={R['Sp2']:+.3f}")
cg=[g for g,R in rows if crit(R)[0]]; eg=[g for g,R in rows if crit(R)[1]]
print(f"g que passam nucleo: {len(cg)} ({min(cg) if cg else '-'}..{max(cg) if cg else '-'}); estendido: {len(eg)} ({min(eg) if eg else '-'}..{max(eg) if eg else '-'})")
