import numpy as np, math, runpy, io, contextlib
from multiprocessing import Pool
from scipy.stats import qmc
Q9=np.array([1,-5,-5,10,4,4])/9
TR={"p":(0,4,5),"n":(0,1,5),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
def z_of(h,R):
    return {"n":abs(R+0.164)/0.0164,"S0":abs(R)/0.15,"L":abs(R-0.15)/0.10,"Sp1":abs(R-0.975)/0.175,"Sp2":abs(R-0.975)/0.175}[h]
def r2(X,q):
    Qt=q.sum(); c=X.mean(0); v=(q*((X-c)**2).sum(1)).sum(); return v if abs(Qt)<1e-12 else v/Qt
def static_ring(N,home):
    th=2*math.pi*np.array(home)/N; out={}
    for h,t in TR.items():
        q=Q9[list(t)]; X=np.stack([np.cos(th[list(t)]),np.sin(th[list(t)])],1); out[h]=r2(X,q)
    return {h:out[h]/out["p"] for h in TR if h!="p"}
with contextlib.redirect_stdout(io.StringIO()):
    H2=runpy.run_path('P14HYB2_calc.py',run_name='lib'); OC=runpy.run_path('P14OCT_calc.py',run_name='lib'); TWm=runpy.run_path('P14TWIST_calc.py',run_name='lib')
g12=H2['gse']; g8=OC['gse']; gtw=TWm['gse']; arr_tw=TWm['arrays']
def dyn_row(args):
    G,x=args; res={}
    try:
        if G=="D12": f=lambda h: g12(h,*x)
        elif G=="D8": f=lambda h: g8(h,*x)
        else:
            t,s,k,l,hh,dl=x; AR={n:arr_tw(tr,hh,dl) for n,tr in TR.items()}; f=lambda h: gtw(h,t,s,k,l,AR)
        rp=f("p")[0]
        if np.isnan(rp) or abs(rp)<1e-9: return None
        for h in ("n","S0","L","Sp1","Sp2"):
            r=f(h)[0]
            if np.isnan(r): return None
            res[h]=z_of(h,r/rp)
    except Exception: return None
    return res
def pts(G,N=3000):
    d=6 if G=="TW" else 4; smp=qmc.LatinHypercube(d=d,seed=61).random(N); rng=np.random.default_rng(61)
    zt=rng.random(N)<0.05; zs=rng.random(N)<0.05; lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=[np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20]
    if G=="TW": X+= [0.2+2.8*smp[:,4],2*math.pi*smp[:,5]]
    return list(np.stack(X,1))
if __name__=="__main__":
    H=("n","S0","L","Sp1","Sp2"); best={}
    for nm,(N,home) in {"S6":(6,range(6)),"S8":(8,(0,1,2,4,5,6)),"S12":(12,tuple(2*k for k in range(6)))}.items():
        R=static_ring(N,home); best[nm]={h:z_of(h,R[h]) for h in H}
        print(nm,"static R:",{h:round(R[h],3) for h in H})
    for G in ("D8","D12","TW"):
        with Pool() as P: rows=[r for r in P.map(dyn_row,[(G,x) for x in pts(G)],chunksize=30) if r]
        best[G]={h:min(r[h] for r in rows) for h in H}; print(G,"valid points",len(rows))
    print("\nz-score minimo por composto e geometria (z<1 = compativel):")
    print("geom  "+"  ".join(f"{h:>6}" for h in H))
    for G,b in best.items(): print(f"{G:5} "+"  ".join(f"{b[h]:6.2f}" for h in H))
    print("\nmelhor geometria por composto:")
    for h in H:
        ok=[G for G in best if best[G][h]<1]; g=min(best,key=lambda G:best[G][h])
        print(f"{h:3}: melhor {g} (z={best[g][h]:.2f}); compativeis: {ok}")
