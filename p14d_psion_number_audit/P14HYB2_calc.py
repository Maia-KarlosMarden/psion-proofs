import io, contextlib, runpy, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from scipy.stats import qmc
with contextlib.redirect_stdout(io.StringIO()): A=runpy.run_path('P14ABL_calc.py',run_name='lib')
DY=A['DY']; AR=DY['AR']; PL=DY['PL']; gse=DY['gse']; qvar=A['qvar']; H_=A['H_']; H_=list(H_); SPEC=A['SPEC']; Lrad=A['Lrad']; TR=DY['TR']
DY['AR']['Sm']=DY['diag_arr']((0,1,2)); TR=dict(TR); TR['Sm']=(0,1,2)
NAT={h:PL.index(tuple(2*k for k in TR[h])) for h in list(H_)+['Sm']}
def static(h):
    a,b,c,r2=AR[h]; i=NAT[h]; return float(r2[i]),np.array([a[i],b[i],c[i]])
HI=np.array([6,6,20,20.]); LO=np.zeros(4)
def full(z):
    t,s,k,l=z
    rp,Ap=static("p"); rn,An=static("n")
    D_={h:gse(h,t,s,k,l) for h in ("S0","L","Sp1","Sp2","Sm")}
    if any(np.isnan(v[0]) for v in D_.values()): return None
    R=[rn/rp,D_["S0"][0]/rp,D_["L"][0]/rp,D_["Sp1"][0]/rp,D_["Sp2"][0]/rp]
    Av={"p":Ap,"n":An}
    for h in D_: Av[h]=np.array(D_[h][2:5])
    Rm=D_["Sm"][0]/rp
    return Lrad(R),Av,R+[Rm]
def Lc(z,ver):
    out=full(np.clip(z,LO,HI))
    if out is None: return 50.0
    Lr,Av,_=out; q=qvar(Av,SPEC[ver])
    return 40.0 if q is None else Lr+20*max(0,1-q/0.336)
U=np.random.default_rng(40).uniform(-1,1,(6,4))
def pert(z,u): return np.clip(z+0.10*np.maximum(np.abs(z),0.05)*u,LO,HI)
def Lrob(z,ver):
    z=np.clip(z,LO,HI); return max(Lc(z,ver),max(Lc(pert(z,u),ver) for u in U))
def tab(arg): z,ver=arg; return Lc(z,ver)
def job(arg):
    ver,z0=arg; r=minimize(lambda z:Lrob(z,ver),np.clip(z0,LO,HI),method='Nelder-Mead',options={'maxiter':60,'xatol':1e-4,'fatol':1e-5}); return ver,float(r.fun),np.clip(r.x,LO,HI)
def val(arg):
    ver,z=arg; V=np.random.default_rng(41).uniform(-1,1,(100,4)); return float(np.mean([Lc(pert(z,u),ver)<1 for u in V]))
if __name__=="__main__":
    rng=np.random.default_rng(43); N=8000; smp=qmc.LatinHypercube(d=4,seed=43).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    rn,_=static("n"); rp,_=static("p"); r0,_=static("S0")
    print(f"estaticos: apenas p, n; R_n = {rn/rp:+.3f}")
    with Pool() as P:
        out=P.map(full,list(X),chunksize=80)
        for ver in ("V0","V2"):
            L=np.array(P.map(tab,[(x,ver) for x in X],chunksize=80))
            Lr=np.array([np.nan if o is None else o[0] for o in out])
            print(f"[1] {ver}: tabela min L_c = {L.min():.3f}; L_c<1: {int((L<1).sum())}; raios isolados: min L_r = {np.nanmin(Lr):.3f}, L_r<1: {int(np.nansum(Lr<1))}")
            st=[(ver,X[i]) for i in np.argsort(L)[:6]]
            res=P.map(job,st,chunksize=1); res.sort(key=lambda a:a[1])
            vals=P.map(val,[(ver,r[2]) for r in res[:3]])
            print(f"[2] {ver}: min L_rob = {res[0][1]:.3f} ({'VIAVEL ROBUSTO' if res[0][1]<1 else 'inviavel'}); viaveis {sum(1 for r in res if r[1]<1)}/{len(res)}")
            for (v,f,z),vv in zip(res[:3],vals):
                o=full(z); print(f"   L_rob={f:.3f} (t,s,k,l)=({z[0]:.3f},{z[1]:.3f},{z[2]:.3f},{z[3]:.3f}) R_n={o[2][0]:+.3f} R_S0={o[2][1]:+.3f} R_L={o[2][2]:+.3f} R_Sp1={o[2][3]:+.3f} R_Sp2={o[2][4]:+.3f} R_Sm={o[2][5]:+.3f}; validacao +-10%: {100*vv:.0f}%")
