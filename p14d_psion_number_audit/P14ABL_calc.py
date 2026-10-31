import io, contextlib, runpy, math, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
with contextlib.redirect_stdout(io.StringIO()):
    TW=runpy.run_path('P14TWIST_calc.py',run_name='lib')
    DY=runpy.run_path('P14ESC_calc.py',run_name='lib')
rngD=np.random.default_rng(36); Dp=rngD.dirichlet([1,1,1],2000); DIRS=np.vstack([Dp,-Dp]); MX=np.abs(DIRS).max(1)
H_=["p","n","S0","L","Sp1","Sp2"]
SPEC={"V0":(1,1,0,1),"V1":(1,0,0,1),"V2":(0,1,0,1),"V3":(1,1,0,0),"V4":(0,0,1,1),"V5":(0,0,0,1)}   # (r1,r2,signs,scale)
def qvar(Av,spec):
    r1u,r2u,sg,sc=spec; M={h:DIRS@Av[h] for h in H_}; D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        ok=D>0
        if r1u: ok&=np.abs((M["S0"]-M["L"])/D/0.437-1)<0.3
        elif sg: ok&=(M["S0"]-M["L"])>0
        if r2u: ok&=(np.abs((M["Sp1"]-M["p"])/D/1.425-1)<0.3)|(np.abs((M["Sp2"]-M["p"])/D/1.425-1)<0.3)
        elif sg: ok&=((M["Sp1"]-M["p"])>0)|((M["Sp2"]-M["p"])>0)
    if not ok.any(): return None
    return float((D/MX)[ok].max())
def Lrad(R):
    return max(abs(R[0]+0.164)/0.0164,abs(R[1])/0.15,abs(R[2]-0.15)/0.10,min(abs(R[3]-0.975),abs(R[4]-0.975))/0.175)
def full_dy(x):   # x=(t,s,k,l) -> Lr, Av
    t,s,k,l=x; r={h:DY['gse'](h,t,s,k,l) for h in H_}
    if any(np.isnan(r[h][0]) for h in H_): return None
    rp=r["p"][0]; R=[r[h][0]/rp for h in ("n","S0","L","Sp1","Sp2")]
    return Lrad(R),{h:np.array(r[h][2:5]) for h in H_}
def full_tw(x):
    t,s,k,l,h,dl=x; AR={n:TW['arrays'](tr,h,dl) for n,tr in TW['TR'].items()}
    r={n:TW['gse'](n,t,s,k,l,AR) for n in H_}
    if any(np.isnan(r[n][0]) for n in H_): return None
    rp=r["p"][0]; R=[r[n][0]/rp for n in ("n","S0","L","Sp1","Sp2")]
    return Lrad(R),{n:np.array(r[n][2:5]) for n in H_}
MODELS={"DYN12":(full_dy,np.array([0,0,0,0.]),np.array([6,6,20,20.])),"TWIST":(full_tw,np.array([0,0,0,0,0.4,0.]),np.array([6,6,20,20,3.,math.pi]))}
def lossv(args,z=None):
    pass
def evalv(arg):
    model,spec_name,x=arg; f,lo,hi=MODELS[model]; out=f(np.clip(x,lo,hi))
    if out is None: return 50.0
    Lr,Av=out
    if spec_name=="V6": return Lr
    sp=SPEC[spec_name]; q=qvar(Av,sp)
    if q is None: return 40.0
    return Lr+20*max(0,1-q/0.336) if sp[3] else Lr
def table_row(arg):
    model,x=arg; f,lo,hi=MODELS[model]; out=f(x)
    if out is None: return [np.nan]*7
    Lr,Av=out; row=[Lr]
    for nm in ("V0","V1","V2","V3","V4","V5"):
        q=qvar(Av,SPEC[nm]); row.append(np.nan if q is None else q)
    return row
def nm(arg):
    model,spec_name,x0=arg; f,lo,hi=MODELS[model]
    r=minimize(lambda z: evalv((model,spec_name,np.clip(z,lo,hi))),np.clip(x0,lo,hi),method='Nelder-Mead',options={'maxiter':400,'xatol':1e-4,'fatol':1e-5}); return float(r.fun),np.clip(r.x,lo,hi)
if __name__=="__main__":
    from scipy.stats import qmc
    # DYN12 tabela SGN (lambda>0)
    d=np.load('P14SGN_table.npz'); Xd=d['X']; Xd=Xd[Xd[:,3]>0]
    # TWIST tabela (mesma amostragem de WP-TWIST)
    rng=np.random.default_rng(35); N=30000
    smp=qmc.LatinHypercube(d=6,seed=35).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u,a,b: np.exp(np.log(a)+u*(np.log(b)-np.log(a)))
    Xt=np.stack([np.where(zt,0,lg(smp[:,0],0.01,6)),np.where(zs,0,lg(smp[:,1],0.01,6)),smp[:,2]*20,smp[:,3]*20,lg(smp[:,4],0.4,3),smp[:,5]*math.pi],1)
    with Pool() as P:
        Rd=np.array(P.map(table_row,[("DYN12",x) for x in Xd],chunksize=100)); Rt=np.array(P.map(table_row,[("TWIST",x) for x in Xt],chunksize=100))
    for model,X_,Rw in (("DYN12",Xd,Rd),("TWIST",Xt,Rt)):
        Lr=Rw[:,0]; print(f"=== {model}: {len(X_)} pontos; raios puros (V6): tabela min L_r = {np.nanmin(Lr):.3f}")
        starts={}
        for j,nm_ in enumerate(("V0","V1","V2","V3","V4","V5"),1):
            q=Rw[:,j]; thr=0.336 if SPEC[nm_][3] else 1e-9
            ok=(q>=thr); v=np.where(ok,Lr,np.inf)
            print(f"  {nm_}: pontos que satisfazem {int(ok.sum())}; tabela min L_r = {np.nanmin(v):.3f}")
            sc=np.where(np.isnan(Lr),50.0,np.where(np.isnan(q),40.0,Lr+(20*np.maximum(0,1-np.nan_to_num(q)/0.336) if SPEC[nm_][3] else 0)))
            starts[nm_]=[X_[i] for i in np.argsort(sc)[:10]]
        starts["V6"]=[X_[i] for i in np.argsort(np.where(np.isnan(Lr),99,Lr))[:10]]
        jobs=[(model,nm_,x) for nm_,xs in starts.items() for x in xs]
        with Pool() as P: res=P.map(nm,jobs,chunksize=1)
        k=0
        for nm_,xs in starts.items():
            rr=res[k:k+len(xs)]; k+=len(xs); best=min(rr,key=lambda a:a[0])
            print(f"  [refino] {nm_}: min = {best[0]:.3f} ({'VIAVEL' if best[0]<1 else 'inviavel'}); viaveis {sum(1 for r in rr if r[0]<1)}/{len(rr)}; z={np.round(best[1],3)}")
