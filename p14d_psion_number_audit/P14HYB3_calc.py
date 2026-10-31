import io, contextlib, runpy, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from scipy.stats import qmc
with contextlib.redirect_stdout(io.StringIO()): H=runpy.run_path('P14HYB2_calc.py',run_name='lib')
A=H['A']; static=H['static']; gse=H['gse']; Lrad=A['Lrad']; qvar=A['qvar']; SPEC=A['SPEC']; DIRS,MX=A['DIRS'],A['MX']
HI=np.array([6,6,20,20.]); LO=np.zeros(4)
HAD=["p","n","S0","L","Sp1","Sp2"]
def full3(z):
    t,s,k,l=z; rp,Ap=static("p"); rn,An=static("n")
    D_={h:gse(h,t,s,k,l) for h in ("S0","L","Sp1","Sp2")}
    if any(np.isnan(v[0]) for v in D_.values()): return None
    R=[rn/rp]+[D_[h][0]/rp for h in ("S0","L","Sp1","Sp2")]
    Av={"p":Ap,"n":An}
    for h in D_: Av[h]=np.array(D_[h][2:5])
    return R,Av
def Lc(z):
    z=np.clip(z,LO,HI); o=full3(z)
    if o is None: return 50.0
    R,Av=o; q=qvar(Av,SPEC["V0"]); Lr=Lrad(R)
    return 40.0 if q is None else Lr+20*max(0,1-q/0.336)
U=np.random.default_rng(40).uniform(-1,1,(6,4))
def pert(z,u): return np.clip(z+0.10*np.maximum(np.abs(z),0.05)*u,LO,HI)
def Lrob(z):
    z=np.clip(z,LO,HI); return max(Lc(z),max(Lc(pert(z,u)) for u in U))
def row(z):
    o=full3(z)
    if o is None: return [np.nan]*20
    R,Av=o; return list(R)+sum([list(Av[h]) for h in HAD],[])
def job(z0):
    r=minimize(Lrob,np.clip(z0,LO,HI),method='Nelder-Mead',options={'maxiter':100,'xatol':1e-4,'fatol':1e-5}); return float(r.fun),np.clip(r.x,LO,HI)
def val(z):
    V=np.random.default_rng(41).uniform(-1,1,(100,4)); return float(np.mean([Lc(pert(z,u))<1 for u in V]))
def mass_q(Avs,r1t,r2t):
    M={h:DIRS@Avs[h] for h in HAD}; D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        ok=(D>0)&(np.abs((M["S0"]-M["L"])/D/r1t-1)<0.3)&((np.abs((M["Sp1"]-M["p"])/D/r2t-1)<0.3)|(np.abs((M["Sp2"]-M["p"])/D/r2t-1)<0.3))&(D/MX>=0.336)
    return bool(ok.any())
if __name__=="__main__":
    rng=np.random.default_rng(43); N=8000; smp=qmc.LatinHypercube(d=4,seed=43).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: T=np.array(P.map(row,list(X),chunksize=80))
    np.savez('P14HYB3_table.npz',X=X,T=T)
    ok=~np.isnan(T[:,0]); R=T[:,:5]; Av=lambda i:{h:T[i,5+3*j:8+3*j] for j,h in enumerate(HAD)}
    def Lr_t(c): 
        d=np.stack([np.abs(R[:,0]-c[0])/0.0164,np.abs(R[:,1]-c[1])/0.15,np.abs(R[:,2]-c[2])/0.10,np.minimum(np.abs(R[:,3]-c[3]),np.abs(R[:,4]-c[3]))/0.175],1); return np.where(np.isnan(R[:,0]),99,d.max(1))
    def hit(c,r1t,r2t):
        idx=np.where(Lr_t(c)<1)[0]
        return any(mass_q(Av(i),r1t,r2t) for i in idx)
    real_c=(-0.164,0.0,0.15,0.975); idx=np.where(Lr_t(real_c)<1)[0]; nreal=[i for i in idx if mass_q(Av(i),0.437,1.425)]
    print(f"tabela {N} pontos; alvo real: raios L_r<1 em {len(idx)} pontos; com massa V0: {len(nreal)} pontos")
    rs=np.random.default_rng(46); ha=hb=0; M_=300
    for _ in range(M_):
        c=rs.uniform(-0.6,2.2,4); r1t=rs.uniform(0.1,1.0); r2t=rs.uniform(0.5,3.0)
        ha+=hit(c,r1t,r2t); cb=(-0.164,c[1],c[2],c[3]); hb+=hit(cb,r1t,r2t)
    print(f"P_acaso(a) (todos aleatorios): {ha}/{M_} = {ha/M_:.3f}; P_acaso(b) (n fixo no real): {hb}/{M_} = {hb/M_:.3f}")
    st=[X[i] for i in np.argsort([Lc(x) for x in X[:0]])] if False else None
    with Pool() as P: L0=np.array(P.map(Lc,list(X),chunksize=80))
    print(f"tabela: min L_c = {L0.min():.3f}; L_c<1: {int((L0<1).sum())}")
    starts=[X[i] for i in np.argsort(L0)[:10]]
    with Pool() as P: res=P.map(job,starts,chunksize=1)
    res.sort(key=lambda a:a[0]); print(f"refino: min L_rob = {res[0][0]:.3f} ({'VIAVEL ROBUSTO' if res[0][0]<1 else 'inviavel'}); viaveis {sum(1 for r in res if r[0]<1)}/{len(res)}")
    with Pool() as P: vals=P.map(val,[r[1] for r in res[:4]])
    for (f,z),v in zip(res[:4],vals):
        o=full3(z); R_,Av_=o; print(f"   L_rob={f:.3f} (t,s,k,l)=({z[0]:.3f},{z[1]:.3f},{z[2]:.3f},{z[3]:.3f}) R={np.round(R_,3)}; validacao +-10%: {100*v:.0f}%")
