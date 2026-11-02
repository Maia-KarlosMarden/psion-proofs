import io, contextlib, runpy, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
with contextlib.redirect_stdout(io.StringIO()): A=runpy.run_path('P14ABL_calc.py',run_name='lib')
evalv=A['evalv']; HI=np.array([6,6,20,20.]); LO=np.zeros(4)
def Lc(z): return evalv(("DYN12","V2",np.clip(z,LO,HI)))
U=np.random.default_rng(40).uniform(-1,1,(6,4))
def pert(z,u): return np.clip(z+0.10*np.maximum(np.abs(z),0.05)*u,LO,HI)
def Lrob(z):
    z=np.clip(z,LO,HI); return max(Lc(z),max(Lc(pert(z,u)) for u in U))
def job(z0):
    r=minimize(Lrob,np.clip(z0,LO,HI),method='Nelder-Mead',options={'maxiter':200,'xatol':1e-4,'fatol':1e-5}); return float(r.fun),np.clip(r.x,LO,HI)
def val(z):
    V=np.random.default_rng(41).uniform(-1,1,(100,4)); return float(np.mean([Lc(pert(z,u))<1 for u in V]))
if __name__=="__main__":
    d=np.load('P14SGN_table.npz'); X=d['X']; X=X[X[:,3]>0]
    with Pool() as P: sc=np.array(P.map(Lc,list(X),chunksize=100))
    starts=[X[i] for i in np.argsort(sc)[:8]]+[np.array([0.478,0.156,0,0.48]),np.array([0.349,0.226,0,0.411]),np.array([0.001,0.069,0.26,0.11])]
    print(f"tabela: melhor L_c nominal (V2) = {sc.min():.3f}; nominal L_c<1: {int((sc<1).sum())}")
    with Pool() as P: res=P.map(job,starts,chunksize=1)
    res.sort(key=lambda a:a[0]); print(f"min L_rob = {res[0][0]:.3f} ({'VIAVEL ROBUSTO' if res[0][0]<1 else 'inviavel'}); viaveis robustos {sum(1 for r in res if r[0]<1)}/{len(res)}")
    with Pool() as P: vals=P.map(val,[r[1] for r in res[:5]])
    for (f,z),v in zip(res[:5],vals):
        print(f"   L_rob={f:.3f} L_c nominal={Lc(z):.3f} (t,s,k,l)=({z[0]:.3f},{z[1]:.3f},{z[2]:.3f},{z[3]:.3f}); validacao +-10% (100 novas): {100*v:.0f}% com L_c<1")
