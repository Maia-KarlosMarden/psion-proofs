import io, contextlib, runpy, numpy as np
from multiprocessing import Pool
from scipy.stats import qmc
with contextlib.redirect_stdout(io.StringIO()): O=runpy.run_path('P14OOS_calc.py',run_name='lib')
full3,extra,checks,Lrad,DIRS,MX=O['full3'],O['extra'],O['checks'],O['A']['Lrad'],O['DIRS'],O['MX']
H3=O['H3']; gse=O['gse']; DY=O['DY']; static=O['static']; rp=O['rp']
def work(z):
    o=full3(np.array(z))
    if o is None: return None
    R,Av=o; Lr=Lrad(R)
    if Lr>=1: return (0,0,0,0,0,0,0)
    ex=extra(np.array(z))
    if ex is None: return (0,0,0,0,0,0,0)
    R,Av,Rsm=ex
    v0,pa,pb,_,_=checks(Av,Rsm,DIRS)
    if not v0.any(): return (1,0,0,0,0,0,0)
    n0=int(v0.sum()); return (1,1,n0,int((v0&pa&pb).sum()),int((v0&pa).sum()),int((v0&pb).sum()),int((v0&pa&pb).any()))
def run(box,seed,label):
    N=16000; smp=qmc.LatinHypercube(d=4,seed=seed).random(N); lo=np.zeros(4); hi=np.array(box)
    X=lo+smp*(hi-lo)
    with Pool() as P: res=P.map(work,[tuple(x) for x in X],chunksize=40)
    res=np.array([r for r in res if r is not None])
    nr=int(res[:,0].sum()); nv=int(res[:,1].sum()); npairs=int(res[:,2].sum()); both=int(res[:,3].sum()); pa=int(res[:,4].sum()); pb=int(res[:,5].sum()); pp=int(res[:,6].sum())
    print(f"{label}: {N} pontos; radios L_r<1: {nr}; com V0 (>=1 direcao): {nv}; pares V0: {npairs}")
    if nv: print(f"   p_ponto (algum par passa P-A&P-B) = {pp}/{nv} = {pp/nv:.3f}; p_par = {both}/{npairs} = {both/max(npairs,1):.3f}; so P-A: {pa/max(npairs,1):.3f}; so P-B: {pb/max(npairs,1):.3f}")
if __name__=="__main__":
    run([0.5,0.8,0.5,0.6],51,"caixa 1")
    run([1.0,1.6,1.0,1.2],52,"caixa 2")
