import numpy as np, io, contextlib, runpy
from scipy.optimize import differential_evolution
import P14MIX_calc as M
with contextlib.redirect_stdout(io.StringIO()): H2=runpy.run_path('P14HYB2_calc.py',run_name='lib')
DY=H2['DY']; g=H2['gse']; DY['AR']['X']=DY['diag_arr']((3,4,5)); DY['AR']['Sm']=DY['diag_arr']((0,1,2))
DIRS,MXD=M.DIRS,M.MXD; RP=g("p",0,0,0,0)[0]
def gs(h,*x):
    try: return g(h,*x)
    except Exception: return [np.nan]*5
def pts(x):
    r={h:gs(h,*x) for h in ("L","S0","Sp1","Sp2","Sm","X")}
    if any(np.isnan(v[0]) for v in r.values()): return None
    return {h:(r[h][0]/RP,np.array(r[h][2:5])) for h in r}
def mscore(Mn,ML,MS0,MSp,Mp):
    D=ML-Mn; q=D/MXD
    with np.errstate(all='ignore'):
        a=(MS0-ML)/D; b=(MSp-Mp)/D
        s=np.maximum.reduce([np.abs(a/0.437-1)/0.3,np.abs(b/1.425-1)/0.3,20*np.maximum(0,1-q/0.336)])
    return np.where(D>0,s,5.0)
zsp=lambda r: abs(r-0.975)/0.175
def terms(x,model):
    o=pts(x)
    if o is None: return None
    Mp=DIRS@M.STp; Mn=DIRS@M.STn; ML=DIRS@o["L"][1]; zL=abs(o["L"][0]-0.15)/0.10
    cands=[]
    if model=="MIXa": cands=[(max(zL,zsp(1.103)),DIRS@M.STs0,DIRS@M.STsp1_oct)]
    elif model=="MIXb": cands=[(max(zL,zsp(o["Sp2"][0])),DIRS@M.STs0,DIRS@o["Sp2"][1])]
    else:
        zs0=abs(o["S0"][0])/0.15
        cands=[(max(zL,zs0,zsp(o[s][0])),DIRS@o["S0"][1],DIRS@o[s][1]) for s in ("Sp1","Sp2")]
    return o,Mn,ML,Mp,cands
def loss(x,model):
    t=terms(x,model)
    if t is None: return 50.0
    o,Mn,ML,Mp,c=t; return min(max(zr,mscore(Mn,ML,MS0,MSp,Mp).min()) for zr,MS0,MSp in c)
def oos(x,model):
    o,Mn,ML,Mp,c=terms(x,model); ok=False
    for zr,MS0,MSp in c:
        if zr>=1: continue
        v0=mscore(Mn,ML,MS0,MSp,Mp)<1
        D=ML-Mn; u=176.1/D
        msm=(DIRS@o["Sm"][1]-Mn)*u; mx=(DIRS@o["X"][1]-Mn)*u; Rsm=o["Sm"][0]
        pa=(0.65<=Rsm<=1.07)&(msm>=181)&(msm<=335); pb=(mx>=205)&(mx<=381)
        ok|=bool((v0&pa&pb).any())
    return ok
if __name__=="__main__":
    for model in ("HYB","MIXa","MIXb"):
        v=[]
        for sd in range(72,82):
            r=differential_evolution(loss,[(0,6),(0,6),(0,20),(0,20)],args=(model,),popsize=15,maxiter=40,seed=sd,tol=1e-6,workers=-1,updating='deferred',polish=False)
            via=r.fun<1; p=oos(r.x,model) if via else None; v.append((r.fun,via,p))
            print(f"{model} seed {sd}: loss {r.fun:.3f} x={np.round(r.x,3)} viavel={via} OOS(P-A&P-B)={p}",flush=True)
        nv=sum(x[1] for x in v); npass=sum(1 for x in v if x[2])
        print(f"=> {model}: otimos viaveis {nv}/10; passam OOS {npass} (fracao {npass/max(nv,1):.2f})",flush=True)
