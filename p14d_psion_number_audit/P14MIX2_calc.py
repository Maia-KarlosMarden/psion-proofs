import numpy as np
from scipy.optimize import differential_evolution
import P14MIX_calc as M
DIRS,MXD,pt=M.DIRS,M.MXD,M.pt
def mscore(Mn,ML,MS0,MSp,Mp):
    D=ML-Mn; q=D/MXD
    with np.errstate(all='ignore'):
        a=(MS0-ML)/D; b=(MSp-Mp)/D
        s=np.maximum.reduce([np.abs(a/0.437-1)/0.3,np.abs(b/1.425-1)/0.3,20*np.maximum(0,1-q/0.336)])
    return np.where(D>0,s,5.0).min()
def loss(x,model):
    o=pt(x)
    if o is None: return 50.0
    Mp=DIRS@M.STp; Mn=DIRS@M.STn; ML=DIRS@o["L"][1]; zL=abs(o["L"][0]-0.15)/0.10
    zsp=lambda r: abs(r-0.975)/0.175
    if model=="MIXa": return max(zL,zsp(1.103),mscore(Mn,ML,DIRS@M.STs0,DIRS@M.STsp1_oct,Mp))
    if model=="MIXb": return max(zL,zsp(o["Sp2"][0]),mscore(Mn,ML,DIRS@M.STs0,DIRS@o["Sp2"][1],Mp))
    zs0=abs(o["S0"][0])/0.15; best=1e9
    for s in ("Sp1","Sp2"): best=min(best,max(zsp(o[s][0]),mscore(Mn,ML,DIRS@o["S0"][1],DIRS@o[s][1],Mp)))
    return max(zL,zs0,best)
if __name__=="__main__":
    OPT={"O1":[0.258,0.224,0.047,0.229],"O2":[0.080,0.439,0.0,0.083],"O3":[0.035,0.456,0.095,0.0]}
    for k,z in OPT.items(): print("controle HYB",k,"loss=%.3f"%loss(np.array(z),"HYB"))
    for model in ("HYB","MIXa","MIXb"):
        r=differential_evolution(loss,[(0,6),(0,6),(0,20),(0,20)],args=(model,),popsize=15,maxiter=40,seed=72,tol=1e-6,workers=-1,updating='deferred',polish=False)
        print(f"{model}: loss minima = {r.fun:.3f} em (t,s,k,l)={np.round(r.x,3)} -> {'VIAVEL' if r.fun<1 else 'inviavel'}")
