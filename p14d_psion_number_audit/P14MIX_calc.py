import io, contextlib, runpy, numpy as np
from multiprocessing import Pool
from scipy.stats import qmc
with contextlib.redirect_stdout(io.StringIO()):
    H3=runpy.run_path('P14HYB3_calc.py',run_name='lib'); H2=runpy.run_path('P14HYB2_calc.py',run_name='lib'); OC=runpy.run_path('P14OCT_calc.py',run_name='lib')
g12=H2['gse']; g8=OC['gse']; DIRS=H3['A']['DIRS']; MXD=np.abs(DIRS).max(1)
S=lambda g,h: np.array(g(h,0,0,0,0)[2:5])
STp,STn,STs0=S(g12,"p"),S(g12,"n"),S(g12,"S0"); STsp1_oct=S(g8,"Sp1")
def pt(x):
    r={h:g12(h,*x) for h in ("p","L","S0","Sp1","Sp2")}
    if any(np.isnan(v[0]) for v in r.values()) or abs(r["p"][0])<1e-9: return None
    rp=r["p"][0]; return {h:(r[h][0]/rp,np.array(r[h][2:5])) for h in r}
def v0(M,Msp,Mp):
    D=M["L"]-M["n"]; q=D/MXD
    with np.errstate(all='ignore'):
        a=(M["S0"]-M["L"])/D; b=(Msp-Mp)/D
        return (D>0)&(np.abs(a/0.437-1)<0.3)&(np.abs(b/1.425-1)<0.3)&(q>=0.336)
def row(x):
    o=pt(x)
    if o is None: return (0,0,0,0,0,0)
    M={"n":DIRS@STn,"S0s":DIRS@STs0,"L":DIRS@o["L"][1],"S0d":DIRS@o["S0"][1]}; Mp=DIRS@STp
    out=[]
    # MIXa
    ok=0
    if 0.05<=o["L"][0]<=0.25 and 0.80<=1.103<=1.15:
        m=dict(M); m["S0"]=m["S0s"]; ok=v0(m,DIRS@STsp1_oct,Mp)
    out+=[int(np.sum(ok))]
    # MIXb
    ok=0
    if 0.05<=o["L"][0]<=0.25 and 0.80<=o["Sp2"][0]<=1.15:
        m=dict(M); m["S0"]=m["S0s"]; ok=v0(m,DIRS@o["Sp2"][1],Mp)
    out+=[int(np.sum(ok))]
    # HYB
    tot=0
    if 0.05<=o["L"][0]<=0.25 and abs(o["S0"][0])<=0.15:
        m=dict(M); m["S0"]=m["S0d"]; ok=np.zeros(len(DIRS),bool)
        for s in ("Sp1","Sp2"):
            if 0.80<=o[s][0]<=1.15: ok|=v0(m,DIRS@o[s][1],Mp)
        tot=int(ok.sum())
    out+=[tot]
    return tuple([out[0]>0,out[1]>0,out[2]>0]+out)
if __name__=="__main__":
    N=4000; rng=np.random.default_rng(71); smp=qmc.LatinHypercube(d=4,seed=71).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: res=np.array(P.map(row,list(X),chunksize=40),dtype=float)
    print("static <A,B,C>: p",STp.round(3),"n",STn.round(3),"S0",STs0.round(3),"Sp1(oct)",STsp1_oct.round(3))
    for i,nm in enumerate(("MIXa (S0 hex, L dyn, S+(1,3,5) oct)","MIXb (S0 hex, L & S+(2,3,4) dyn)","HYB (S0,L,S+ dyn; p,n static)")):
        print(f"{nm}: pontos viaveis {int(res[:,i].sum())}; pares (ponto,direcao) {int(res[:,3+i].sum())}")
