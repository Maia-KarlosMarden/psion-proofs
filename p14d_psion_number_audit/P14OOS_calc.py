import io, contextlib, runpy, numpy as np
from multiprocessing import Pool
from scipy.stats import qmc
with contextlib.redirect_stdout(io.StringIO()): H3=runpy.run_path('P14HYB3_calc.py',run_name='lib')
H2=runpy.run_path('P14HYB2_calc.py',run_name='lib'); DY=H2['DY']; gse=H2['gse']; static=H2['static']
DY['AR']['X']=DY['diag_arr']((3,4,5)); DY['AR']['Sm']=DY['diag_arr']((0,1,2))
A=H3['A']; DIRS,MX=A['DIRS'],A['MX']; qvar=A['qvar']; SPEC=A['SPEC']; full3=H3['full3']
rp=static("p")[0]
def extra(z):
    t,s,k,l=z; o=full3(z)
    sm=gse("Sm",t,s,k,l); xx=gse("X",t,s,k,l)
    if o is None or np.isnan(sm[0]) or np.isnan(xx[0]): return None
    R,Av=o; Av=dict(Av); Av["Sm"]=np.array(sm[2:5]); Av["X"]=np.array(xx[2:5]); return R,Av,sm[0]/rp
def checks(Av,Rsm,dirs,mask_scale=True):
    M={h:dirs@Av[h] for h in Av}; D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        r1=(M["S0"]-M["L"])/D; a=(M["Sp1"]-M["p"])/D; b=(M["Sp2"]-M["p"])/D; q=D/np.abs(dirs).max(1)
        v0=(D>0)&(np.abs(r1/0.437-1)<0.3)&((np.abs(a/1.425-1)<0.3)|(np.abs(b/1.425-1)<0.3))&(q>=0.336)
        u=176.1/D; msm=(M["Sm"]-M["n"])*u; mx=(M["X"]-M["n"])*u
        pa=(0.65<=Rsm)&(Rsm<=1.07)&(msm>=181)&(msm<=335); pb=(mx>=205)&(mx<=381)
    return v0,pa,pb,msm,mx
def rowbase(z):
    o=extra(z)
    if o is None: return (0,0,0)
    R,Av,Rsm=o; v0,pa,pb,_,_=checks(Av,Rsm,POS)
    return int(v0.sum()),int((v0&pa&pb).sum()),int((v0&pa).sum())
POS=H3['DIRS'][:2000]
if __name__=="__main__":
    OPT={"O1":[0.258,0.224,0.047,0.229],"O2":[0.080,0.439,0.0,0.083],"O3":[0.035,0.456,0.095,0.0]}
    npass=0
    for nm,z in OPT.items():
        z=np.array(z); R,Av,Rsm=extra(z); v0,pa,pb,msm,mx=checks(Av,Rsm,DIRS)
        n0=int(v0.sum()); both=v0&pa&pb
        print(f"{nm} {z}: R_Sm={Rsm:+.2f}; direcoes V0 aprovadas {n0}; Sigma- passa (P-A) em {int((v0&pa).sum())}; (3,4,5) passa (P-B) em {int((v0&pb).sum())}; ambos em {int(both.sum())}")
        print(f"   massa rel. n nas aprovadas: Sigma- [{msm[v0].min():.0f}, {msm[v0].max():.0f}] (258); (3,4,5) [{mx[v0].min():.0f}, {mx[v0].max():.0f}] (293)")
        npass+= bool(both.any())
    print(f"otimos que passam em (P-A) E (P-B): {npass}/3 -> {'SUCESSO GLOBAL' if npass>=2 else 'sem sucesso global'}")
    rng=np.random.default_rng(47); N=4000; smp=qmc.LatinHypercube(d=4,seed=47).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: res=np.array(P.map(rowbase,list(X),chunksize=40))
    tot=res[:,0].sum(); both=res[:,1].sum(); pa=res[:,2].sum()
    print(f"linha de base: {int((res[:,0]>0).sum())} pontos com alguma direcao V0; pares (ponto,direcao) V0: {int(tot)}; passam (P-A)&(P-B): {int(both)} => p0 = {both/max(tot,1):.3f}; passam so (P-A): {pa/max(tot,1):.3f}")
