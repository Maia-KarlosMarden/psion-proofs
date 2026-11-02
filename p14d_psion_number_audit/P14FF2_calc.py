import io, contextlib, math, runpy, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.optimize import differential_evolution
with contextlib.redirect_stdout(io.StringIO()): H2=runpy.run_path('P14HYB2_calc.py',run_name='lib')
DY=H2['DY']; AR=DY['AR']; T1=DY['T1']; S1=DY['S1']; PA=DY['PA']; XY=DY['XY']; Q9=DY['Q9']; TR=DY['TR']
def dists(t):
    X=XY[PA[:, :]]; X=XY[PA]; c=X.mean(1,keepdims=True); return np.linalg.norm(X-c,axis=2)   # (1320,3) distance to unweighted centroid
D={h:dists(TR[h]) for h in ("p","n")}; QQ={h:Q9[list(TR[h])] for h in ("p","n")}
hc=0.19733; mN=0.9389; mun=-1.913
def j0(z): z=np.asarray(z,float); out=np.ones_like(z); m=np.abs(z)>1e-8; out[m]=np.sin(z[m])/z[m]; return out
def prob(h,t,s,k,l):
    A,B,C,R2=AR[h]; dg=A+k*B+l*C
    if t<1e-9 and s<1e-9: 
        p=np.zeros(len(dg)); p[int(np.argmin(dg))]=1; return p
    Hm=sp.diags(dg)-t*T1-s*S1
    w,v=eigsh(Hm.tocsc(),k=1,which='SA'); return v[:,0]**2
def FF(h,p,a,q): 
    return float(sum(p@(QQ[h][i]*j0(q*a*D[h][:,i])) for i in range(3)))
def R2m(h,p): return float(p@AR[h][3])
Q2P=(0.1,0.2,0.3,0.5,0.8,1.0); Q2X=(0.4,0.6)
def curves(x,Q2s):
    t,s,k,l=x; pp=prob("p",*x); pn=prob("n",*x); rp=R2m("p",pp)
    if not np.isfinite(rp) or rp<=1e-9: return None
    a=0.8414/math.sqrt(rp); Rn=R2m("n",pn)/rp; out=[]
    for Q2 in Q2s:
        q=math.sqrt(Q2)/hc; GD=(1+Q2/0.71)**-2; tau=Q2/(4*mN**2); GEn=-mun*tau/(1+5.6*tau)*GD
        out.append((Q2,FF("p",pp,a,q),GD,FF("n",pn,a,q),GEn))
    return Rn,out
def loss(x,extra=False):
    try:
        o=curves(x,Q2P)
        if o is None: return 50.0
        Rn,rows=o; mx=max(abs(r[4]) for r in rows)
        lp=max(abs(r[1]/r[2]-1) for r in rows)/0.15; ln=max(abs(r[3]-r[4]) for r in rows)/(0.15*mx); lr=abs(Rn+0.164)/0.0164
        return max(lp,ln,lr)
    except Exception: return 50.0
def extra(x):
    o=curves(x,Q2X); Rn,rows=o; mx=max(abs(r[4]) for r in rows)
    return max(abs(r[1]/r[2]-1) for r in rows)/0.15, max(abs(r[3]-r[4]) for r in rows)/(0.15*mx)
if __name__=="__main__":
    print("referencia estatica L = 5.2 (WP-FF)")
    res=[]
    for sd in range(72,77):
        r=differential_evolution(loss,[(0,6),(0,6),(0,20),(0,20)],popsize=15,maxiter=40,seed=sd,tol=1e-6,workers=-1,updating='deferred',polish=False)
        via=r.fun<1; ex=extra(r.x) if r.fun<50 else None
        print(f"seed {sd}: L = {r.fun:.3f} x={np.round(r.x,3)} viavel={via}; pontos extras (Q2=0.4,0.6) perdas relativas: {None if ex is None else tuple(round(v,2) for v in ex)}",flush=True)
        res.append((r.fun,ex))
    nv=sum(1 for f,e in res if f<1); nx=sum(1 for f,e in res if f<1 and e and max(e)<1)
    print(f"=> viaveis {nv}/5; viaveis que tambem passam nos pontos extras {nx}")
    print("veredito:", "KILL" if nv==0 else ("SUPPORT" if nx>=2 else "inconclusivo"))
