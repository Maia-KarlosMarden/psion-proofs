import math, itertools, numpy as np
import P14SPIN_calc as S
import P14FF3_calc as F
sz=S.sz; Q9=S.Q9; hc=0.19733
sizes={0:0.413,1:0.442,4:0.691,5:0.340}; a=0.940
isL=lambda k:k<=2
def chord(i,j): return 2*math.sin(math.pi*abs(i-j)/6)
def Hs(t,rule):
    H=0
    for x in range(3):
        for y in range(x+1,3):
            same=isL(t[x])==isL(t[y]); w=0.0
            if rule=="R1": w=1.0 if not same else 0
            elif rule=="R2": w=1.0
            elif rule=="R3": w=1.0 if same else 0
            elif rule=="R4": w=-1.0 if same else 0
            elif rule=="R5": w=1/chord(t[x],t[y]) if not same else 0
            elif rule=="R6": w=1/chord(t[x],t[y])
            elif rule=="R7": w=1.0 if {t[x],t[y]}=={2,3} else 0
            if w!=0: H=H+w*S.sdot(x,y)
    return H if not isinstance(H,int) else np.zeros((8,8),complex)
Sz=S.Stot[2].real
def state(t,rule):
    H=Hs(t,rule); w,v=np.linalg.eigh(H); g=np.where(w<w[0]+1e-9)[0]; V=v[:,g]
    # restrict to S_z=+1/2 inside the ground manifold
    M=V.conj().T@Sz@V; ew,ev=np.linalg.eigh((M+M.conj().T)/2)
    # sub-manifold with Sz=+1/2
    idx=np.where(abs(ew-0.5)<1e-8)[0]
    # spin-3/2-only ground (no Sz=+1/2 doublet component) is excluded
    # check total spin: S^2 expectation on those states
    S2=sum(Stot_c@Stot_c for Stot_c in S.Stot)
    if len(idx)==0: return None
    W=V@ev[:,idx]; s2=np.real(np.trace(W.conj().T@S2@W))/len(idx)
    if abs(s2-0.75)>1e-6 and len(idx)==1: return "spin32"
    sig=[]
    for k in range(3):
        O=S.op(S.sz,k); sig.append(float(np.real(np.trace(W.conj().T@O@W))/len(idx)))
    return sig
def mom(t,rule):
    s=state(t,rule)
    if s is None or s=="spin32": return None
    return float(sum(Q9[l]*x for l,x in zip(t,s))),s
def shapes(t,rule):
    s=state(t,rule)
    lab,qq,d=F.geo(t)
    def G(q): return sum(qi*x*math.exp(-q*q*sizes[l]**2/6)*F.j0(q*a*di) for l,qi,di,x in zip(lab,qq,d,s))
    g0=G(0.0)
    if abs(g0)<1e-9: return None
    out=[]
    for Q2 in F.Q2P:
        q=math.sqrt(Q2)/hc; out.append(G(q)/g0/(1+Q2/0.71)**-2)
    return np.array(out)
if __name__=="__main__":
    surv=[]
    for rule in ("R1","R2","R3","R4","R5","R6","R7"):
        mp=mom((0,4,5),rule); mn=mom((0,1,5),rule)
        if mp is None or mn is None or abs(mp[0])<1e-9: print(f"{rule}: estado de spin invalido (spin 3/2 ou momento nulo do proton): excluida/falha"); continue
        ratio=mn[0]/mp[0]
        spp=shapes((0,4,5),rule); spn=shapes((0,1,5),rule)
        okm=abs(ratio+0.685)<=0.20
        okp=spp is not None and np.all(np.abs(spp-1)<=0.15); okn=spn is not None and np.all(np.abs(spn-1)<=0.15)
        print(f"{rule}: mu_n/mu_p={ratio:+.3f} ({'ok' if okm else 'falha'}); forma p {'ok' if okp else 'falha'} (max desvio {np.abs(spp-1).max():.2f}); forma n {'ok' if okn else ('falha' if spn is not None else 'G_M(0)=0')}"+("" if spn is None else f" (max desvio {np.abs(spn-1).max():.2f})"))
        if okm and okp and okn: surv.append(rule)
    print("estagio 1 aprovadas:",surv)
    for rule in surv:
        mp=mom((0,4,5),rule)[0]; ml=mom((1,2,3),rule); msa=mom((1,3,5),rule); msb=mom((2,3,4),rule)
        rl=None if ml is None else ml[0]/mp; ra=None if msa is None else msa[0]/mp; rb=None if msb is None else msb[0]/mp
        okl=rl is not None and abs(rl+0.219)<=0.30; oks=any(r is not None and abs(r-0.880)<=0.30 for r in (ra,rb))
        print(f"{rule}: Lambda/p = {rl}, Sigma+/p = {ra}, {rb}; estagio 2: {'passa' if okl and oks else 'falha'}")
    print("sucesso da familia:", any(True for _ in surv) and 'ver estagio 2')
