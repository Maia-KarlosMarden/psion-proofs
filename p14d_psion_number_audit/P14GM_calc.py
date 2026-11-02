import math, numpy as np
import P14SPIN_calc as S
import P14FF3_calc as F
Q9=S.Q9; sz=S.sz
def spin_expect(t,J=+1):
    isL=[k<=2 for k in t]
    H=sum(J*S.sdot(a,b) for a in range(3) for b in range(a+1,3) if isL[a]!=isL[b])
    w,v=np.linalg.eigh(H); g=np.where(w<w[0]+1e-9)[0]; V=v[:,g]
    Sz=V.conj().T@S.Stot[2]@V; ew,ev=np.linalg.eigh((Sz+Sz.conj().T)/2); psi=V@ev[:,np.argmax(ew)]
    return [float((psi.conj()@S.op(sz,a)@psi).real) for a in range(3)]
sizes={0:0.413,1:0.442,4:0.691,5:0.340}; a=0.940; hc=0.19733
res={}
for nm,t in (("p",(0,4,5)),("n",(0,1,5))):
    sg=spin_expect(t); lab,qq,d=F.geo(t); print(f"{nm}: <sigma_z> por psion {np.round(sg,3)} (rotulos {lab}); pesos q*<sigma> = {np.round(qq*np.array(sg),3)}; soma (momento ~) = {float((qq*np.array(sg)).sum()):+.4f}")
    res[nm]=(lab,qq,d,sg)
def GM(nm,q):
    lab,qq,d,sg=res[nm]; return sum(qi*s*math.exp(-q*q*sizes[l]**2/6)*F.j0(q*a*di) for l,qi,di,s in zip(lab,qq,d,sg))
mup=GM("p",0.0); mun=GM("n",0.0); print(f"mu_n/mu_p = {mun/mup:+.3f} (obs -0.685)")
bad={"p":0,"n":0}
print("Q2    p: G_M/G_M(0)/G_D | n: G_M/G_M(0)/G_D")
for Q2 in F.Q2P:
    q=math.sqrt(Q2)/hc; GD=(1+Q2/0.71)**-2
    rp=GM("p",q)/mup/GD; rn=GM("n",q)/mun/GD
    bad["p"]+=abs(rp-1)>0.15; bad["n"]+=abs(rn-1)>0.15
    print(f"{Q2:4.1f}  {rp:+.3f}             | {rn:+.3f}")
print("pontos fora de 15%:",bad)
print("veredito:", "SUPPORT" if bad["p"]==0 and bad["n"]==0 else ("KILL" if bad["p"]>=3 or bad["n"]>=3 else "INCONCLUSIVO"))
