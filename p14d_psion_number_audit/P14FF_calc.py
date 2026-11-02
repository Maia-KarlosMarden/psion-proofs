import math, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
def X(t): th=2*math.pi*np.array(t)/6; return np.stack([np.cos(th),np.sin(th)],1)
def parts(t):
    x=X(t); c=x.mean(0); d=np.linalg.norm(x-c,axis=1); return Q9[list(t)],d
j0=lambda z: np.where(np.abs(z)<1e-8,1.0,np.sin(z)/np.where(np.abs(z)<1e-8,1,z))
def F(t,q,a): qq,d=parts(t); return sum(qi*j0(q*a*di) for qi,di in zip(qq,d))
qp,dp=parts((0,4,5)); r2unit=sum(qi*di**2 for qi,di in zip(qp,dp))/qp.sum()
a=0.8414/math.sqrt(r2unit); print(f"<r2>_p (anel unitario) = {r2unit:.4f}; raio do anel a = {a:.4f} fm")
qn,dn=parts((0,1,5)); r2n=sum(qi*di**2 for qi,di in zip(qn,dn)); print(f"<r2>_n previsto = {r2n*a*a:+.4f} fm^2 (obs -0.1161)")
hc=0.19733; mN=0.9389; mun=-1.913
rows=[]
for Q2 in (0.1,0.2,0.3,0.5,0.8,1.0):
    q=math.sqrt(Q2)/hc; GD=(1+Q2/0.71)**-2; tau=Q2/(4*mN**2); GEn=-mun*tau/(1+5.6*tau)*GD
    rows.append((Q2,float(F((0,4,5),q,a)),GD,float(F((0,1,5),q,a)),GEn))
mxn=max(abs(r[4]) for r in rows)
print("Q2    F_p     G_D    F_p/G_D | F_n     G_En    |dF_n|/max|G_En|")
bad_p=[]; ok_n=True
for Q2,Fp,GD,Fn,GEn in rows:
    rp=Fp/GD; dn_=abs(Fn-GEn)/mxn; print(f"{Q2:4.1f} {Fp:+.4f} {GD:.4f} {rp:+.3f}   | {Fn:+.4f} {GEn:+.4f} {dn_:.2f}")
    if abs(rp-1)>0.15: bad_p.append(Q2)
    if dn_>0.15: ok_n=False
print("pontos do proton fora de 15%:",bad_p," | neutron dentro de 15% em todos:",ok_n)
nbad=sum(1 for q in bad_p if q>=0.5)
print("veredito:", "KILL" if nbad>=2 else ("SUPPORT" if (not bad_p and ok_n) else "inconclusivo"))
