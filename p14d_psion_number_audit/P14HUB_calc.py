import math, numpy as np
from scipy.optimize import brentq
Q9=np.array([1,-5,-5,10,4,4])/9
def mom(pts,q):
    X=np.array(pts); c=X.mean(0); v=float(np.sum(q*np.sum((X-c)**2,1))); Q=q.sum(); return v/Q if abs(Q)>1e-12 else v
print("Parte 1: identidade do arco <r^2>/a^2 = 7/9 Q_tot - 2/3 Q_c")
mx=0
for c in range(6):
    t=[(c-1)%6,c,(c+1)%6]; pts=[[math.cos(math.pi*k/3),math.sin(math.pi*k/3)] for k in t]; q=Q9[t]
    Q=q.sum(); direto=float(np.sum(q*np.sum((np.array(pts)-np.array(pts).mean(0))**2,1)))
    ident=7/9*Q-2/3*Q9[c]; mx=max(mx,abs(direto-ident))
print(f"   max |direto - identidade| = {mx:.1e}")
for nm,t,c in (("p (0,4,5) centro 5",(0,4,5),5),("n (0,1,5) centro 0",(0,1,5),0),("Sigma- (0,1,2) centro 1",(0,1,2),1),("Lambda (1,2,3) centro 2",(1,2,3),2),("(2,3,4) centro 3",(2,3,4),3)):
    print(f"   {nm}: {7/9*Q9[list(t)].sum()-2/3*Q9[c]:+.4f}")
print("   nao-arcos (1,3,5) e (0,2,4) sao triangulos equilateros: <r^2>/a^2 = Q_tot * (1) -> neutra = 0 exato (Sigma0=0)")
# Parte 2
CH=lambda k:'L' if k<=2 else 'R'
def hub(t):
    L=[k for k in t if CH(k)=='L']; R=[k for k in t if CH(k)=='R']
    if len(L)==2 and len(R)==1: return R[0]
    if len(L)==1 and len(R)==2: return L[0]
    return None
def r2(t,phi):
    h=hub(t)
    if h is None: return None
    ends=[k for k in t if k!=h]; b=1.0; s=math.radians(phi)/2
    pts={h:[0.0,0.0],ends[0]:[b*math.cos(s),b*math.sin(s)],ends[1]:[b*math.cos(s),-b*math.sin(s)]}
    return mom([pts[k] for k in t],Q9[list(t)])
rp=lambda phi:r2((0,4,5),phi)
Rn=lambda phi:r2((0,1,5),phi)/rp(phi)
phi0=brentq(lambda f:Rn(f)+0.164,61,179)
print(f"Parte 2: phi de calibracao (R_n = -0,164) = {phi0:.2f} graus; b = {math.sqrt(0.7071/rp(phi0)):.3f} fm")
print("   previsoes fora da amostra:")
res={}
for nm,t,tgt in (("Sigma0 (0,2,4)",(0,2,4),"|R|<=0,15"),("Lambda (1,2,3)",(1,2,3),"(0;0,25]"),("Sigma+ (1,3,5)",(1,3,5),"[0,80;1,15]"),("Sigma+ (2,3,4)",(2,3,4),"[0,80;1,15]")):
    R=r2(t,phi0)/rp(phi0); res[nm]=R
    ok={"|R|<=0,15":abs(R)<=0.15,"(0;0,25]":0<R<=0.25,"[0,80;1,15]":0.8<=R<=1.15}[tgt]
    print(f"     {nm}: hub {hub(t)}  R = {R:+.3f}  alvo {tgt} -> {'passa' if ok else 'falha'}")
print("   inertes (sem hub):",[t for t in ((0,1,2),(3,4,5)) if hub(t) is None])
print("   varredura de phi (informativa):")
for f in (62,75,90,105,120,150,180):
    print(f"     phi={f:3d}: R_n={Rn(f):+.3f} Sigma0={r2((0,2,4),f)/rp(f):+.3f} Lambda={r2((1,2,3),f)/rp(f):+.3f} Sigma+={r2((2,3,4),f)/rp(f):+.3f}")
