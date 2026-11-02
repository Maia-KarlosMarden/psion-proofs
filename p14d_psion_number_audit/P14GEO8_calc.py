import math, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
def r2(t,X):
    Y=np.array([X[k] for k in t]); q=Q9[list(t)]; c=Y.mean(0)
    v=float(np.sum(q*np.sum((Y-c)**2,1))); Q=q.sum(); return v/Q if abs(Q)>1e-12 else v
hexa={k:np.array([math.cos(math.pi*k/3),math.sin(math.pi*k/3)]) for k in range(6)}
T={"p":(0,4,5),"n":(0,1,5),"Sigma-":(0,1,2),"Lambda":(1,2,3),"Sigma0":(0,2,4),"Sigma+ (1,3,5)":(1,3,5),"Sigma+ (2,3,4)":(2,3,4),"(3,4,5) Q=+2":(3,4,5)}
rp=r2(T["p"],hexa); a=math.sqrt(0.7071/rp)
print(f"A) hexagono (a = {a:.3f} fm; r_p^2 fixado = 0,7071 fm^2):")
for k,t in T.items(): R=r2(t,hexa)/rp; print(f"   {k:16s} R = {R:+.4f}  <r^2> = {R*0.7071:+.3f} fm^2")
okS=lambda R:0.80<=R<=1.15; okL=lambda R:0<R<=0.25
print("B) psion 3 livre:")
for cand in ((1,3,5),(2,3,4)):
    good=[]; 
    for r in np.linspace(0,2,81):
        for ph in np.arange(0,360,3):
            X=dict(hexa); X[3]=r*np.array([math.cos(math.radians(ph)),math.sin(math.radians(ph))])
            S=r2(cand,X)/rp; L=r2((1,2,3),X)/rp
            if okS(S) and okL(L): good.append((round(r,3),ph))
    X=dict(hexa); nat=(okS(r2(cand,X)/rp), okL(r2((1,2,3),X)/rp))
    rs=[g[0] for g in good]; ps=[g[1] for g in good]
    print(f"   Sigma+={cand}: area valida {len(good)}/{81*120} pontos"+(f"; r em [{min(rs):.2f},{max(rs):.2f}], phi em [{min(ps)},{max(ps)}]" if good else "")+f"; posicao natural (1,180): Sigma+ {'ok' if nat[0] else 'falha'}, Lambda {'ok' if nat[1] else 'falha'}")
