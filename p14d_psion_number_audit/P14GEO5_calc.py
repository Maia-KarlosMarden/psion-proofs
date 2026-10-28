import math, numpy as np
Q9 = np.array([1,-5,-5,10,4,4])/9
def r2(t,X):
    X=np.array([X[k] for k in t]); q=Q9[list(t)]; c=X.mean(0)
    v=float(np.sum(q*np.sum((X-c)**2,1))); Q=q.sum(); return v/Q if abs(Q)>1e-12 else v
def rep(name,X):
    rp=r2((0,4,5),X); Rn=r2((0,1,5),X)/rp; RS=r2((0,1,2),X)/rp
    print(f"{name}: R_n={Rn:+.4f} R_S={RS:+.4f} raio-psion/a_hex={np.linalg.norm(X[0][:2]):.3f} escala a(fm)={math.sqrt(0.7071/rp):.3f}")
hexa={k:np.array([math.cos(math.pi*k/3),math.sin(math.pi*k/3),0]) for k in range(6)}
rep("hexagono   ",hexa)
star={k:math.sqrt(3)*np.array([math.cos(math.pi*k/3+math.pi/6),math.sin(math.pi*k/3+math.pi/6),0]) for k in range(6)}
rep("estrela    ",star)
for ph in (45,90,135,180):
    # dobra sobre a aresta externa (apotema sqrt3/2): vertice central gira em torno da aresta
    ap=math.sqrt(3)/2; f=math.radians(ph)
    P={k:np.array([(ap-ap*math.cos(f))*math.cos(math.pi*k/3+math.pi/6),(ap-ap*math.cos(f))*math.sin(math.pi*k/3+math.pi/6),ap*math.sin(f)]) for k in range(6)}
    rep(f"petala {ph:3d}",P)
