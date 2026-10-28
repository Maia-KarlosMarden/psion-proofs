import math, itertools, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
def r2(t,X):
    Y=np.array([X[k] for k in t]); q=Q9[list(t)]; c=Y.mean(0)
    v=float(np.sum(q*np.sum((Y-c)**2,1))); Q=q.sum(); return v/Q if abs(Q)>1e-12 else v
def ev(X):
    rp=r2((0,4,5),X)
    if rp<=0: return None
    return r2((0,1,5),X)/rp, r2((0,1,2),X)/rp
ok=lambda r: r is not None and abs(r[0]/-0.164-1)<0.1 and abs(r[1]-0.86)<0.21
P={"N":[(5,0),(1,2),(3,4)],"E":[(0,5),(1,4),(2,3)]}
def star(pairs,d=0.0,beta=0.0,s=1.0):
    X={}
    for j,(a,b) in enumerate(pairs):
        th=math.radians(120*j); u=np.array([math.cos(th),math.sin(th),0]); w=np.array([-math.sin(th),math.cos(th),0]); z=np.array([0,0,1.0])
        wb=math.cos(beta)*w+math.sin(beta)*z  # torcao em torno do eixo radial
        for k,sg in ((a,-1),(b,1)): X[k]=(d+s*math.sqrt(3)/2)*u+sg*(s/2)*wb
    return X
def band(lbl,xs,fn):
    vals=[(x,ev(fn(x))) for x in xs]; g=[x for x,r in vals if ok(r)]
    Rn=[r[0] for _,r in vals if r]
    print(f"  {lbl}: R_n em [{min(Rn):+.3f},{max(Rn):+.3f}]; passam {len(g)}/{len(xs)}"+(f" em [{min(g):.2f},{max(g):.2f}]" if g else ""))
    return vals
print("V1 afastamento d/s:")
for p in "NE":
    v=band(p,np.linspace(0,3,301),lambda d:star(P[p],d=d))
    print("     amostras:",", ".join(f"d={x:.1f}:{r[0]:+.3f}/{r[1]:.2f}" for x,r in v[::50] if r))
print("V2 helice beta (graus):")
for p in "NE":
    v=band(p,np.linspace(0,90,181),lambda b:star(P[p],beta=math.radians(b)))
    print("     amostras:",", ".join(f"b={x:.0f}:{r[0]:+.3f}/{r[1]:.2f}" for x,r in v[::36] if r))
print("V3 permutacoes no anel equiespacado:")
ring=[np.array([math.cos(math.pi*i/3),math.sin(math.pi*i/3)]) for i in range(6)]
passed=[]; vals=set()
for perm in itertools.permutations(range(6)):
    X={k:ring[perm[k]] for k in range(6)}; r=ev(X)
    if r: vals.add(round(r[0],3))
    if ok(r): passed.append(perm)
print(f"  passam {len(passed)}/720; R_n distintos: {len(vals)}")
# classificar: ordem ciclica dos psions ao redor do anel
cyc=set()
for perm in passed:
    order=[k for _,k in sorted((perm[k],k) for k in range(6))]
    i=order.index(0); o=order[i:]+order[:i]
    if o[1]>o[-1]: o=[o[0]]+o[1:][::-1]
    cyc.add(tuple(o))
print("  ordens ciclicas que passam (a menos de rotacao/reflexao):",sorted(cyc))
print("V4 estrela de 9 pontas {9/3}, 6 de 9 pontas em ordem natural:")
pts=[np.array([math.cos(2*math.pi*i/9),math.sin(2*math.pi*i/9)]) for i in range(9)]
good=[]; tot=0
for occ in itertools.combinations(range(9),6):
    tot+=1; X={k:pts[occ[k]] for k in range(6)}; r=ev(X)
    if ok(r): good.append((occ,round(r[0],3),round(r[1],3)))
print(f"  passam {len(good)}/{tot}:",good[:9])
