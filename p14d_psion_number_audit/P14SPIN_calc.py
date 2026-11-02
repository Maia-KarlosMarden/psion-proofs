import numpy as np, itertools
Q9=np.array([1,-5,-5,10,4,4])/9
I2=np.eye(2); sx=np.array([[0,1],[1,0]]); sy=np.array([[0,-1j],[1j,0]]); sz=np.diag([1.,-1.])
def op(m,i):
    ms=[I2,I2,I2]; ms[i]=m; return np.kron(np.kron(ms[0],ms[1]),ms[2])
S=[[op(m/2,i) for m in (sx,sy,sz)] for i in range(3)]
def sdot(i,j): return sum(S[i][c]@S[j][c] for c in range(3))
Stot=[sum(S[i][c] for i in range(3)) for c in range(3)]
def moment(t,q,Jsgn):
    isL=[k<=2 for k in t]
    H=sum(Jsgn*sdot(a,b) for a in range(3) for b in range(a+1,3) if isL[a]!=isL[b])
    if not isinstance(H,np.ndarray) or not np.any(H): return np.nan
    # diagonalise H and Sz together: ground multiplet, max Sz
    w,v=np.linalg.eigh(H+1e-6*Stot[2].real*0)
    g=np.where(w<w[0]+1e-9)[0]; V=v[:,g]
    Sz=V.conj().T@Stot[2]@V; ew,ev=np.linalg.eigh((Sz+Sz.conj().T)/2); k=np.argmax(ew); psi=V@ev[:,k]
    mu=sum(q[a]*(psi.conj()@op(sz,a)@psi).real for a in range(3))
    return mu
def ratios(Q,Jsgn,sp):
    qq=lambda t: Q[list(t)]
    mp=moment((0,4,5),qq((0,4,5)),Jsgn)
    return {"n":moment((0,1,5),qq((0,1,5)),Jsgn)/mp,"L":moment((1,2,3),qq((1,2,3)),Jsgn)/mp,"Sp":moment(sp,qq(sp),Jsgn)/mp,"S0":moment((0,2,4),qq((0,2,4)),Jsgn)/mp}
OBS={"n":-0.685,"L":-0.219,"Sp":0.880}
def ok(r): return all(abs(r[k]-OBS[k])<=0.2 for k in OBS)
def best(Q):
    out=[]
    for J in (+1,-1):
        for sp in ((1,3,5),(2,3,4)):
            r=ratios(Q,J,sp); out.append((sum(abs(r[k]-OBS[k]) for k in OBS),J,sp,r))
    return out
if __name__=="__main__":
    for tot,J,sp,r in best(Q9):
        print(f"J={J:+d} Sigma+={sp}: n {r['n']:+.3f} (obs -0.685), Lambda {r['L']:+.3f} (-0.219), Sigma+ {r['Sp']:+.3f} (0.880), Sigma0 {r['S0']:+.3f}; soma|dev|={tot:.3f}; criterio {'OK' if ok(r) else 'falha'}")
    passed=any(ok(r) for _,_,_,r in best(Q9)); print("criterio no modelo:", "PASSA" if passed else "FALHA")
    cnt=0; n=0
    for perm in itertools.permutations(range(6)):
        Q=Q9[list(perm)]
        try:
            if any(ok(r) for _,_,_,r in best(Q)): cnt+=1
        except Exception: pass
        n+=1
    print(f"nulo: {cnt}/{n} permutacoes de Q9 passam (fracao {cnt/n:.3f})")
