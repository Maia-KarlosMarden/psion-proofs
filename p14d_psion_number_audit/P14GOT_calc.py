from fractions import Fraction as F
Q={k:F(q,9) for k,q in enumerate((1,-5,-5,10,4,4))}
UNC=[(0,1,2),(0,1,5),(0,2,4),(0,4,5),(1,2,3),(1,3,5),(2,3,4),(3,4,5)]
qT=lambda t: sum(Q[k] for k in t); q2=lambda t: sum(Q[k]**2 for k in t)
for t in UNC: assert (sum(t)%6) in (0,3)
def channels(NQ, Nt):
    out=[]
    for T in UNC:
        for k in range(6):
            for l in {k,(k+3)%6}:
                if qT(T)+Q[k]-Q[l]!=NQ: continue
                if T==Nt and l==k: continue   # canal neutro simetrico (cancela)
                d=q2(T)+Q[k]**2+Q[l]**2-q2(Nt)
                out.append((T,k,l,d))
    return out
p=(0,4,5); n=(0,1,5)
Cp=channels(1,p); Cn=channels(0,n)
print("T1 canais carregados do proton (Q=+1): T', meson (k, l_bar), Delta Q^2")
for T,k,l,d in Cp: print(f"   p -> {T} + ({k},{l}bar): Delta = {d} = {float(d):+.3f}")
print("canais carregados do neutron (Q=0):")
for T,k,l,d in Cn: print(f"   n -> {T} + ({k},{l}bar): Delta = {d} = {float(d):+.3f}")
Sval=q2(p)-q2(n)
mp=max(d for *_,d in Cp); mn=min(d for *_,d in Cn)
Pmax=F(1,2)
Smar_max=Pmax*mp-Pmax*mn
print(f"T2 S_val = {Sval} = {float(Sval):+.4f}; melhor canal p: {float(mp):+.3f}, pior canal n: {float(mn):+.3f}")
print(f"   S_mar,max (P<=0,5 em cada) = {float(Smar_max):+.3f}; S_G,max = {float(Sval+Smar_max):+.3f}")
need=F(235,1000)-Sval
print(f"   S_mar necessario = {float(need):+.3f}; P necessario (cada nucleon, canal otimo) = {float(need/(mp-mn)):.3f}")
Pn=float(need/(mp-mn)); print(f"   -> {'passa' if float(Sval+Smar_max)>=0.209 else 'falha'} (criterio: S_G,max >= 0,209 com P <= 0,5)")
print(f"T3 razao mar necessario / mar QCD = {float(need)/0.098:.2f}")
# sensibilidade: sem simetria dos canais neutros
print(f"   [fora do registro] se os pares neutros de rotulo 3 (Q^2={float(Q[3]**2):.2f}) nao cancelassem: dS = 2 Q3^2 dP; "
      f"dP necessario = {float(need/(2*Q[3]**2)):.3f} (sem mecanismo que o gere)")
