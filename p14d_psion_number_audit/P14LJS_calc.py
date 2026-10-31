import math
from mpmath import mp, mpf, log
mp.dps=60
LAM=1.2609e12; T=0.163; b0=13/6; be=325/1152; Nc=6; g2=2*math.pi**2/3; S=212.677; C=2.5*3.5
x=mpf(T)/mpf(LAM)
D=x**2*log(1+1/x**2)
print(f"A) M/Lambda = {float(x):.3e}; sintonia exigida Delta = 1 - G_c/G = {float(D):.2e}")
print(f"   (para comparacao: M/Lambda = 0,1 exigiria Delta = {0.01*math.log(1+100):.2e})")
forms={"B1 um laco 4D, b0":8*math.pi**2/(b0*g2),"B2 exp(-S/g^2)":S/g2,"B3 degrau l=8":8/be,"B4 Kondo, C=j(j+1)":8*math.pi**2/(C*g2)}
tgt=math.log(LAM/T)
print(f"B) alvo: expoente {tgt:.2f}")
for k,e in forms.items():
    L=LAM*math.exp(-e); d=math.log(L/T)
    print(f"   {k:22s} expoente {e:7.2f} -> Lambda_JS = {L:.3e} GeV | ln(prev/alvo) = {d:+.2f} -> {'ACERTA' if abs(d)<=0.69 else 'falha'}")
