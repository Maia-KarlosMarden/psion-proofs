import math, numpy as np
hc=197.327
K=0.921*1.5*hc   # MeV*fm : <|p|> = K/r
Mp,Mn=938.272,939.565
def Tmin(rs,t): return sum(K/rs[i] for i in t)
sets={"FF3 floor":{0:0.41,1:0.44,4:0.69,5:0.34},"JNT2 R1":{0:0.762,1:0.73,4:1.021,5:0.573},"JNT2 R4":{0:0.732,1:0.732,4:1.055,5:0.588}}
for nm,rs in sets.items():
    for hn,t,M in (("p",(0,4,5),Mp),("n",(0,1,5),Mn)):
        T=Tmin(rs,t); print(f"U1 {nm} {hn}: T_min={T:.0f} MeV, T/M={T/M:.2f}, binding needed B=M-T={M-T:.0f} MeV")
print(f"U2 common r threshold (p): r>={3*K/Mp:.3f} fm ; neutron-radius bound r<~0.17 fm -> T_min(0.17 fm, 3 psions)={3*K/0.17:.0f} MeV = {3*K/0.17/Mp:.1f} M_p")
for nm,rs,a in (("FF3 floor",sets["FF3 floor"],0.94),):
    for k,r in rs.items(): print(f"U3 {nm} psion {k}: sigma={r/math.sqrt(3):.3f} fm, a/sigma={a/(r/math.sqrt(3)):.1f}")
print(f"U3 hbar c/a (a=0.94)={hc/0.94:.0f} MeV ; (a=1.212)={hc/1.212:.0f} MeV")
for r in (0.2,0.1,0.05): print(f"U4 r={r}: T_min(3 psions)={3*K/r:.0f} MeV")
