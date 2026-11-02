import math, numpy as np
import P14FF3_calc as F
sz={0:0.413,1:0.442,4:0.691,5:0.340}
o=F.evalr(sz); a=o[0]; print(f"a = {a:.4f} fm; R_n = {o[1]:+.4f}; perdas FF (p,n,Rn) = {o[2]:.2f}, {o[3]:.2f}, {o[4]:.2f}")
rp2=0.8414**2
def R_of(t,rm):
    lab,qq,d=F.geo(t); return sum(qi*(a*a*di*di+rm[l]**2) for l,qi,di in zip(lab,qq,d))/qq.sum()
grid=np.linspace(0.340,0.691,25)
# Sigma-
vals=[]
for r2 in grid:
    rm=dict(sz); rm[2]=r2; vals.append(R_of((0,1,2),rm)/rp2)
vals=np.array(vals); Robs,sR=0.862,0.212
print(f"Sigma-: R banda [{vals.min():.3f}, {vals.max():.3f}] (obs {Robs} +- {sR}); desvios [{(vals.min()-Robs)/sR:+.2f}, {(vals.max()-Robs)/sR:+.2f}] sigma")
lo,hi=vals.min(),vals.max()
if lo>=Robs-sR and hi<=Robs+sR: v="SUPPORT"
elif lo>=Robs-2*sR and hi<=Robs+2*sR: v="CONSISTENT (banda inteira dentro de 2 sigma)"
elif hi<Robs-2*sR or lo>Robs+2*sR: v="KILL"
else: v="INCONCLUSIVO"
print("veredito Sigma-:",v)
# Sigma+ bands (prediction, no data)
for t,need in (((1,3,5),(3,)),((2,3,4),(2,3))):
    vs=[]
    for r2 in grid:
        for r3 in grid:
            rm=dict(sz); rm[2]=r2; rm[3]=r3; vs.append(math.sqrt(R_of(t,rm)))
    print(f"Sigma+ {t}: raio rms banda [{min(vs):.2f}, {max(vs):.2f}] fm")
