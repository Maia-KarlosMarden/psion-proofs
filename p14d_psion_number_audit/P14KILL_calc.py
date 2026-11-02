import math, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
def R2(t,N=6):
    th=2*math.pi*np.array(t)/N; q=Q9[list(t)]; X=np.stack([np.cos(th),np.sin(th)],1); c=X.mean(0)
    v=(q*((X-c)**2).sum(1)).sum(); Qt=q.sum(); return v if abs(Qt)<1e-12 else v/Qt
rp=R2((0,4,5)); Rn=R2((0,1,5))/rp; Rs=R2((0,1,2))/rp
print(f"K1: R_n pred {Rn:+.3f} (obs -0.164 +-0.003/0.0164): dev {abs(Rn+0.164):.3f}; R_Sigma- pred {Rs:+.3f} (obs 0.862 +-0.297): dev {abs(Rs-0.862):.3f} = {abs(Rs-0.862)/0.297:.2f} sigma")
print("   veredito K1:", "KILL" if abs(Rs-0.862)>0.594 else ("SUPPORT" if abs(Rs-0.862)<=0.297 else "inconclusivo"))
M={"pi+-":139.570,"pi0":134.977,"K+-":493.677,"K0":497.611,"eta":547.862,"rho":775.26,"omega":782.66,"eta'":957.78,"phi":1019.461,"p":938.272,"n":939.565,"Lambda":1115.683,"Sigma+":1189.37,"Sigma0":1192.642,"Sigma-":1197.449,"Xi0":1314.86,"Xi-":1321.71,"Delta":1232.0,"Sigma*":1383.7,"Xi*":1531.8,"Omega-":1672.45,"D0":1864.84,"J/psi":3096.9,"B+-":5279.34,"Upsilon":9460.3}
lg=math.log(34.6)
leff=np.array([10-math.log(m/0.5089)/lg for m in M.values()])
fr=lambda l: np.abs(l-np.round(l))
Nin=int((fr(leff)<=0.15).sum())
rng=np.random.default_rng(81); off=rng.random(100000); Nnull=np.array([(fr(leff+o)<=0.15).sum() for o in off])
p=(Nnull>=Nin).mean()
print("K2: l_eff:",{k:round(v,2) for k,v in zip(M,leff)})
print(f"    N_in = {Nin}/{len(M)}; nulo (offset da grade): media {Nnull.mean():.2f}; p = {p:.3f}")
print("    veredito K2:", "KILL" if p>0.20 else ("SUPPORT" if p<=0.05 else "inconclusivo"))
