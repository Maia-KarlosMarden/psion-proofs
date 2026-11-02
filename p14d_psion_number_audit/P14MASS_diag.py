# DIAGNOSTICO POS-HOC (nao pre-registrado): sinais das diferencas de energia do estado fundamental na tabela de WP-MASS
import numpy as np
d=np.load('P14MASS_table.npz'); X=d['X']; T=d['TAB']; H=["p","n","S0","L","Sp1","Sp2"]; E={h:T[:,5+i] for i,h in enumerate(H)}
ok=~np.isnan(T).any(1); print("pontos validos:",int(ok.sum()))
for nm,v in (("E_L-E_n",E["L"]-E["n"]),("E_S0-E_L",E["S0"]-E["L"]),("E_Sp1-E_p",E["Sp1"]-E["p"]),("E_Sp2-E_p",E["Sp2"]-E["p"]),("E_n-E_p",E["n"]-E["p"]),("E_S0-E_n",E["S0"]-E["n"])):
    v=v[ok]; print(f"{nm}: >0 em {int((v>0).sum())}/{len(v)}; mediana {np.median(v):+.3f}; min {v.min():+.3f}; max {v.max():+.3f}")
# ordem observada: n < Lambda < Sigma0 ; p < Sigma+ ; (n ~ p)
ordem=(E["L"]>E["n"])&(E["S0"]>E["L"])&((E["Sp1"]>E["p"])|(E["Sp2"]>E["p"]))
print("ordem observada (n<L<S0 e p<Sp) em", int((ordem&ok).sum()))
for k in (0,1,2,3,4,5):
    pass
# energias medias por hadron
print("energia media do fundamental:",{h:round(float(np.nanmean(E[h])),3) for h in H})
# dependencia em lambda (acoplamento J.Sigma): corr
for h in ("n","L","S0","p"):
    print(f"corr(E_{h}, lambda)={np.corrcoef(E[h][ok],X[ok,3])[0,1]:+.2f}, corr(E_{h}, kappa)={np.corrcoef(E[h][ok],X[ok,2])[0,1]:+.2f}")
