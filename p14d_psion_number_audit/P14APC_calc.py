import io, contextlib, runpy, math
with contextlib.redirect_stdout(io.StringIO()): B=runpy.run_path('P14BETA_calc.py',run_name='x')
pole,T=B['pole'],0.163
S=212.677; g2=2*math.pi**2/3; LUV=S/(2*g2)
c={}
c["D1 alpha_meta = pi/6"]=math.pi/6
c["D2 alpha_psi = 1-9/pi^2"]=1-9/math.pi**2
c["D3 epsilon = 1/6"]=1/6
c["D4 gauge induzido (n_f=1)"]=1/(LUV/(3*math.pi))
ai=1/(math.pi/6)-(31/3)/(2*math.pi)*LUV
c["D5 alpha_meta em Lambda_UV, corrido"]=(1/ai if ai>0 else None)
c["D6 limiar NJL com gauge pi/(3 C_F)"]=math.pi/4
print(f"ln(Lambda_UV/Lambda_meta) = S/2g^2 = {LUV:.3f}; exigido 1/alpha ~ 48,8-55,2")
for k,a in c.items():
    if a is None: print(f"  {k:38s}: 1/alpha(Lambda_meta) = {ai:+.1f} < 0 -> polo ACIMA de Lambda_meta -> falha"); continue
    out=[]; ok=False
    for nf in (1,0):
        for lp in (1,2):
            try: L,_=pole(a,nf,lp); d=math.log(L/T)
            except Exception: d=float('nan')
            out.append(f"nf{nf}/{lp}l:{d:+.1f}")
            if nf==1 and lp==2: ok=abs(d)<=0.69
    print(f"  {k:38s}: alpha = {a:.4f} (1/alpha = {1/a:6.2f}) ln(prev/alvo) {' '.join(out)} -> {'ACERTA' if ok else 'falha'}")
