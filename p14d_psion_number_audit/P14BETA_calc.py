import math, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
LAM=1.2609e12; T=0.163; L0=math.log(LAM/T)
print("== Parte 1: J.Sigma em 4D ==")
G,mu,Lc,N,lam=sp.symbols('G mu Lambda N lambda',positive=True)
# soma de bolhas: G_eff(mu) = G/(1 - G Pi(mu)), Pi = N lam (Lc^2 - mu^2)/(4 pi^2) (corte 4D, dominante)
Pi=N*lam*(Lc**2-mu**2)/(4*sp.pi**2)
Geff=G/(1-G*Pi); gh=Geff*mu**2*N*lam/(4*sp.pi**2)
beta=sp.simplify(mu*sp.diff(gh,mu))
ghs=sp.symbols('ghat')
# eliminar G: expressar beta em funcao de ghat
Gsol=sp.solve(sp.Eq(gh,ghs),G)[0]
b=sp.simplify(beta.subs(G,Gsol))
print("  beta(ghat) =",sp.factor(sp.simplify(b)), "(em mu << Lambda: 2 ghat - 2 ghat^2 apos reescalar ghat* = 1)")
print("  ponto fixo ghat* = 1 (G_c = 4 pi^2/(N lam Lambda^2)); linearizado: beta ~ -2 (ghat-ghat*) => M ~ Lambda |dG/G|^(1/2), nu = 1/2")
d=(T/LAM)**2*math.log(1+(LAM/T)**2)
print(f"  => para M/Lambda = {T/LAM:.2e} exige-se |G-G_c|/G_c ~ {d:.1e} (como WP-LJS). Sem log: nao ha transmutacao exponencial.")
print("== Parte 2: proto-cor SU(3) ==")
def coeffs(nf): return 11-2*nf/3, 102-38*nf/3
def pole(alpha0,nf,loops):
    b0,b1=coeffs(nf)
    if loops==1:
        e=2*math.pi/(b0*alpha0); return LAM*math.exp(-e), e
    f=lambda t,a: [-(b0/(2*math.pi))*a[0]**2-(b1/(8*math.pi**2))*a[0]**3]
    ev=lambda t,a: a[0]-math.pi; ev.terminal=True
    s=solve_ivp(f,[0,-200],[alpha0],events=ev,rtol=1e-10,atol=1e-14,max_step=0.5)
    t=s.t_events[0][0]; return LAM*math.exp(t), -t
a1=(2*math.pi**2/3)/(4*math.pi)
for nf in (1,0):
    b0,b1=coeffs(nf)
    for loops in (1,2):
        Lp,e=pole(a1,nf,loops); dd=math.log(Lp/T)
        print(f"  C1 alpha=pi/6={a1:.3f}, n_f={nf}, {loops} laco(s): Lambda_pc = {Lp:.3e} GeV (expoente {e:.2f}) ln(prev/alvo) = {dd:+.2f} -> {'ACERTA' if abs(dd)<=0.69 else 'falha'}")
print("  C2 composicao: 1/alpha(mu) = 1/alpha(Lambda) - (b0/2pi) ln(Lambda/mu); com 1/alpha(Lambda) = 0 e b0 > 0,")
print("     1/alpha < 0 logo abaixo de Lambda: inconsistente (gauge induzido exige b0 < 0, triagem). -> nao se aplica")
print("== Exigido (nao e teste) ==")
for nf in (1,0):
    b0,b1=coeffs(nf)
    req1=2*math.pi/(b0*L0)
    req2=brentq(lambda a: math.log(pole(a,nf,2)[0]/T),0.005,0.2)
    print(f"  n_f={nf}: alpha_pc(Lambda_meta) exigido = {req1:.4f} (1 laco), {req2:.4f} (2 lacos)  [1/alpha = {1/req1:.1f}, {1/req2:.1f}]")
ai=1/0.1180 + (23/3)/(2*math.pi)*math.log(172.7/91.19)
ai=ai + 7/(2*math.pi)*math.log(LAM/172.7)
print(f"  referencia MP: alpha_s(Lambda_meta) 1 laco = {1/ai:.4f} [1/alpha = {ai:.1f}]")
