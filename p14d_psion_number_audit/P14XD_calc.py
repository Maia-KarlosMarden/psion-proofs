import math
hbarc = 1.97327e-16          # GeV*m
MP = 1.2209e19               # GeV (massa de Planck, 1/sqrt(G))
lP = hbarc / MP              # m
be = 325 / 1152; Ne = 10; S = 6 * Ne / be; g2 = 2 * math.pi ** 2 / 3; K = 1.16335
Lam = math.sqrt(K) * MP * math.exp(-S / (2 * g2))
print(f"l_P = {lP:.3e} m; Lambda_meta = {Lam:.4e} GeV")
print("T1 literal: 3 l_P =", f"{3*lP:.2e} m; proton 8.4e-16 m (razao {8.4e-16/(3*lP):.1e}); limite eletron ~1e-20 m (razao {1e-20/(3*lP):.1e})")
fac = math.exp(S / (2 * g2)) / math.sqrt(K)
print(f"T2: e^(S/2g^2)/sqrt(K) = {fac:.4e}; l_P*fator = {lP*fac:.4e} m; hbar c/Lambda = {hbarc/Lam:.4e} m (identidade: {abs(lP*fac/(hbarc/Lam)-1):.1e})")
for l, nome, dado in ((8, "proton (raio de carga)", 8.4e-16), (10, "eletron (limite de tamanho)", 1e-20)):
    m = Lam * math.exp(-l / be); R = hbarc / m
    print(f"T3: R({l}) = 1/m({l}) = {R:.2e} m (m = {m*1e3:.3f} MeV); {nome}: {dado:.1e} m; razao R/dado = {R/dado:.1e}")
print(f"    expoente total Planck->eletron: S/2g^2 + 10/b_e = {S/(2*g2)+10/be:.2f}; ln(M_P/m_e) = {math.log(MP/0.511e-3):.2f}")
A = math.exp(-7 / be)
print(f"T4: sobreposicao e^(-7/b_e) = {A:.1e} contra exigencia 5,5e-30 -> {'passa' if A <= 5.5e-30 else 'FALHA'}")
