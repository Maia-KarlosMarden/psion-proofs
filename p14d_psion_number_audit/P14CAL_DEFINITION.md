# P14 / WP-CAL: calibração do universo PH (pré-registro)

Pedido do autor: dosar as variáveis para reproduzir o universo observado.
Isto é CALIBRAÇÃO: os alvos entram como dados, não contam como previsão. O valor
científico está nas verificações guardadas (não usadas no ajuste) e no custo da
receita (quantos botões, quão finos).

## Regras adicionadas (fora da eq. (1), declaradas)
R1 Violação de CP explícita e global durante a montagem (não domínios: WP-DOM
   mostrou que domínios chegam tarde, são pequenos e dão média zero). Forma: a de
   WP-DOM (CPT e unitariedade por estado inicial), com uma direção v fixa no
   universo inteiro e tamanho ε.
R2 Topo da cascata (matéria escura) como relíquia térmica simétrica, com nível ℓ_DM
   e ⟨σv⟩ = κ/m² livres (κ ≤ 4π, limite de unitariedade).

## Botões
Montagem: c_Γ ∈ {1e-5, 1e-3, 1}, lei ∈ {L1, L2} (6 combinações); direção v no
espaço de 69 dimensões (3 circulações de transição + 66 de rearranjo, WP-DOM V2′);
ε. Escuro: ℓ_DM, κ ∈ {1, 4π}.

## Alvos (ajustados)
A1 η_B = 6,1e-10 em prótons: Δp = D = η/(7,04 Y0), Y0 = 0,208·12/106,75.
A2 Neutralidade: Δe = Δp.
A3 Sem relíquias exóticas: Δ = 0 para as outras 17 tuplas (inclusive (0,1,5)) e
   para os 6 psions livres (cargas fracionárias de Q9 seriam relíquias proibidas).
A4 N_ψ fecha com antineutrinos: Δ(1,2,3) + Δ(0,2,4) = −2D (o balanço de N_ψ fica
   em ν̄_e e ν̄_μ, cargas zero; é a saída que WP-D0′ apontou).
A5 Ω_DM h² = 0,12 pelo topo (R2), resolvendo ℓ_DM para cada κ.
Método: Δ é linear em ε v (verificado em WP-DOM, desvio 4e-6). Matriz resposta G
(26 saídas × 69 direções) por corrida com cada vetor da base, ε = 1e-6. Para cada
combinação (c_Γ, lei), programa linear: minimizar ε_ef = max(|ε x|, |ε c|) sujeito a
A1–A4. Escolher a combinação de menor ε_ef ("a pitada mínima").

## Verificações guardadas (não ajustadas)
V1 Corrida não linear completa com a receita: A1–A4 reproduzidos a 1%?
V2 Sintonia fina: η com c_Γ × 0,5 e × 2, e com v perturbado por ruído de 1% e 10%.
V3 Degenerescência de neutrinos implícita ξ (excesso de ν̄ = 2 n_B, dividido em dois
   sabores) contra o limite |ξ| ≲ 0,03 e contra μ_ν/T = 0,123 de §13.5.
V4 Positividade das taxas (ε_ef < 1).
V5 Massa do topo calibrada contra a âncora do v85 (nível 0,8, 7,3e10 GeV).

## Antecipação registrada
O ajuste deve ser possível (69 botões contra 25 restrições), com ε_ef ~1e-4–1e-3
na melhor combinação. Espero sintonia fina forte em c_Γ (a assimetria depende de
estar perto da transição de equilíbrio), ξ ~1e-8 (trivialmente permitido, mas
incompatível com 0,123 se §13.5 for assimetria leptônica) e ℓ_DM ≈ 4–5 (~100 TeV),
longe do nível 0,8 do paper.
