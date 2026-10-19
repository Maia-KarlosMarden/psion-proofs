# P1 — Passo 5: o numerador 2^(2j) de Z_NP (pré-registro)

Alvo: Z_NP(j) = 2^(2j) / (2j)^2, com n = 2j estados livres. Pergunta: o numerador
2^n e o denominador n^2 saem do mesmo objeto (o traço no espaço de Fock de n modos
fermiônicos), com uma só normalização?

Não há parâmetro livre neste passo; o cálculo é exato.

## Testes

T1. Espaço de Fock de n modos (Jordan–Wigner, n = 1..9). Verificar Tr 1 = 2^n,
    <N> = n/2 e <N(N−1)> = n(n−1)/4 no estado maximamente misto, e a decomposição
    de Wick: Hartree n²/4, Fock −n/4. Identificar (2j)² = 4<N>² (Hartree puro).
T2. Normalizações. Para o quociente numerador/denominador, calcular com
    (a) numerador e denominador ambos normalizados (Tr ρ = 1),
    (b) ambos não normalizados (ρ = 1),
    (c) híbrido: numerador não normalizado, denominador normalizado.
    Ver qual reproduz 2^n/n².
T3. Estado de Gibbs com energia ε uniforme nos n modos: Z(βε)/(2<N>)² em função de
    βε. Mostra para que condição a fórmula (c) vale.
T4. Conferir que a fórmula do paper be(j) = (4j+3)(2j)²/[2^(2j)(2j+1)²] é b0/(Z_NP·N_c).

## Critério

"Numerador derivado" exige que 2^n e n² resultem da MESMA normalização, e que a
condição sobre o estado seja uma consequência das hipóteses da PH, não uma
hipótese nova. Se só o híbrido (c) o reproduz, ou se a condição é nova, o
resultado é "numerador não derivado; condicional a ...", e a condição é
registrada como postulado.
