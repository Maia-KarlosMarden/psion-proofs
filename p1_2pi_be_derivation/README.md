# P1: derivar b_e do núcleo 2PI (Hipótese do Psion)

Resultado: **negativo** nas escolhas pré-registradas. b_e = 325/1152 não sai do
núcleo 2PI; segue como postulado calibrado. Ver o relatório conclusivo e a
especificação (documentos do projeto) para a análise completa.

| Passo | Definição | Script | Saída |
| --- | --- | --- | --- |
| 1, 2 | (especificação, seção 7) | `P_2pi_step1.py`, `P_2pi_step2.py` | impressão em tela |
| 4 (D3) | `P1_STEP3_DEFINITION.md` | `P_2pi_step4.py` | `P_2pi_step4_output.txt` |
| 4b (D3′) | `P1_STEP3_D3PRIME_DEFINITION.md` | `P_2pi_step4b.py` | `P_2pi_step4b_output.txt` |
| 4c (D1′) | `P1_D1PRIME_DEFINITION.md` | `P_2pi_step4c.py` | `P_2pi_step4c_output.txt` |
| 5 | `P1_STEP5_DEFINITION.md` | `P_2pi_step5.py` | `P_2pi_step5_output.txt` |

Os arquivos `P1_*_REGISTRY.txt` guardam os SHA-256 e carimbos UTC registrados
antes e depois de cada execução (D3′, D1′, passo 5). Os testes cegos de
look-elsewhere são `blind_test_look_elsewhere.py` e `blind_test_extended.py`.

Os scripts dos passos 4b e 4c importam `P_2pi_step4.py`; rode-os de dentro desta pasta.
Pré-requisitos: Python 3 e numpy.

Tentativas acumuladas: 12 nos ramos principais (D3, D3′, D1′ × 4 ramos σ, N_D),
mais 2 com âncora k = 0. Nenhum ramo foi adotado.
