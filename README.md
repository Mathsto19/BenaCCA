# BenaCCA - Repositório Oficial

O **BenaCCA** é um software desktop para dimensionamento e análise de
crossovers passivos Butterworth de segunda ordem.

O programa calcula filtros passa-baixas e passa-altas, seleciona componentes
comerciais, compara as respostas ideal e real em gráficos de Bode e gera um
relatório técnico em PDF.

![Interface do BenaCCA](<Documentação Academica/Imagens/interface_benacca.png>)

## Sobre o projeto

O BenaCCA foi desenvolvido para facilitar o projeto de sistemas de áudio de
duas vias. O filtro passa-baixas direciona as frequências menores ao woofer,
enquanto o filtro passa-altas direciona as frequências maiores ao tweeter.

A ferramenta reúne os cálculos, a seleção dos componentes e a análise das
respostas em uma interface única, sem depender de serviços externos.

## Objetivos e especificações acadêmicas

O projeto atende ao caso obrigatório do enunciado:

- filtro passa-baixas (LPF) Butterworth de 2ª ordem para o woofer;
- filtro passa-altas (HPF) Butterworth de 2ª ordem para o tweeter;
- frequência de corte de `2 kHz`;
- carga resistiva de `8 ohm`;
- seleção de indutores e capacitores exclusivamente a partir das tabelas
  comerciais fornecidas.

Embora a interface permita testar outros valores, os resultados de referência
do trabalho usam `fc = 2000 Hz` e `R = 8 ohm`.

## Funções de transferência e fórmulas

As topologias implementadas usam o denominador comum:

```text
D(s) = R + Ls + RLCs²
```

Para o filtro passa-baixas:

```text
                 R
H_LPF(s) = -----------------
           R + Ls + RLCs²
```

Para o filtro passa-altas:

```text
               RLCs²
H_HPF(s) = -----------------
           R + Ls + RLCs²
```

Para uma resposta Butterworth de segunda ordem:

```text
wc = 2*pi*fc
L = sqrt(2)R / wc
C = 1 / (sqrt(2)Rwc)
```

## Lógica do programa

O arquivo [`BenaCCA.py`](<Codigo Fonte/BenaCCA.py>) é dividido em seções
comentadas. As bibliotecas ficam no início do arquivo e, em seguida, aparecem:

1. tabelas de componentes comerciais;
2. estruturas de dados;
3. cálculos do crossover;
4. gráficos de Bode;
5. desenhos do circuito e análises extras;
6. geração do relatório PDF;
7. interface gráfica em PyQt6.

O fluxo principal valida os parâmetros, calcula os valores ideais, busca o valor
comercial mais próximo, calcula as respostas complexas dos filtros, localiza os
pontos de `-3,01 dB`, plota as curvas e permite exportar um relatório em PDF.

## Principais recursos

- cálculo dos valores ideais de indutância e capacitância;
- projeto de LPF e HPF Butterworth de segunda ordem;
- seleção do componente comercial mais próximo;
- comparação entre componentes ideais e comerciais;
- gráficos de Bode de magnitude e fase;
- consulta interativa dos pontos das curvas;
- análise da frequência de corte e da tolerância dos componentes;
- geração de relatório técnico em PDF;
- interface desktop local em PyQt6.

## Compatibilidade

| Sistema | Código-fonte | Aplicativo |
|---|---|---|
| Windows 10/11 | `BenaCCA.py` | Preparado para distribuição futura |

## Estrutura do repositório

```text
BenaCCA/
├── Aplicativo/
│   └── README.md
├── Codigo Fonte/
│   ├── BenaCCA.py
│   ├── executar_benacca.bat
│   ├── requirements.txt
│   └── README.md
├── Documentação Academica/
│   ├── Imagens/
│   ├── Enunciado - Projeto Final.pdf
│   └── RELATORIO_ACADEMICO.md
├── CHANGELOG.md
├── LICENSE
└── README.md
```

## Componentes

### Código-fonte

A pasta [`Codigo Fonte`](<Codigo Fonte/README.md>) contém o programa completo,
as dependências e as instruções para execução com Python.

### Aplicativo

A pasta [`Aplicativo`](Aplicativo/README.md) está reservada para a futura
distribuição empacotada do BenaCCA para Windows. O executável não faz parte
desta versão do repositório.

### Documentação acadêmica

O desenvolvimento original também foi apresentado como trabalho da disciplina
de Circuitos de Corrente Alternada da UTFPR.

O relatório exigido para a atividade foi preservado em
[`RELATORIO_ACADEMICO.md`](<Documentação Academica/RELATORIO_ACADEMICO.md>),
junto do enunciado e dos resultados usados na entrega.

## Executar pelo código-fonte

No PowerShell:

```powershell
cd "Codigo Fonte"
python -m pip install -r requirements.txt
python BenaCCA.py
```

Também é possível abrir `executar_benacca.bat` dentro da mesma pasta.

## Uso básico

1. Informe a frequência de corte, a impedância e a tolerância.
2. Clique em **Calcular crossover**.
3. Consulte os componentes e os desvios na aba **Resultados**.
4. Abra **Bode: Magnitude e Fase** para comparar as curvas.
5. Gere o relatório em PDF quando precisar registrar os resultados.

## Resultados do caso obrigatório

Para `2 kHz` e `8 ohm`, o programa calcula:

| Componente | Valor ideal | Valor comercial | Diferença |
|---|---:|---:|---:|
| Indutor | `0,9003 mH` | `0,82 mH` | `8,92%` |
| Capacitor | `7,0337 uF` | `6,8 uF` | `3,32%` |

Com os componentes comerciais escolhidos:

| Filtro | Corte ideal | Corte comercial | Desvio |
|---|---:|---:|---:|
| LPF / woofer | `2000 Hz` | `2193,9 Hz` | `+9,70%` |
| HPF / tweeter | `2000 Hz` | `2070,6 Hz` | `+3,53%` |

O fator de qualidade comercial fica em `Q = 0,7285`.

![Bode comparativo](<Documentação Academica/Imagens/bode_comparativo.png>)

## Análise crítica e conclusão

A principal diferença vem do indutor, já que `0,82 mH` é o valor comercial mais
próximo de `0,9003 mH` dentro da tabela exigida. Essa substituição desloca a
frequência de corte, principalmente no LPF, e pode alterar a região de transição
entre woofer e tweeter.

O impacto audível depende dos alto-falantes reais, da caixa, do ambiente e das
tolerâncias dos componentes. O modelo usado no projeto considera a carga como
resistência constante de `8 ohm`, enquanto um alto-falante real tem impedância
variável com a frequência.

O objetivo do projeto foi atingido: os filtros foram dimensionados, os
componentes comerciais foram escolhidos a partir das tabelas obrigatórias e a
resposta ideal foi comparada com a resposta usando componentes reais.

## Tecnologias

- Python 3.10 ou superior;
- PyQt6;
- NumPy;
- Matplotlib;
- ReportLab.

## Versão

Versão inicial documentada: **1.0.0** — 20 de junho de 2026.

O histórico de alterações está disponível em [CHANGELOG.md](CHANGELOG.md).

## Autor

**Matheus Augusto**

Contato: [matheusaugustooliveira@alunos.utfpr.edu.br](mailto:matheusaugustooliveira@alunos.utfpr.edu.br)

## Licença

Copyright © 2026 Matheus Augusto. Todos os direitos reservados.

O código é disponibilizado para consulta e avaliação. Uso, modificação,
redistribuição ou exploração comercial dependem de autorização prévia do
autor. Consulte o arquivo [LICENSE](LICENSE).
