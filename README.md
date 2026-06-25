# BenaCCA — Repositório Oficial

O **BenaCCA** é um software desktop para dimensionamento e análise de
crossovers passivos Butterworth de segunda ordem.

O programa calcula filtros passa-baixas e passa-altas, seleciona componentes
comerciais, compara as respostas ideal e real em gráficos de Bode e gera um
relatório técnico em PDF.

![Interface do BenaCCA](<Documentação Acadêmica/Imagens/interface_benacca.png>)

## Sobre o projeto

O BenaCCA foi desenvolvido para facilitar o projeto de sistemas de áudio de
duas vias. O filtro passa-baixas direciona as frequências menores ao woofer,
enquanto o filtro passa-altas direciona as frequências maiores ao tweeter.

A ferramenta reúne os cálculos, a seleção dos componentes e a análise das
respostas em uma interface única, sem depender de serviços externos.

## Principais recursos

- cálculo dos valores ideais de indutância e capacitância;
- projeto de LPF e HPF Butterworth de segunda ordem;
- seleção do componente comercial mais próximo;
- comparação entre componentes ideais e comerciais;
- gráficos de Bode de magnitude e fase;
- consulta interativa dos pontos das curvas;
- análise da frequência de corte e da tolerância dos componentes;
- geração de relatório técnico em PDF;
- interface desktop local em Tkinter.

## Compatibilidade

| Sistema | Código-fonte | Aplicativo |
|---|---|---|
| Windows 10/11 | `BenaCCA.py` | Preparado para distribuição futura |

## Estrutura do repositório

```text
BenaCCA/
├── Aplicativo/
│   └── README.md
├── Código Fonte/
│   ├── BenaCCA.py
│   ├── executar_benacca.bat
│   ├── requirements.txt
│   └── README.md
├── Documentação Acadêmica/
│   ├── Imagens/
│   ├── Enunciado - Projeto Final.pdf
│   └── RELATORIO_ACADEMICO.md
├── CHANGELOG.md
├── LICENSE
└── README.md
```

## Componentes

### Código-fonte

A pasta [`Código Fonte`](<Código Fonte/README.md>) contém o programa completo,
as dependências e as instruções para execução com Python.

### Aplicativo

A pasta [`Aplicativo`](Aplicativo/README.md) está reservada para a futura
distribuição empacotada do BenaCCA para Windows. O executável não faz parte
desta versão do repositório.

### Documentação acadêmica

O desenvolvimento original também foi apresentado como trabalho da disciplina
de Circuitos de Corrente Alternada da UTFPR.

O relatório exigido para a atividade foi preservado em
[`RELATORIO_ACADEMICO.md`](<Documentação Acadêmica/RELATORIO_ACADEMICO.md>),
junto do enunciado e dos resultados usados na entrega.

## Executar pelo código-fonte

No PowerShell:

```powershell
cd "Código Fonte"
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

## Tecnologias

- Python 3.10 ou superior;
- Tkinter;
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
