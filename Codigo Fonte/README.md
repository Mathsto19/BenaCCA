# BenaCCA — Código-fonte

Esta pasta contém a implementação completa do BenaCCA.

## Conteúdo

- `BenaCCA.py`: arquivo principal do software;
- `requirements.txt`: dependências necessárias;
- `executar_benacca.bat`: atalho de inicialização no Windows;
- `README.md`: instruções desta pasta.

## Requisitos

- Windows, Linux ou macOS;
- Python 3.10 ou superior;
- PyQt6 (interface gráfica).

## Instalar as dependências

Abra o PowerShell nesta pasta e execute:

```powershell
python -m pip install -r requirements.txt
```

## Executar

```powershell
python BenaCCA.py
```

Como alternativa, abra:

```text
executar_benacca.bat
```

## Dependências

- `PyQt6`: interface gráfica;
- `numpy`: cálculos numéricos;
- `matplotlib`: gráficos de Bode;
- `reportlab`: geração de relatórios em PDF.

Todas as dependências são instaladas via `pip` (veja `requirements.txt`).

## Organização do arquivo principal

O arquivo `BenaCCA.py` está dividido em:

1. tabelas de componentes comerciais;
2. estruturas de dados;
3. cálculos do crossover;
4. geração dos gráficos;
5. geração do relatório;
6. interface gráfica;
7. inicialização do programa.

## Saídas

O usuário escolhe onde salvar os relatórios PDF. Arquivos gerados localmente
não precisam ser adicionados ao repositório.

## Documentação dos cálculos

As fórmulas, os resultados do caso de `2 kHz / 8 Ω` e a análise crítica estão
em
[`RELATORIO_ACADEMICO.md`](<../Documentação Academica/RELATORIO_ACADEMICO.md>).
