# BenaCCA — Código-fonte

Esta pasta contém a implementação completa do BenaCCA.

## Conteúdo

- `BenaCCA.py`: arquivo principal do software;
- `requirements.txt`: dependências necessárias;
- `executar_benacca.bat`: atalho de inicialização no Windows;
- `README.md`: instruções desta pasta.

## Requisitos

- Windows 10 ou Windows 11;
- Python 3.10 ou superior;
- suporte ao Tkinter.

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

- `numpy`: cálculos numéricos;
- `matplotlib`: gráficos de Bode;
- `reportlab`: geração de relatórios em PDF.

O Tkinter faz parte da instalação padrão do Python para Windows.

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
[`RELATORIO_ACADEMICO.md`](<../Documentação Acadêmica/RELATORIO_ACADEMICO.md>).
