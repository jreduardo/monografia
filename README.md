<img src = "https://gitlab.c3sl.ufpr.br/eerj12/tccDocument/raw/master/docs/images/tccDocument.png" width=150px align="right" display="block">

# Extensões e Aplicações dos Modelos de Regressão COM-Poisson
-------------------------------------------

 > [Eduardo Elias Ribeiro Junior](https://gitlab.c3sl.ufpr.br/eerj12) &
   [Walmes Marques Zeviani](https://gitlab.c3sl.ufpr.br/walmes)

Este repositório tem por objetivo gerenciar os arquivos textuais
elaborados em meu trabalho de conclusão de curso na Universidade Federal
do Paraná para a obtenção do título de bacharel em Estatística. O
trabalho compreenderá a avaliação e comparação dos modelos de regressão
**COM-Poisson**, nome dado em homenagem ao seus autores Richard Conway e
Willian Maxwell, com aplicações e extensões (inclusão de efeitos
aleatórios e acomodação de inflação de zeros). A elaboração dos
documentos neste repositório é realizada em conjunto com o repositório
[tccPackage], que possui a implementação computacional para ajuste
desses modelos em formato de pacote R.

***

### Introdução ###

Dados de contagem não raramente são analisados utilizando a distribuição
Poisson como distribuição associada a um modelo linear generalizado, no
caso de modelos de regressão. A distribuição Poisson possui somente um
parâmetro, que representa a média e variância. Essa propriedade, chamada
de equidispersão, é particularidade do modelo Poisson e pode não ser
adequada a diversas situações. Quando o modelo Poisson é aplicado sob
negligência desta suposição, o modelo apresenta erros padrões
inconsistentes para as estimativas dos parâmetros e, por consequência,
para toda função desses parâmetros.

A distribuição COM-Poisson é uma distribuição que contempla os casos de
sub e superdispersão devido a adição de mais um parâmetro e pode ser
adapatada em modelos de regressão. Este trabalho apresenta uma revisão
bibliográfica desses modelos, métodos baseados em verossimilhança para
ajuste e comparação de modelos COM-Poisson, extensões para os casos de
excesso de zeros e estrutura com efeitos aleatórios, aplicações dos
modelos COM-Poisson e uma discussão sobre o emprego desses modelos na
análise de dados de contagem.

### Organização do repositório ###

Os arquivos de documentação do projeto de pesquisa e da pesquisa
propriamente dita são mantidos no diretório `docs/` onde estão os
arquivos `00-projeto` e `01-tcc`. Neste diretório todos são elaborados
utilizando o [LaTeX] com o pacote [abnTeX2], para diagramação de texto,
em arquivos de extensão `.Rnw`, que mesclam macros LaTeX com códigos
R. Assim são armazenados somente os arquivos fontes, `.Rnw`. Ainda para o
relatório da pesquisa, em `01-tcc`, decidiu-se trabalhar com partições
do documento sendo essas partições os capítulos do arquivo que são
incluídos no arquivo geral via recursos do [knitr]. 

Nesse repositório de documentação os slides elaborados para
apresentações sobre a pesquisa também são mantido. O diretório
`presentations/` é dedidaco a este fim. Neste diretório os fontes tem
extensão `.Rmd`, que compilados geram um `.pdf`. Os arquivos
`abntcite.csl` e `preambulo-beamer.tex` são utilizados para configuração
do estilo de bibliografia e _layout_/estrutura dos slides respectivamente.


```

.
├── docs
│   ├── 00-projeto.pdf
│   ├── 00-projeto.Rnw
│   ├── 01-tcc.pdf
│   ├── 01-tcc.Rnw
│   ├── cap01_introducao.Rnw
│   ├── cap02_revisao-de-literatura.Rnw
│   ├── cap03_materiais-e-metodos.Rnw
│   ├── cap04_resultados-e-discussao-concordance.tex
│   ├── cap04_resultados-e-discussao.log
│   ├── cap04_resultados-e-discussao.Rnw
│   ├── cap04_resultados-e-discussao.tex
│   ├── cap05_consideracoes-finais.Rnw
│   ├── compois.bib
│   ├── images
│   │   ├── tccDocument.png
│   │   ├── tccDocument.svg
│   │   └── ufpr_bg.pdf
│   └── _setup.R
├── examples
│   ├── moscaBranca.txt
│   ├── ovos.txt
│   ├── Ultimate_Gamma_Count.html
│   └── Ultimate_Gamma_Count.Rmd
├── presentations
│   ├── abntcite.csl
│   ├── images
│   │   ├── logoR.png
│   │   └── regressao.pdf
│   ├── preambulo-beamer.tex
│   ├── slides-projeto.pdf
│   └── slides-projeto.Rmd
├── README.md
└── vignettes
    ├── 1-consRendGraosSoja12-13.csv
    ├── CG_COM.R
    ├── defoliation.pdf
    ├── defoliation.Rnw
    ├── defoliation.tex
    ├── mexilhao.R
    ├── mexilhao.txt
    └── SojaComponentesProd.Rmd

```


Para uma visualização rápida dos documentos resultantes, os arquivos
`.pdf` (resultados de compilação) também são mantidos sob versionamento.

[tccPackage]: https://gitlab.c3sl.ufpr.br/eerj12/tccPackage/
[LaTeX]: https://www.latex-project.org/
[abnTex2]: http://www.abntex.net.br/
[knitr]: http://yihui.name/knitr/
