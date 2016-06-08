##======================================================================
## Função para usar wrapfigure
##   - tentativas usando knitr hooks falharam, os ambientes latex do
##     knitr (knitrout, kframe) não trabalham bem com wrapfigure
##   - Decidiu-se fazer uma função retornar o código latex (utilizar com
##     resultis = "asis")
wrapfigure <- function(align = "L", width = "0.5") {
    with(opts_current$get(), {
        path <- paste0(fig.path, label, "-1")
        label <- paste0("\\label{", fig.lp, label, "}")
        caption <- paste0("\\caption{", fig.cap, "}")
        wraf <- paste0("\\begin{wrapfigure}{",
                       align, "}{", width,
                       "\\textwidth}")
        graf <- paste0("\\includegraphics[width=\\maxwidth]{",
                       path, "}")
        end <- paste0("\\end{wrapfigure}")
        cat(wraf,
            graf,
            caption,
            label,
            end,
            sep = "\n")
    })
    invisible()
}

##======================================================================
## Opções gerais do knitr
library(knitr)
opts_chunk$set(
    warning = FALSE,
    message = FALSE,
    echo = FALSE,
    tidy = TRUE,
    cache = TRUE,
    fig.width = 7,
    fig.height = 5,
    fig.align = "center",
    dev.args = list(family = "Palatino"))

##======================================================================
## Configura opções de gráficos do knitr
library(lattice)
library(latticeExtra)
mycol <- c(1, "#377EB8", "#E41A1C", "#4DAF4A",
           "#ff00ff", "#FF7F00", "#984EA3", "#FFFF33")

## Trellis graphical style.
ps <- list(
    box.rectangle = list(col = 1, fill = c("gray70")),
    box.umbrella = list(col = 1, lty = 1),
    box.dot = list(pch = "|"),
    dot.symbol = list(col = 1, pch = 19),
    dot.line = list(col = "gray50", lty = 3),
    plot.symbol = list(col = 1),
    plot.line = list(col = 1),
    plot.polygon = list(col = "gray95"),
    superpose.line = list(col = mycol, lty = 1),
    superpose.symbol = list(col = mycol, pch = 1),
    superpose.polygon = list(col = mycol),
    strip.background = list(col = c("gray90", "gray70")))

##======================================================================
## Para incluir fonte (source) nas figuras

## da lattice (use o argumento sub = "texto" nos graficos)
ps$par.sub.text <- list(font = 1, just = "left", cex = 0.9,
                        x = grid::unit(5, "mm"))
trellis.par.set(ps)

## da lattice combinados (1. use ps.sub como par.settings, nos gráficos
## bottom da combinação e 2. use fonte.xy("texto") após o print dos
## gráficos)
library(grid)
library(gridExtra)
ps.sub <- list(layout.heights = list(bottom.padding = 5))
fonte.xy <- function(texto, ...) {
    grid::grid.text(texto, x = 0.01, y = 0.01,
                    default.units = "npc",
                    just = c("left", "bottom"), ...)
    invisible()
}


## da graphics (use fonte("texto") após os gráficos)
fonte <- function(texto, side = 1, line = -1, adj = 0,
                  outer = TRUE, ...) {
    mtext("Fonte: Elaborado pelo autor.", cex = 0.9,
          side = side, line = line, adj = adj, outer = outer, ...)
    invisible()
}

##======================================================================
## Configura opções de output no documento
options(digits = 3, OutDec = ",",
        xtable.caption.placement = "top",
        xtable.include.rownames = FALSE,
        xtable.booktabs = TRUE)
