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
    ## dev = "tikz",
    fig.width = 7,
    fig.height = 5,
    fig.align = "center",
    fig.pos = "!htb",
    dev.args = list(family = "Palatino"))

##======================================================================
## Configura opções de output no documento
options(digits = 3, OutDec = ",",
        xtable.caption.placement = "top",
        xtable.include.rownames = FALSE,
        xtable.booktabs = TRUE)

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
## Para padronizar as analises de deviance de todos os modelos
##   -
myanova <- function(...) {
    model.list <- list(...)
    test <- "Chisq"
    cls <- sapply(model.list, function(x) class(x)[1])
    if (all(cls == "glm")) {
        fam <- sapply(model.list, function(x) x$family$family)
        if (all(fam == "quasipoisson")) {
            test <- "F"
        }
    }
    nps <- sapply(model.list, function(x) attr(logLik(x), "df"))
    aic <- sapply(model.list, function(x) AIC(x))
    if (test == "Chisq") {
        lls <- sapply(model.list, function(x) logLik(x))
        cst <- 2 * diff(lls)
        pvs <- pchisq(cst, df = diff(nps), lower.tail = FALSE)
        tab <- cbind(
            "np" = nps,
            "ll" = lls,
            "AIC" = aic,
            "2*dif" = c(NA, cst),
            "dif np" = c(NA, diff(nps)),
            "P(>Chisq)" = c(NA, pvs))
    }
    if (test == "F") {
        sigma <- summary(model.list[[which.max(nps)]])$dispersion
        df.sigma <- model.list[[which.max(nps)]]$df.residual
        dev <- sapply(model.list, function(x) deviance(x))
        dfs <- c(NA, -diff(dev))
        dnp <- c(NA, diff(nps))
        F <- c(NA, -diff(dev))/(dnp * sigma)
        pvs <- pf(F, df1 = dnp, df2 = df.sigma, lower.tail = FALSE)
        tab <- cbind(
            "np" = nps,
            "dev" = dev,
            "AIC" = aic,
            "F" = F,
            "dif np" = dnp,
            "P(>F)" = pvs)
    }
    return(tab)
}

##======================================================================
## Definições de plot.profile.mle2 para perfil de log-verossimilhança

cols <- trellis.par.get("superpose.line")$col[1:2]
myplot <- function(perfil) {
    par(mar = c(5, 5, 2, 2) + 0.1)
    plot(perfil, conf = c(0.99, 0.95),
         col.minval = 1, lty.minval = 3,
         col.conf = "gray20", lty.conf = 2,
         col.prof = 1, lty.prof = 1,
         ylab = expression(sqrt(-2*(l(hat(phi)) - l(phi)))),
         xlab = expression(phi), las = 1, main = "")
    abline(v = 0, col = cols[2])
    invisible()
}

myprof <- function(prof, conf = c(0.9, 0.95, 0.99),
                   ylab = expression(abs(z)~~(sqrt(~Delta~"deviance"))),
                   xlab = expression(phi),
                   namestrip = expression("Perfil para"~phi),
                   subset = 3.8,
                   ...) {
    ##-------------------------------------------
    conf <- conf[order(conf, decreasing = TRUE)]
    vx <- sapply(conf, function(x) confint(prof, level = x))
    vz <- sqrt(qchisq(pmax(0, pmin(1, conf)), 1))
    ##-------------------------------------------
    ylab <- ylab
    xlab <- xlab
    ##-------------------------------------------
    mle <- prof@summary@coef["phi", 1]
    ##-------------------------------------------
    da <- as.data.frame(prof)
    xyplot(abs(z) ~ focal | param, data = da,
           subset = abs(z) < subset,
           type = c("l", "g"),
           strip = strip.custom(
               factor.levels = namestrip
           ),
           xlab = xlab, ylab = ylab,
           ..., panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.arrows(c(vx), 0, c(vx), rep(vz, each = 2),
                            code = 1, length = 0.1, lty = 2,
                            col = "gray40")
               panel.segments(vx[1, ], vz, vx[2, ], vz, lty = 2,
                              col = "gray40")
               panel.abline(h = 0, v = mle, lty = 3)
               panel.text(x = rep(mle, 2), y = vz+0.1,
                          labels = paste(conf*100, "%"),
                          col = "gray20")
               panel.abline(v = 0, col = cols[2])
           })
    ##-------------------------------------------
}

##======================================================================
## Gráfico de correlação entre os parâmetros

library(corrplot)
mycorrplot <- function(Corr, ...) {
    corrplot.mixed(
        Corr, lower = "number", upper = "ellipse",
        diag = "l", tl.pos = "lt", tl.col = "black",
        tl.cex = 0.8, col = brewer.pal(9, "Greys")[-(1:3)], ...)
}
