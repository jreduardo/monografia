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
library(xtable)
opts_chunk$set(
    warning = FALSE,
    message = FALSE,
    echo = FALSE,
    tidy = FALSE,
    cache = TRUE,
    results = "hide",
    ## dev = "tikz",
    fig.width = 7,
    fig.height = 5,
    fig.align = "center",
    fig.pos = "h",
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
                  cex = 0.9, outer = TRUE, ...) {
    mtext(texto, cex = cex,
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
        ## ver pág. 187 expresssão 7.10 em Modern Applied, Ripley
        sigma <- sapply(model.list, function(x) summary(x)$dispersion)
        df.sigma <- sapply(model.list, function(x) x$df.residual)
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
## Testa o parametro de dispersão de um modelo quasi poisson
quasitest <- function(...) {
    quasi.list <- list(...)
    cls <- sapply(quasi.list, function(model) model$family$family)
    if (sum(!cls %in% "quasipoisson") > 0) {
        stop("O modelo não é quasipoisson")
    }
    ##-------------------------------------------
    dev <- sapply(quasi.list, function(model) deviance(model))
    dfs <- sapply(quasi.list, function(model) model$df.residual)
    sig <- sapply(quasi.list, function(model) summary(model)$dispersion)
    pvs <- pchisq(q = dev, df = dfs)*2
    out <- cbind("sigma2" = sig, `P(>Chisq)` = pvs)
    rownames(out) <- rep("sigma == 1", length(quasi.list))
    return(out)
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
                   namestrip = NULL,
                   subset = 4,
                   ylab = expression(abs(z)~~(sqrt(~Delta~"deviance"))),
                   xlab = expression(phi),
                   ...) {
    ##-------------------------------------------
    conf <- conf[order(conf, decreasing = TRUE)]
    da <- as.data.frame(prof)
    if (!is.null(subset)) {
        da <- subset(da, abs(z) <= subset)
    }
    ##-------------------------------------------
    fl <- levels(da$param)
    if (!is.null(namestrip)) {
        fl <- namestrip
    }
    xyplot(abs(z) ~ focal | param,
           data = da,
           layout = c(NA, 1),
           xlab = xlab,
           ylab = ylab,
           scales = list(x = "free"),
           type = c("l", "g"),
           strip = strip.custom(
               factor.levels = fl),
           panel = function(x, y, subscripts, ...) {
               conf <- c(0.9, 0.95, 0.99)
               hl <- sqrt(qchisq(conf, 1))
               ##-------------------------------------------
               mle <- x[y == 0]
               xl <- x[x < mle]; yl <- y[x < mle]
               xr <- x[x > mle]; yr <- y[x > mle]
               ##-------------------------------------------
               funleft <- approxfun(x = yl, y = xl)
               funright <- approxfun(x = yr, y = xr)
               vxl <- funleft(hl)
               vxr <- funright(hl)
               vz <- c(hl, hl)
               ##-------------------------------------------
               panel.xyplot(x, y, ...)
               panel.arrows(c(vxl, vxr), 0, c(vxl, vxr), vz,
                            code = 1, length = 0.1, lty = 2,
                            col = "gray40")
               panel.segments(vxl, vz, vxr, vz, lty = 2,
                              col = "gray40")
               panel.abline(h = 0, v = mle, lty = 3)
               panel.text(x = rep(mle, 2), y = vz+0.1,
                          labels = paste(conf*100, "%"),
                          col = "gray20")
           }, ...)
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

##======================================================================
## Tabela de comparação dos coeficientes dos modelos para dados de
## contagem

##-------------------------------------------
## Extrai os coeficientes dos modelos
getCoefs <- function(model, rownames = NULL) {
    cls <- class(model)[1]
    if (!cls %in% c("glm", "negbin", "mle2")) {
        stop("Classe de modelos não contemplada")
    }
    est <- switch(
        cls,
        "glm" = {
            fam <- model$family$family
            if (!fam %in% c("poisson", "quasipoisson")) {
                stop("Família de distribuições não contemplada")
            }
            switch(
                fam,
                "poisson" = {
                    rbind(c(NA, NA), summary(model)$coef[, c(1, 3)])
                },
                "quasipoisson" = {
                    rbind(c(summary(model)$dispersion, NA),
                          summary(model)$coef[, c(1, 3)])
                }
            )
        },
        "negbin" = {
            rbind(c(model$theta, model$SE.theta),
                  summary(model)$coef[, c(1, 3)])
        },
        "mle2" = {
            summary(model)@coef[, c(1, 3)]
        }
    )
    colnames <- c("estimate", "est/std.error")
    if (!is.null(rownames)) {
        rownames(est) <- rownames
    }
    return(est)
}
##-------------------------------------------
## Faz a tabela de comparação
coeftab <- function(..., rownames = NULL) {
    model.list <- list(...)
    coefs <- lapply(model.list, getCoefs, rownames = rownames)
    tab <- do.call(cbind, coefs)
    return(tab)
}
