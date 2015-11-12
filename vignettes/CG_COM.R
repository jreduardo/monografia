##======================================================================
## Análise de dados de contagem com subdispersão pelo modelo de
## regressão Count-Gamma e COM-Poisson
## 
##                                                        Walmes Zeviani
##======================================================================

##----------------------------------------------------------------------
## Definições da sessão.

options(help_type="html")
require(lattice)
require(latticeExtra)
require(gdata)
require(plyr)
require(ellipse)
require(rootSolve)
require(splancs)

##----------------------------------------------------------------------
## Função de log-verossimilhança para dados com distribuição count
## gamma.

ll.cg <- function(theta, y, X){
    ## theta: vetor de parâmetros/estimativas
    ## y: vetor de dados observados
    ## X: matrix de covariáveis associadas aos valores observados
    eXb <- exp(X%*%theta[-1]) #*theta[1]
    ## Retorna a log verossimilhança para a amostra
    sum(log(pgamma(1, theta[1]*y, theta[1]*eXb)-
                pgamma(1, theta[1]*y+theta[1], theta[1]*eXb)))
}

##----------------------------------------------------------------------
## Função densidade para a distribuição Count-Gamma.

dgammacount <- function(n, alpha, beta){
    pgamma(1, n*alpha, beta)-pgamma(1, n*alpha+alpha, beta)
}

##----------------------------------------------------------------------
## Perfil de log-verossimilhança para o parâmetro alpha do modelo
## count-gamma.

ll.a.cg <- function(theta, alpha, y, X){
    ## theta: vetor de parâmetros/estimativas
    ## alpha: escalar do parâmetro alpha fixo
    ## y: vetor de dados observados
    ## X: matrix de covariáveis associadas aos valores observados
    eXb <- exp(X%*%theta) #*theta[1]
    ## Retorna a log verossimilhança para a amostra com um valor FIXO de
    ## ALPHA
    sum(log(pgamma(1, alpha*y, alpha*eXb)-
                pgamma(1, alpha*y+alpha, alpha*eXb)))
}

ll.a.b0.cg <- function(theta, alpha, b0, y, X){
    ## theta: vetor de parâmetros/estimativas
    ## alpha: escalar do parâmetro alpha fixo
    ## y: vetor de dados observados
    ## X: matrix de covariáveis associadas aos valores observados
    eXb <- exp(X%*%c(b0,theta)) #*theta[1]
    ## Retorna a log verossimilhança para a amostra com um valor FIXO de
    ## ALPHA
    sum(log(pgamma(1, alpha*y, alpha*eXb)-
                pgamma(1, alpha*y+alpha, alpha*eXb)))
}

##----------------------------------------------------------------------
## Função de log-verossimilhança para o modelo COM-Poisson.

ll.com <- function(theta, y, X, slfy=sum(lfactorial(y))){
    ## theta: vetor de parâmetros/estimativas
    ## y: vetor de dados observados
    ## X: matrix de covariáveis associadas aos valores observados
    eXb <- exp(X%*%theta[-1])
    nu <- theta[1]
    ## nu <- 1/log(theta[1])
    ## nu <- 1/theta[1]
    ## nu <- sqrt(theta[1])
    ## nu <- exp(theta[1])
    ## nu <- exp(theta[2])*theta[1]
    ## eXb <- exp(X%*%theta[-1])
    ## t1 <- sum(y*log(eXb))
    ## Xb <- (X%*%theta[-1])^(1/nu)
    ## eXb <- exp(X%*%theta[-1])^(1/nu)
    t1 <- sum(y*log(eXb))
    ## Tem como melhorar via estatística suficiente
    ## t2 <- -nu*sum(log(factorial(y)))
    t2 <- -nu*slfy
    ## Truncando a soma em 30, os dados estão sempre abaixo desse valor,
    ## tem que começar em 0
    s <- 0:30
    ## factorial(s) pode ser passado fora, pode ser feito exp()
    Z <- sapply(eXb, function(lam) sum(lam^s/factorial(s)^nu))
    ## Z <- sapply(Xb,
    ##             function(lam) sum(exp(s*lam)/exp(nu*lfactorial(s))))
    ## factorial(s) pode ser passado fora, pode ser feito exp()
    t3 <- -sum(log(Z))
    ## Retorna a log-verossimilhança para a amostra
    t1+t2+t3
}

## y <- 4:6; X <- matrix(1, nrow=3); theta <- c(1,1)
## sum(dpois(y, lambda=exp(X%*%theta[-1]), log=TRUE))
## ll.com(theta, y, X)

##----------------------------------------------------------------------
## Função de densidade da distribuição COM-Poisson.

dcom <- function(y, beta, nu){
    sapply(y, function(yi) exp(ll.com(c(nu,beta), y=yi, X=1)))
}

##----------------------------------------------------------------------
## Perfil de log-verossimilhança para o parâmetro alpha do modelo
## COM-Poisson.

ll.a.com <- function(theta, alpha, y, X, slfy=sum(lfactorial(y))){
    ## theta: vetor de parâmetros/estimativas
    ## alpha: escalar do parâmetro alpha fixo
    ## y: vetor de dados observados
    ## X: matrix de covariáveis associadas aos valores observados
    eXb <- exp(X%*%theta)
    nu <- alpha
    t1 <- sum(y*log(eXb))
    t2 <- -nu*slfy
    s <- 0:30
    Z <- sapply(eXb, function(lam) sum(lam^s/factorial(s)^nu))
    t3 <- -sum(log(Z))
    t1+t2+t3
}

ll.a.b0.com <- function(theta, alpha, b0, y, X,
                        slfy=sum(lfactorial(y))){
    ## theta: vetor de parâmetros/estimativas
    ## alpha: escalar do parâmetro alpha fixo
    ## y: vetor de dados observados
    ## X: matrix de covariáveis associadas aos valores observados
    eXb <- exp(X%*%c(b0,theta))
    nu <- alpha
    t1 <- sum(y*log(eXb))
    t2 <- -nu*slfy
    s <- 0:30
    Z <- sapply(eXb, function(lam) sum(lam^s/factorial(s)^nu))
    t3 <- -sum(log(Z))
    t1+t2+t3
}

##----------------------------------------------------------------------
## Função que pega o objeto de classe glm Poisson e estima pelo
## Count-Gamma ou COM-Poisson.

poi2other <- function(m0, model){
    ## m0: objeto de classe glm family=poisson
    ## model: função de log-verossmilhança
    ## Extrai a fórmula
    form <- formula(m0)
    form.split <- strsplit(as.character(form), "~")
    ## Extrai a resposta
    yobs <- m0$data[,form.split[[2]]]
    ## Extrai a matriz do modelo
    X <- model.matrix(m0)
    ## Os chutes iniciais
    alpha <- 0.9; gamma <- coef(m0); 
    ## Passa para optim()
    op <- optim(c(alpha, gamma), model, y=yobs, X=X,
                method="L-BFGS-B",
                ## method="SANN",
                hessian=TRUE,
                lower=c(0,rep(-Inf,length(gamma))),
                upper=c(Inf,rep(Inf,length(gamma))),
                control=list(fnscale=-1))
    ## Retorna resultado da optim()
    return(op)
}

##----------------------------------------------------------------------
## Faz o LRT para resultados do Count Gamma, similar a anova.lm.

anova.other <- function(...){
    ## ... sequência de modelos ajustados pela optim()
    model.list <- list(...)
    lls <- sapply(model.list, function(x) x$value)
    nps <- sapply(model.list, function(x) length(x$par))
    cst <- 2*diff(lls)
    pvs <- pchisq(cst, df=diff(nps), lower.tail=FALSE)
    data.frame(ll=lls, npar=nps, two.ll.dif=c(NA,cst),
               npar.dif=c(NA,diff(nps)), pvalue=c(NA,pvs))
}

##----------------------------------------------------------------------
## Retorna tabela de estimativas com erro padrão e p valor a partir da
## saida da optim().

summary.other <- function(m0, conf.int=0.95){
    ## Tabela t para estimativas dos parâmetros
    est <- m0$par
    sdt <- sqrt(diag(-solve(m0$hessian)))
    zval <- est/sdt
    pval <- 2*pnorm(abs(zval), lower=FALSE)
    ci <- est+outer(sdt, c(1,-1)*qnorm((1-conf.int)/2), "*")
    cbind(Estimate=est, Std.Error=sdt, "z value"=zval,
          "P(>|z|)"=pval, lwr=ci[,1], upr=ci[,2])
}

##----------------------------------------------------------------------
## Painel lattice para o gráfico beeswarm.

panel.beeswarm <- function(x, y, subscripts, r, ...){
    xx <- x; yy <- y
    aux <- by(cbind(yy, xx, subscripts), xx,
              function(i){
                  or <- order(i[,1])
                  ys <- i[or,1]
                  yt <- table(ys)
                  dv <- sapply(unlist(yt),
                               function(j) seq(1,j,l=j)-(j+1)/2)
                  if(!is.list(dv)){ dv <- as.list(dv) }
                  xs <- i[or,2]+r*do.call(c, dv)
                  cbind(x=xs, y=ys, subscripts[or])
              })
    aux <- do.call(rbind, aux)
    panel.xyplot(aux[,1], aux[,2], subscripts=aux[,3], ...)
}

##----------------------------------------------------------------------
## Para fazer a representação gráfica de matrizes de correlação.

myplotcorr <- function(corr, outline = TRUE, col = "grey",
                       numbers = FALSE,
                       type = c("full", "lower", "upper"),
                       diag = (type == "full"), bty = "n",
                       axes = FALSE, xlab = "", ylab = "", asp = 1,
                       cex.lab = par("cex.lab"), 
                       cex = 0.75 * par("cex"),
                       mar = 0.1 + c(2, 2, 4, 2),
                       margimnames, ...) ## Argumento que inclui
{
    savepar <- par(pty = "s", mar = mar)
    on.exit(par(savepar))
    if (is.null(corr)) 
        return(invisible())
    if ((!is.matrix(corr)) ||
        (round(min(corr, na.rm = TRUE), 6) < -1) ||
        (round(max(corr, na.rm = TRUE), 6) > 1)) 
        stop("Need a correlation matrix")
    plot.new()
    par(new = TRUE)
    rowdim <- dim(corr)[1]
    coldim <- dim(corr)[2]
    ##   rowlabs <- dimnames(corr)[[1]]
    ##   collabs <- dimnames(corr)[[2]]
    rowlabs <- collabs <- margimnames ## Linha que inclui
    if (is.null(rowlabs)) 
        rowlabs <- 1:rowdim
    if (is.null(collabs)) 
        collabs <- 1:coldim
    ##   rowlabs <- as.character(rowlabs)
    ##   collabs <- as.character(collabs)
    col <- rep(col, length = length(corr))
    dim(col) <- dim(corr)
    type <- match.arg(type)
    cols <- 1:coldim
    rows <- 1:rowdim
    xshift <- 0
    yshift <- 0
    if (!diag) {
        if (type == "upper") {
            cols <- 2:coldim
            rows <- 1:(rowdim - 1)
            xshift <- 1
        }
        else if (type == "lower") {
            cols <- 1:(coldim - 1)
            rows <- 2:rowdim
            yshift <- -1
        }
    }
    maxdim <- max(length(rows), length(cols))
    plt <- par("plt")
    xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", 
                              cex = cex.lab))/(plt[2] - plt[1])
    xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
    ylabwidth <- max(strwidth(collabs[cols], units = "figure", 
                              cex = cex.lab))/(plt[4] - plt[3])
    ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
    plot(c(-xlabwidth - 0.5, maxdim + 0.5),
         c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty,
         axes = axes, xlab = "", 
         ylab = "", asp = asp, cex.lab = cex.lab, ...)
    text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], 
         adj = 1, cex = cex.lab)
    text(cols - xshift, rep(length(rows) + 1, length(cols)), 
         labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
    mtext(xlab, 1, 0)
    mtext(ylab, 2, 0)
    mat <- diag(c(1, 1))
    plotcorrInternal <- function() {
        if (i == j && !diag) 
            return()
        if (!numbers) {
            mat[1, 2] <- corr[i, j]
            mat[2, 1] <- mat[1, 2]
            ell <- ellipse(mat, t = 0.43)
            ell[, 1] <- ell[, 1] + j - xshift
            ell[, 2] <- ell[, 2] + length(rows) + 1 - i - yshift
            polygon(ell, col = col[i, j])
            if (outline) 
                lines(ell)
        }
        else {
            text(j + 0.3 - xshift, length(rows) + 1 - i - yshift, 
                 round(10 * corr[i, j], 0), adj = 1, cex = cex)
        }
    }
    for (i in 1:dim(corr)[1]) {
        for (j in 1:dim(corr)[2]) {
            if (type == "full") {
                plotcorrInternal()
            }
            else if (type == "lower" && (i >= j)) {
                plotcorrInternal()
            }
            else if (type == "upper" && (i <= j)) {
                plotcorrInternal()
            }
        }
    }
    invisible()
}

##======================================================================
## Análise dos dados de número de capulhos.

##----------------------------------------------------------------------
## Descrição dos dados experimentais:
## des: nível desfolha artifical em relação a àrea foliar total da
##   planta.
## est: estádio fenológico da planta.
## rept: repetição, unidade experimental são duas plantas.
## nc: número de capulhos produzidos por 2 plantas.

cap <- expand.grid(rept=1:5,
                   des=seq(0,100,l=5)/100,
                   est=factor(
                       c("vegetativo","botao","flor","maca","capulho"),
                       levels=c("vegetativo","botao","flor","maca",
                                "capulho")))
cap$nc <- c(10,9,8,8,10,11,9,10,10,10,8,8,10,8,9,9,7,7,8,9,8,6,6,5,6,7,
            8,8,9, 10,9,12,7,10,9,8,9,9,10,8,11,10,7,8,8,7,7,7,7,8,10,9,
            8,12,8,7,5, 5,7,5,6,5,7,4,7,8,5,7,6,4,5,5,4,4,5,8,10,7,8,10,
            9,6,6,8,6,9,7,11, 8,9,6,6,6,6,7,3,3,2,4,3,11,7,9,12,11,9,13,
            8,10,10,9,7,7,9,9,8,8, 10,8,10,9,8,10,8,10)
est.ing <- c("vegetative","cotton bud","cotton flower","boll","harvest")
est.por <- c("vegetativo","botão floral","florescimento","maça",
             "capulho")

##----------------------------------------------------------------------
## Ver.

dev.off()
x11(width=9, height=3)

xyplot(nc~des|est, data=cap, layout=c(5,1), col=1,
       xlim=extendrange(c(0:1),f=0.15), xlab="Defoliation level",
       ylab="Number of bolls in two plants",
       strip=strip.custom(bg="gray90", factor.levels=est.ing), r=0.05,
       panel=panel.beeswarm)

##df("cotton.pdf", width=9, height=3)
xyplot(nc~des|est, data=cap, layout=c(5,1), col=1,
       xlim=extendrange(c(0:1),f=0.15),
       xlab="Nível de desfolha artificial", ylab="Número de capulhos",
       type=c("p","smooth"), col.line="gray50",
       strip=strip.custom(bg="gray90", factor.levels=est.por), r=0.05,
       panel=panel.beeswarm)
##ev.off()

##----------------------------------------------------------------------
##gráficos de média vs variâcia

aux1 <- ddply(cap, .(des, est), summarise, m=mean(nc), v=var(nc))
str(aux1)

dev.off()

summary(lm(v~m, aux1))

##df("mv_pt.pdf", width=5, height=5)
xyplot(v~m, aux1, aspect="iso", col=1, type=c("p","g","r"), lty=1,
       col.line="gray50", jitter.x=TRUE,
       xlab=expression(média~amostral~(bar(x))),
       ylab=expression(variância~amostral~(s^2)), xlim=c(0,11),
       ylim=c(0,11), panel=function(...){
           panel.xyplot(...)
           panel.abline(a=0, b=1, lty=2)
       })
##ev.off()

##======================================================================
## Ajuste dos modelos:
##   Poisson (P),
##   Quasi-Poisson (Q),
##   Count Gamma (G) e
##   COM-Poisson (C).
## 
## Modelos aninhados:
##   nulo (0),
##   ~desfolha (1),
##   ~estágio:desfolha (2),
##   ~estágio:(desfolha,d2) (3).

cpP0 <- glm(nc~1, data=cap, family=poisson)
cpP1 <- glm(nc~des+I(des^2), data=cap, family=poisson)
cpP2 <- glm(nc~est:des+I(des^2), data=cap, family=poisson)
cpP3 <- glm(nc~est:(des+I(des^2)), data=cap, family=poisson)
par(mfrow=c(2,2)); plot(cpP3); layout(1)
anova(cpP0, cpP1, cpP2, cpP3, test="Chisq")
logLik(cpP3)

cpQ0 <- glm(formula(cpP0), data=cap, family=quasipoisson)
cpQ1 <- glm(formula(cpP1), data=cap, family=quasipoisson)
cpQ2 <- glm(formula(cpP2), data=cap, family=quasipoisson)
cpQ3 <- glm(formula(cpP3), data=cap, family=quasipoisson)
anova(cpQ0, cpQ1, cpQ2, cpQ3, test="F")

cpG0 <- poi2other(cpP0, model=ll.cg)
cpG1 <- poi2other(cpP1, model=ll.cg)
cpG2 <- poi2other(cpP2, model=ll.cg)
cpG3 <- poi2other(cpP3, model=ll.cg)
summary.other(cpG3)
anova.other(cpG0, cpG1, cpG2, cpG3)

cpC0 <- poi2other(cpP0, model=ll.com)
cpC1 <- poi2other(cpP1, model=ll.com)
cpC2 <- poi2other(cpP2, model=ll.com)
cpC3 <- poi2other(cpP3, model=ll.com)
summary.other(cpC3)
anova.other(cpC0, cpC1, cpC2, cpC3)

## O CP demora mais que o CG, apresenta problemas de estourar devido a
## soma truncada gerar valores -Inf/Inf, tava s=0:50 e troquei para
## s=0:30 para ter sucesso em verossimilhança CP e CG são equivalentes

cbind("Count-Gamma"=cpG3$par, "COM-Poisson"=cpC3$par)
anova.other(cpG0, cpC0, cpG1, cpC1, cpG2, cpC2, cpG3, cpC3)
## A diferença em verossimilhança é muito pequena, praticamente nula.

##----------------------------------------------------------------------
## Imagem da matriz de covariâncias.

vcov <- list(cg=-solve(cpG3$hessian), com=-solve(cpC3$hessian))
vcor <- lapply(vcov, cov2cor)
vcor[["com"]][1:2,1:2]
lapply(vcov, det)

cp.parnames <- c(expression(alpha), expression(beta[0]),
                 expression(beta[11]), expression(beta[12]),
                 expression(beta[13]), expression(beta[14]),
                 expression(beta[15]), expression(beta[21]),
                 expression(beta[22]), expression(beta[23]),
                 expression(beta[24]), expression(beta[25]))
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")

par(mfrow=c(1,2))
myplotcorr(vcor[["cg"]], margimnames=cp.parnames,
           col=colors[5*vcor[["cg"]]+6], mar=c(0,0,2,0))
title("Count-Gamma")
myplotcorr(vcor[["com"]], margimnames=cp.parnames,
           col=colors[5*vcor[["com"]]+6], mar=c(0,0,2,0))
title("COM-Poisson")
layout(1)

##----------------------------------------------------------------------
## Predição intervalar.

pred <- expand.grid(est=levels(cap$est), des=seq(0,1,l=20))
form <- as.formula(do.call(paste,
                           strsplit(as.character(formula(cpP3)), "nc")))
Xf <- model.matrix(form, pred)

tau <- list(cg=cpG3$par[-1], com=cpC3$par[-1])
delta <- list(cg=cpG3$par[1], com=cpC3$par[1])

pred$cg.eta <- Xf%*%tau[["cg"]]
pred$com.eta <- Xf%*%tau[["com"]]

dcom(3, 1, 1)
dpois(3, exp(1))

pred$cg.mu <- sapply(pred$cg.eta,
                     function(m){
                         sum(pgamma(1, shape=delta[["cg"]]*(1:50),
                                    rate=delta[["cg"]]*exp(m)))
                     })
pred$com.mu <- sapply(pred$com.eta,
                      function(m){
                          sum(0:50*dcom(0:50, m, delta[["com"]]))
                      })
pred ## São coincidentes

xyplot(cg.mu+com.mu~des|est, data=pred, type="l") ## Ok

##----------------------------------------------------------------------
##obtendo os IC

##preditos pelo Poisson
aux <- predict(cpP3, newdata=pred, se.fit=TRUE)
aux <- exp(aux$fit+outer(aux$se.fit,
                         qnorm(0.975)*c(lwr=-1, fit=0, upr=1), "*"))
colnames(aux) <- c("Plwr","Pfit","Pupr")
pred <- cbind(pred, aux)

##preditos pelo Quasi-Poisson
aux <- predict(cpQ3, newdata=pred, se.fit=TRUE)
aux <- exp(aux$fit+outer(aux$se.fit,
                         qnorm(0.975)*c(lwr=-1, fit=0, upr=1), "*"))
colnames(aux) <- c("Qlwr","Qfit","Qupr")
pred <- cbind(pred, aux)
str(pred)

##----------------------------------------------------------------------
## Mmatriz de covariâncias condicionais (expressões em
## Daniel-Multivariada, pg 123, teo 3.6).

round(vcov[["com"]][1,],3)
## Covariâncias entre betas e dispersão não nulas
round(vcov[["cg"]][1,],3)
## Covariâncias entre betas e dispersão quase nulas
lapply(vcov, det)
## Determinante próximo de zero indica não ortogonalidade/colinearidade

## O choleski tá assumindo ortogonalidade... esse é um ponto
## importante!!!  mas as estimativas apontam efeito enquando os IC são
## amplos e permitem passar uma horizontal

## Matriz de covariância marginal dos betas
vcov[["cg"]][-1,-1]
vcov[["cg"]][-1,-1]-vcov[["cg"]][-1,1]%*%
    solve(vcov[["cg"]][1,1])%*%vcov[["cg"]][1,-1]

## Matriz de covariância marginal dos betas
vcov[["com"]][-1,-1]
vcov[["com"]][-1,-1]-vcov[["com"]][-1,1]%*%
    solve(vcov[["com"]][1,1])%*%vcov[["com"]][1,-1]

vcovc <- lapply(vcov,
                function(x){ x[-1,-1]-x[-1,1]%*%
                                 solve(x[1,1])%*%x[1,-1] })
vcovc

round(vcov[["cg"]][-1,-1]-vcovc[["cg"]],3)   ## Ortogonalidade
round(vcov[["com"]][-1,-1]-vcovc[["com"]],3) ## Não ortogonalidade

vcorc <- lapply(vcovc, cov2cor)

##df("corrc.pdf", w=10, h=6)
x11()
par(mfrow=c(1,2))
myplotcorr(vcorc[["cg"]], margimnames=cp.parnames[-1],
           col=colors[5*vcorc[["cg"]]+6], mar=c(0,0,2,0))
title("Count-Gamma")
myplotcorr(vcorc[["com"]], margimnames=cp.parnames[-1],
           col=colors[5*vcorc[["com"]]+6], mar=c(0,0,2,0))
title("COM-Poisson")
layout(1)
##ev.off()

## As matrizes de correlação são praticamente as mesmas
round(vcorc[["cg"]]-vcorc[["com"]], 3)

##----------------------------------------------------------------------
## Agora podemos fazer as bandas de predição.

## Preditos pelo count-gamma.
U.cg <- chol(vcov[["cg"]][-1,-1]) ## Marginal
U.cg <- chol(vcovc[["cg"]])       ## Condicional
se.cg <- sqrt(apply(Xf%*%t(U.cg), 1,
                    function(x) sum(x^2)))

aux <- c(pred$cg.eta)+outer(se.cg,
                            qnorm(0.975)*c(lwr=-1, fit=0, upr=1), "*")
aux <- matrix(sapply(aux,
                     function(m){
                         sum(pgamma(1, shape=delta[["cg"]]*(1:50),
                                    rate=delta[["cg"]]*exp(m)))
                     }),
              ncol=3, dimnames=list(NULL, c("Glwr","Gfit","Gupr")))
pred <- cbind(pred, aux)

## Preditos pelo COM-Poisson.
U.com <- chol(vcov[["com"]][-1,-1]) ## Marginal
U.com <- chol(vcovc[["com"]])       ## Condicional
se.com <- sqrt(apply(Xf%*%t(U.com), 1, function(x) sum(x^2)))

aux <- c(pred$com.eta)+outer(se.com,
                             qnorm(0.975)*c(lwr=-1, fit=0, upr=1), "*")
aux <- matrix(sapply(aux,
                     function(m){
                         sum(0:50*dcom(0:50, m, delta[["com"]]))
                     }),
              ncol=3, dimnames=list(NULL, c("Clwr","Cfit","Cupr")))
pred <- cbind(pred, aux)
head(pred)

##----------------------------------------------------------------------
## Gráfico

dev.off()
x11(width=9, height=3)
xyplot(nc~des|est, data=cap, layout=c(5,1), col=1,
       xlim=extendrange(c(0:1),f=0.15),
       xlab="Nível de desfolha artificial",
       ylab="Número de capulhos produzidos",
       strip=strip.custom(bg="gray90", factor.levels=est.por),
       key=list(columns=2,
                text=list(c("Poisson","Quasi-Poisson","Gamma-count",
                            "COM-Poisson")),
                lines=list(lty=1:4, lwd=2, col=c(1,1,2,3))), r=0.05,
       panel=panel.beeswarm)+
    as.layer(xyplot(Pfit+Plwr+Pupr~des|est, data=pred, type="l", col=1,
                    lty=1, lwd=c(2,1,1)))+
    as.layer(xyplot(Qfit+Qlwr+Qupr~des|est, data=pred, type="l", col=1,
                    lty=2, lwd=c(2,1,1)))+
    as.layer(xyplot(Gfit+Glwr+Gupr~des|est, data=pred, type="l", col=2,
                    lty=3, lwd=c(2,1,1)))+
    as.layer(xyplot(Cfit+Clwr+Cupr~des|est, data=pred, type="l", col=3,
                    lty=4, lwd=c(2,1,1)))

##----------------------------------------------------------------------
## Predizendo probabilidade para o nível zero de desfolha.

dgammacount(3, 1, 1)
dpois(3, 1)
dcom(3, log(1), 1)

x <- 0:20
px.po <- dpois(x, lambda=exp(coef(cpP3)[1]))
plot(px.po~x, type="h")

px.cg <- dgammacount(x, delta[["cg"]], delta[["cg"]]*exp(cpG3$par[2]))
plot(px.cg~x, type="h")

px.com <- dcom(x, cpC3$par[2], delta[["com"]])
plot(px.com~x, type="h", col=2)

plot(px.po~x, type="h", ylim=c(0, max(c(px.cg, px.com))), lwd=2)
points(px.cg~c(x-0.05), type="h", col=2, lwd=2)
points(px.com~c(x+0.05), type="h", col=3, lwd=2)

## Mesma distribuição de probabilidades condicional.

##----------------------------------------------------------------------
## Perfil para dispersão.

cg.list <- list(alpha=seq(3.5,7.2,l=30))
cg.list$perfil <- sapply(cg.list$alpha,
                         function(a){
                             op <- optim(cpG3$par[-1], ll.a.cg, alpha=a,
                                         y=cap$nc, X=model.matrix(cpP3),
                                         method="BFGS", hessian=TRUE,
                                         control=list(fnscale=-1))
                             c(op$value, op$par[1])
                         })

with(cg.list, plot(perfil[1,]~alpha))

cg.list$coef <- c(ll=cpG3$value, alpha=cpG3$par[1],
                  sdalpha=summary.other(cpG3)[1,2])
cg.list$dev.perf <- with(cg.list,
                         2*(coef["ll"]-perfil[1,]))
cg.list$dev.quad <- with(cg.list,
                         ((coef["alpha"]-alpha)/coef["sdalpha"])^2)

qchi <- qchisq(0.95, df=1)
cg.list$fperf <- approxfun(cg.list$alpha, cg.list$dev.perf-qchi)
cg.list$limll <- uniroot.all(cg.list$fperf, c(0, 10))
diff(cg.list$limll)
cg.list$limas <- cg.list$coef["alpha"]+
    c(-1,1)*1.96*cg.list$coef["sdalpha"]
diff(cg.list$limas)

dev.off()
x11(width=5, height=5)
keytext <- c("perfil de verossimilhança","aproximação quadrática")
xyplot(cg.list$dev.perf~cg.list$alpha, type="l", lty=1, col=1,
       xlab=expression(alpha),
       ylab=expression(2*(L(hat(alpha))-L(alpha))),
       key=list(columns=1, title="Gamma-Count Model",
                lines=list(lty=1:2, col=1),
                text=list(keytext)))+
    as.layer(xyplot(cg.list$dev.quad~cg.list$alpha, type="l", lty=2,
                    col=1))
trellis.focus("panel", 1,1, highlight=FALSE)
panel.abline(h=qchi, lty=3)
panel.arrows(cg.list$limll, qchi, cg.list$limll, 0, lty=1, length=0.1)
panel.arrows(cg.list$limas, ((cg.list$coef["alpha"]-cg.list$limas)/
                             cg.list$coef["sdalpha"])^2, cg.list$limas,
             0, lty=2, length=0.1)
panel.abline(v=cg.list$coef["alpha"])
trellis.unfocus()

com.list <- list(alpha=seq(3.5,6.5,l=30))
com.list$perfil <- sapply(com.list$alpha,
                          function(a){
                              op <- optim(cpC3$par[-1], ll.a.com,
                                          alpha=a, y=cap$nc,
                                          X=model.matrix(cpP3),
                                          method="BFGS", hessian=TRUE,
                                          control=list(fnscale=-1))
                              c(op$value, op$par[1])
                          })

with(com.list, plot(perfil[1,]~alpha))

com.list$coef <- c(ll=cpC3$value, alpha=cpC3$par[1],
                   sdalpha=summary.other(cpC3)[1,2])
com.list$dev.perf <- with(com.list, 2*(coef["ll"]-perfil[1,]))
com.list$dev.quad <- with(com.list,
                          ((coef["alpha"]-alpha)/coef["sdalpha"])^2)

com.list$fperf <- approxfun(com.list$alpha, com.list$dev.perf-qchi)
com.list$limll <- uniroot.all(com.list$fperf, c(0, 10))
diff(com.list$limll)
com.list$limas <- com.list$coef["alpha"]+
    c(-1,1)*1.96*com.list$coef["sdalpha"]
diff(com.list$limas)

dev.off()
x11(width=5, height=5)
keytext <- c("perfil de verossimilhança","aproximação quadrática")
xyplot(com.list$dev.perf~com.list$alpha, type="l", lty=1, col=1,
       xlab=expression(alpha),
       ylab=expression(2*(L(hat(alpha))-L(alpha))),
       key=list(columns=1, title="COM-Poisson Model",
                lines=list(lty=1:2, col=1),
                text=list(keytext)))+
    as.layer(xyplot(com.list$dev.quad~com.list$alpha, type="l", lty=2,
                    col=1))
trellis.focus("panel", 1,1, highlight=FALSE)
panel.abline(h=qchi, lty=3)
panel.arrows(com.list$limll, qchi, com.list$limll, 0, lty=1, length=0.1)
panel.arrows(com.list$limas, ((com.list$coef["alpha"]-com.list$limas)/
                              com.list$coef["sdalpha"])^2,
             com.list$limas, 0, lty=2, length=0.1)
panel.abline(v=com.list$coef["alpha"])
trellis.unfocus()

## Gráficos lado a lado
models <- c("Gamma-Count Model","COM-Poisson Model")
perfis <- data.frame(dev.perf=c(cg.list$dev.perf, com.list$dev.perf),
                     dev.quad=c(cg.list$dev.quad, com.list$dev.quad),
                     alpha=c(cg.list$alpha, com.list$alpha),
                     model=rep(models, c(length(cg.list$alpha),
                                         length(com.list$alpha))))
str(perfis)
perfis$model <- factor(perfis$model, levels=c(models))

dev.off()
x11(width=9, height=5)
xyplot(dev.perf+dev.quad~alpha|model, data=perfis, type="l", lty=1:2,
       col=1, xlab=expression(alpha),
       ylab=expression(2*(L(hat(alpha))-L(alpha))),
       scales=list(x=list(relation="free")), ylim=c(-0.25,6),
       key=list(columns=2,# title="COM-Poisson Model",
                lines=list(lty=1:2, col=1),
                text=list(keytext)))
trellis.focus("panel", 1, 1, highlight=FALSE)
panel.abline(h=qchi, lty=3)
panel.arrows(cg.list$limll, qchi, cg.list$limll, 0, lty=1, length=0.1)
panel.arrows(cg.list$limas, ((cg.list$coef["alpha"]-cg.list$limas)/
                             cg.list$coef["sdalpha"])^2,
             cg.list$limas, 0, lty=2, length=0.1)
panel.abline(v=cg.list$coef["alpha"])
trellis.unfocus()
trellis.focus("panel", 2, 1, highlight=FALSE)
panel.abline(h=qchi, lty=3)
panel.arrows(com.list$limll, qchi, com.list$limll, 0, lty=1, length=0.1)
panel.arrows(com.list$limas, ((com.list$coef["alpha"]-com.list$limas)/
                              com.list$coef["sdalpha"])^2,
             com.list$limas, 0, lty=2, length=0.1)
panel.abline(v=com.list$coef["alpha"])
trellis.unfocus()

## Ótimo!

##----------------------------------------------------------------------
## Perfil conjunto entre dispersão e intercepto.

niveis <- c(0.9,0.95,0.99)     ## Níveis de confiança das regiões
cortes <- qchisq(niveis, df=2) ## Valores de corte na deviance
keytext <- c("perfil de verossimilhança", "aproximação quadrática",
             expression("perfil para "~alpha~e~beta[0]))

cg.list2 <- list()
grid <- expand.grid(alpha=seq(3,8,l=100), b0=seq(2.13,2.33,l=100))
ele <- ellipse(x=vcov[["cg"]][1:2,1:2], t=5, centre=cpG3$par[1:2])
apply(ele, 2, range)
grid$inout <- inout(grid[,c("alpha","b0")], ele)
levelplot(inout~alpha+b0, grid)
grid <- subset(grid, inout)
dim(grid)
levelplot(inout~alpha+b0, grid, col.regions=3)
cg.list2$grid <- grid

cg.list2$perfil <- apply(cg.list2$grid, 1,
                         function(x){
                             op <- optim(cpG3$par[-(1:2)], ll.a.b0.cg,
                                         alpha=x[1], b0=x[2], y=cap$nc,
                                         X=model.matrix(cpP3),
                                         method="BFGS", #hessian=TRUE,
                                         control=list(fnscale=-1))
                             op$value
                         })

cg.list2$dev <- 2*(cg.list$coef["ll"]-cg.list2$perfil)
cg.list2$el <- lapply(niveis, function(l){
    ellipse(x=vcov[["cg"]][1:2,1:2], level=l)
})

dev.off()
x11(width=5, height=5)
levelplot(cg.list2$dev~cg.list2$grid$alpha+cg.list2$grid$b0,
          xlab=expression(alpha), ylab=expression(beta[0]),
          key=list(columns=1,
                   lines=list(lty=c(1,2,1), col=c(1,1,1), lwd=c(1,1,2)),
                   text=list(keytext)),
          panel=function(..., at, contour, region){
              panel.levelplot(..., at=at)
              panel.contourplot(..., at=cortes, contour=TRUE,
                                region=FALSE)
              panel.lines(cg.list$alpha, cg.list$perfil[2,], col=1,
                          lwd=2)
              panel.abline(v=cg.list$coef["alpha"], h=cpG3$par[2],
                           lty=3)
              lapply(cg.list2$el,
                     function(x){
                         panel.points(x[,1]+cg.list$coef["alpha"],
                                      x[,2]+cpG3$par[2], type="l",
                                      col=1, lty=2)
                     })
          })

com.list2 <- list()
grid <- expand.grid(alpha=seq(2,8,l=100), b0=seq(5,17,l=100))
ele <- ellipse(x=vcov[["com"]][1:2,1:2], t=9, centre=cpC3$par[1:2])
apply(ele, 2, range)
grid$inout <- inout(grid[,c("alpha","b0")], ele)
levelplot(inout~alpha+b0, grid)
grid <- subset(grid, inout)
dim(grid)
levelplot(inout~alpha+b0, grid, col.regions=3)
com.list2$grid <- grid

com.list2$perfil <- apply(com.list2$grid, 1,
                          function(x){
                              op <- optim(cpC3$par[-(1:2)], ll.a.b0.com,
                                          alpha=x[1], b0=x[2], y=cap$nc,
                                          X=model.matrix(cpP3),
                                          slfy=sum(log(
                                              factorial(cap$nc))),
                                          method="BFGS", #hessian=TRUE,
                                          control=list(fnscale=-1))
                              op$value
                          })

com.list2$dev <- 2*(com.list$coef["ll"]-com.list2$perfil)
com.list2$el <- lapply(niveis, function(l){
    ellipse(x=vcov[["com"]][1:2,1:2], level=l)
})

dev.off()
x11(width=5, height=5)
levelplot(com.list2$dev~com.list2$grid$alpha+com.list2$grid$b0,
          xlab=expression(alpha), ylab=expression(beta[0]),
          key=list(columns=1,
                   lines=list(lty=c(1,2,1), col=c(1,1,1), lwd=c(1,1,2)),
                   text=list(keytext)),
          panel=function(..., at, contour, region){
              panel.levelplot(..., at=at)
              panel.contourplot(..., at=cortes, contour=TRUE,
                                region=FALSE)
              panel.lines(com.list$alpha, com.list$perfil[2,], col=1,
                          lwd=2)
              panel.abline(v=com.list$coef["alpha"], h=cpC3$par[2],
                           lty=3)
              lapply(com.list2$el,
                     function(x){
                         panel.points(x[,1]+com.list$coef["alpha"],
                                      x[,2]+cpC3$par[2], type="l",
                                      col=1, lty=2)
                     })
          })

## save.image("coutdata.RData")

##----------------------------------------------------------------------
## Daqui em diante não é mais resultados do artigo.

y <- 4:6; X <- matrix(1, nrow=3); theta <- c(1,1)
sum(dpois(y, lambda=exp(X%*%theta[-1]), log=TRUE))
ll.com(theta, y, X)
ll.cg(theta, y, X)

n <- 500
da <- data.frame(x=runif(n, 3, 6))
da$y <- rpois(n, lambda=da$x)
m0 <- glm(y~x, data=da, family=poisson)
c0 <- poi2other(m0, model=ll.com)
g0 <- poi2other(m0, model=ll.cg)
cov2cor(-solve(c0$hessian))
cov2cor(-solve(g0$hessian))

##----------------------------------------------------------------------
## Estudo com rpanel da forma das funções.

require(rpanel)

com.panel <- function(panel){
    y <- panel$y
    beta <- panel$beta
    nu <- panel$nu
    py <- dcom(y, beta, nu)
    plot(py~y, type="h", ylim=panel$ylim)
    panel
}

panel <- rp.control(y=seq(0:40), ylim=c(0,0.4))
rp.slider(panel, beta, 0, 10, initval=0.5, showvalue=TRUE,
          action=com.panel)
rp.slider(panel, nu, 0.005, 3, initval=1, showvalue=TRUE,
          action=com.panel)

cg.panel <- function(panel){
    y <- panel$y
    beta <- panel$beta
    alpha <- panel$alpha
    py <- dgammacount(y, alpha, beta)
    plot(py~y, type="h", ylim=panel$ylim)
    panel
}

panel <- rp.control(y=seq(0:40), ylim=c(0,0.4))
rp.slider(panel, beta, 0, 10, initval=0.5, showvalue=TRUE,
          action=cg.panel)
rp.slider(panel, alpha, 0.005, 3, initval=1, showvalue=TRUE,
          action=cg.panel)

##----------------------------------------------------------------------
