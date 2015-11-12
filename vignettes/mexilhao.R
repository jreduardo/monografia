#==========================================================================================
# Experimento do lactec para desenvolvimento de produtos para controle de mexilhão
# dourado nas tubulações
#==========================================================================================

rm(list=ls()) # limpa área de trabalho

#------------------------------------------------------------------------------------------
# importação e descrição dos dados

da <- read.table("mexilhao.txt", header=TRUE, sep="\t")
str(da)

## Coleta: index da coleta, começa na 5 pois antes disso não oservou-se mexilhão
## Data.de.Coleta: data da coleta
## Produto: 4 níveis de produtos empregados no manejo/eliminação do mexilhão dourado
## Numero.da.placa: index da placa, unidade experimental
## Lado.da.placa: 3 níveis onde são encontrados os mexilhões (superior, inferior, lama)
## Comprimento: resposta observada, é o comprimento de cada mexilhão encontrado
## Largura: resposta observada, é a largura de cada mexilhão encontrado

## informações da conversa com o Leon
## são 12 tubulações. elas estão arrajadas em grupos de 4 formando 3 andares de tubulação.
## cada tubulação de cada paramar recebe um nível de produto (3 blocos, 4 produtos, 12 tubos)
## cada tubo possui uma placa contendo 12 chapas. uma chama é retirada a cada data e avaliada.
## são previstas 12 coletas, em cada usa-se uma chapa. elas são removidas e não substituídas.
## as chamas são retiradas na ordem e não aleatóriamente. os níveis de produto foram atribuídos
## sem aleatorização aos tubos. cada tubo tem 12 placas, cada patamar tem 4*12=48 placas. então
## as placas de 1 à 48 são do primeiro patamar. com isso recuparar os níveis de patamar.
## com o produto dos níveis de produto com patamar recuperar a id dos tubos, de 1 à 12 e também
## de placa dentro do tudo de 1 à 12. o número da coleta é o index da placa observada.
## com os dados do controle é possível fazer a curva de crescimento do mexilhão, embora não
## tenhamos informação sobre dias de vida, apenas sobre dias dentro da tubulação.

# mudando os nomes
names(da) <- c("clt", "data", "prod", "placa", "lado", "comp", "larg")

# convertendo data para data
da$data <- as.Date(as.character(da$data), format="%m/%d/%Y")
str(da)

#==========================================================================================
# análise exploratória dos dados

require(lattice)

xyplot(comp~data, data=da)
xyplot(comp~data, data=da, groups=prod, auto.key=TRUE, jitter.x=TRUE)
xyplot(comp~data|prod, data=da, auto.key=TRUE, jitter.x=TRUE)
xyplot(comp~data|prod, data=da, groups=lado, auto.key=TRUE, jitter.x=TRUE, type=c("p","a"))
xyplot(comp~data, data=subset(da, prod="Controle"), groups=lado, auto.key=TRUE, jitter.x=TRUE, type=c("p","a"))
xyplot(comp~data|lado, data=subset(da, prod="Controle"), auto.key=TRUE, jitter.x=TRUE, type=c("p","a"))

xyplot(comp~larg, data=da)
xyplot(comp~larg, data=da, groups=prod, auto.key=TRUE)
xyplot(comp~larg, data=da, groups=data, auto.key=TRUE)
xyplot(comp~larg|data, data=da, groups=prod, auto.key=TRUE)
xyplot(comp~larg|prod, data=da, auto.key=TRUE)
xyplot(comp~larg|prod, data=da, groups=lado, auto.key=TRUE)

histogram(~larg|prod, data=da)
histogram(~larg|data, data=da)

da$controle <- "tratato"
da$controle[da$prod=="Controle"] <- "não tratado"
da$comp.class <- cut(da$comp, seq(0,20,1))

aux <- with(da, aggregate(comp, list(contr=controle, dt=data, cc=comp.class), length))
barchart(x~cc|dt, groups=contr, data=aux, box.ratio=5,
         scales=list(x=list(rot=90)), auto.key=list(columns=2),
         ylab="Número de mexilhões encontrados", xlab="Classe de comprimento (cm)")

#------------------------------------------------------------------------------------------
# tabelas de frequência

with(da, tapply(comp, list(prod, data), length))
with(subset(da, prod="Controle"), tapply(comp, list(lado, data), length))

#------------------------------------------------------------------------------------------
# gráfico de barras dos totais

aux <- with(da, aggregate(comp, list(prod=prod, data=data, lado=lado), length))
barchart(x~prod|data, groups=lado, data=aux)

#==========================================================================================
# analisar a contagem

#------------------------------------------------------------------------------------------
# transformar os números das placas nos indices de bloco e tubo
# não existe efeito de bloco

da$bl <- cut(da$placa, seq(1, by=12*4, length.out=4), include.lowest=TRUE, right=FALSE) # 3 blocos
da$tb <- cut(da$placa, seq(1, by=12, length.out=13), include.lowest=TRUE, right=FALSE) # 12 tubos
str(da)

#------------------------------------------------------------------------------------------
# o fator lado de fixação foi ignorado

db <- as.data.frame(with(da, table(prod, bl, data)))
names(db) <- tolower(names(db))
db$dt <- as.Date(as.character(db$data))
str(db)

#------------------------------------------------------------------------------------------
# gráfico com as contagens

barchart(freq~prod|data, groups=bl, data=db)
xyplot(freq~dt, groups=prod, data=db, type=c("p","a"))

#------------------------------------------------------------------------------------------
# a seguir virão ajustes com diversas especificações de modelos na tentativa de encontrar
# o que melhor represente o comportamento dos dados. Teste de razão de verossimilhanças
# entre modelos encaixados serão aplicados. Com o modelo final serão aplicadas inferências
# como predições e contrastes. Os contrastes são:
# 1. controle vs resto
# 2. cloro+mdx vs NaOH
# 3. cloro vs mdx

#------------------------------------------------------------------------------------------
# análise com glm(family=poisson)

g0 <- glm(freq~prod*data, data=db, family=quasipoisson)
summary(g0)
anova(g0, test="F")

par(mfrow=c(2,2))
plot(g0)
layout(1)

#------------------------------------------------------------------------------------------
# usando um modelo inflacionado de zeros

require(pscl)   # modelos inflacionados de zero
require(lmtest) # teste de razão de verossimilhanças

g2 <- zeroinfl(freq~prod*data|prod, data=db, dist="poisson")
g3 <- zeroinfl(freq~prod*data|prod, data=db, dist="negbin")
lrtest(g3, g2) # poisson vs negbinom, negbin superior

g4 <- zeroinfl(freq~prod*data|data, data=db, dist="negbin")
g5 <- zeroinfl(freq~prod*data|1, data=db, dist="negbin")
lrtest(g3, g5) # sem efeito de prod no excesso de zeros
lrtest(g4, g5) # sem efeito de data no excesso de zeros

g6 <- zeroinfl(freq~prod+data|1, data=db, dist="negbin")
lrtest(g6, g5) # tem interação

#------------------------------------------------------------------------------------------
# modelo com exesso de zeros final

summary(g5)

#------------------------------------------------------------------------------------------
# comparar os modelos com excesso de zero com os sem excesso de zero

g0 <- glm(freq~prod*data, data=db, family=poisson)
g7 <- zeroinfl(freq~prod*data|1, data=db, dist="poisson")
lrtest(g7, g0) # não precisa ser um modelo inflacionado poisson

require(MASS) # função glm.nb()
g8 <- glm.nb(freq~prod*data, data=db)
lrtest(g5, g8) # não precisa ser um modelo inflacionado negbin, ver o porque pois excesso de 0 é nitido

#------------------------------------------------------------------------------------------
# muda o nível de referência

db$prod <- relevel(db$prod, ref="Controle") # muda o nível de referência para Controle
levels(db$prod)

#------------------------------------------------------------------------------------------
# entrando com a data de forma númerica, a coleta mais antiga é a origem

db$dtn <- as.numeric(db$dt)-min(as.numeric(db$dt))
g9 <- glm.nb(freq~prod*(dtn+I(dtn^2)), data=db)
lrtest(g8, g9) # modelo com tempo numérico é suficiente

#------------------------------------------------------------------------------------------
# modelo final

summary(g9)

#------------------------------------------------------------------------------------------
# fazendo gráfico com valores observados e preditos

range(db$dtn)
new <- expand.grid(prod=levels(db$prod), dtn=seq(0,220,1))
new$freq <- predict(g9, newdata=new, type="response")
new$tipo <- "predito"
db$tipo <- "observado"
str(new)

#------------------------------------------------------------------------------------------
# gráfico

dc <- rbind(db[,c("prod","dtn","freq","tipo")], new)
str(dc)

xyplot(freq~dtn|prod, groups=tipo, data=dc,
       distribute.type=TRUE, type=c("p","l"), ylim=extendrange(db$freq))

#------------------------------------------------------------------------------------------
# tratamento final
# *para o nível NaOH remover o termo quadrático
# *fazer o gráfico com as bandas de 95% de confiança para E(y|...)
# *não há métodos da contrast::contrast(), multcomp::glht(), gmodels::fit.contrast() para
# classe "negbin", então o trabalho será matricial

#------------------------------------------------------------------------------------------
# no procedimento a seguir será removido o termo quadrático para NaOH
# prodecimento é criar uma anti-indicadora para NaOH

db$ai <- 1
db$ai[db$prod=="NaOH"] <- 0

g10 <- glm.nb(freq~prod*(dtn+I(ai*dtn^2)), data=db)
lrtest(g9, g10) # remoção do termos quadrádico para NaOH

summary(g10) # presença de NA para o termo quadrático do NaOH

#------------------------------------------------------------------------------------------
# fazendo a prediçãom com o modelo g10

new <- expand.grid(prod=levels(db$prod), dtn=seq(0,220,1))
new$ai <- 1; new$ai[new$prod=="NaOH"] <- 0
new$freq <- predict(g10, newdata=new, type="response")
new$tipo <- "predito"
db$tipo <- "observado"
str(new)

#------------------------------------------------------------------------------------------
# gráfico da predição para o modelo g10

dc <- rbind(db[,c("prod","dtn","ai","freq","tipo")], new)
str(dc)

xyplot(freq~dtn|prod, groups=tipo, data=dc,
       distribute.type=TRUE, type=c("p","l"), ylim=extendrange(db$freq))

#------------------------------------------------------------------------------------------
# aprimorar o gráfico com o IC, o erro padrão do ajuste é na escala do preditor linear
# link é log()

aux <- predict(g10, newdata=new, se=TRUE)
new$eta <- aux$fit
new$se <- aux$se.fit
new$lwr <- with(new, exp(eta-1.96*se))
new$upr <- with(new, exp(eta+1.96*se))
new$freq <- predict(g10, newdata=new, type="response")
str(new)

#------------------------------------------------------------------------------------------
# junta os valores observados com predito e IC

require(reshape) # função melt()
new2 <- melt(new[,c("prod","dtn","freq","lwr","upr")], id=c("prod","dtn"))
str(new2)
names(new2)[3:4] <- c("tipo","freq")

dc <- rbind(db[,c("prod","dtn","tipo","freq")], new2)
str(dc)
unique(dc$tipo)

#------------------------------------------------------------------------------------------
# gráfico final

xyplot(freq~dtn|prod, groups=tipo, data=dc, #ylim=extendrange(db$freq),
       distribute.type=TRUE, # ordem alfabética: freq, lwr, observado, upr
       type=c("l","l","p","l"), col=c(2,3,1,3), lty=c(1,2,0,2),
       scales="free", # eixos livres
       jitter.x=TRUE, strip=strip.custom(bg="gray90"),
       key=list(space="top", divide=1, type="o", columns=3,
         lines=list(lty=c(0,1,2), pch=c(1,NA,NA), col=c(1,2,3)),
         text=list(c("Observado","Predito","IC 95%"))))

#------------------------------------------------------------------------------------------
# conduzir os contrastes de interesse em cada uma das 8 datas
# para evitar ter que manipular os coeficientes sob a parametrização usada, etc, vou
# criar um modelo de classe lm com mesma estrutura e usar o vetor do contrastes retornado
# pela contrast::contrast()$X, daí sigo com as operações matriciais "na mão"
# detalhe: não pode ter estimativa NA, então usar termos quadrático para NaOH
# detalhe: não pode usar o operador I(), criar uma nova variável

formula(g10)                            # formula do modelo final
db$dtn2 <- with(db, dtn^2)              # cria a nova variável dtn2 = dtn^2
m0 <- lm(freq~prod*(dtn+dtn2), data=db) # modelo de mesmo preditor
cbind(lm=coef(m0), negbin=coef(g10))    # verifica estruturas

#------------------------------------------------------------------------------------------
# aplicar contrastes pela função e comparar pelo método matricial

require(contrast) # função constrast

# contraste: Controle vs Outros na data 0
c1 <- contrast(m0, type="average",
               list(prod=levels(db$prod)[1], dtn=0, dtn2=0),
               list(prod=levels(db$prod)[-1], dtn=0, dtn2=0))
c1$X # esse procedimento é para obter esse vetor de contraste que será usado adiante

dim(vcov(g10)) # não pode ter NA pois dá diferença de dimensções entre coef() e vcov()

#------------------------------------------------------------------------------------------
# o contraste matricialmente é calculado assim, conferir com os resultados da contrast

X <- matrix(c1$X, nrow=1)     # vetor do contraste
b <- matrix(coef(m0), ncol=1) # vetor de estimativas
V <- vcov(m0)                 # matriz de covariância das estimativas

X%*%b                                        # estimativa do contraste
sqrt(X%*%V%*%t(X))                           # erro padrão dessa estimativa
Fc <- t(X%*%b)%*%solve(X%*%V%*%t(X))%*%X%*%b # valor F do contraste
pf(Fc, df1=nrow(X), df2=df.residual(g10), lower.tail=FALSE)

c1 # os mesmos resultados do procedimento matricial

#------------------------------------------------------------------------------------------
# repetir esse procedimento para as 8 datas e 3 contrastes, ao todo 8*3=24 comparações
# reunião dos ingredientes para fazer os contrastes que são para o modelo negbin
# os contrastes são na escala do preditor, não interpretar o valor em si, só a direção

dtin <- unique(db$dtn); dtin # datas onde serão realidados os contrastes

p <- length(c1$X)
b <- matrix(coef(g10)[-p], ncol=1) # vetor de estimativas
V <- vcov(g10)                     # matriz de covariância das estimativas
df <- df.residual(g10)             # dim(y)-dim(b)

do.contrast <- function(X, b, V, df){ # função para realizar contrastes
  ## X: *vetor* de contrastes
  ## b: vetor de estivativas
  ## V: matriz de covariâcia das estimativas
  ## df: dimensão(y)-dimensão(b)
  est <- X%*%b                               # valor do contraste
  se <- sqrt(X%*%V%*%X)                      # erro padrão do contraste
  Fc <- t(X%*%b)%*%solve(X%*%V%*%X)%*%X%*%b  # valor F do contraste
  pc <- pf(Fc, df1=1, df2=df, lower.tail=FALSE)
  ## retorna um vetor
  return(c(Contrast=est, S.E.=se, t=sqrt(Fc), "Pr(>|t|)"=pc))
}

#------------------------------------------------------------------------------------------
# contraste 1. controle vs outros em 8 datas de interesse

c1X <- sapply(dtin,
              function(x){
                c0 <- contrast(m0, type="average",
                               list(prod=levels(db$prod)[1], dtn=x, dtn2=x^2),
                               list(prod=levels(db$prod)[-1], dtn=x, dtn2=x^2))
                c0$X
              })
c1X <- t(c1X)[,-p]
rownames(c1X) <- dtin

t(apply(c1X, 1, do.contrast, b=b, V=V, df=df))

#------------------------------------------------------------------------------------------
# contraste 2. cloro+mdx vs NaOH em 8 datas de interesse

c2X <- sapply(dtin,
              function(x){
                c0 <- contrast(m0, type="average",
                               list(prod=levels(db$prod)[2:3], dtn=x, dtn2=x^2),
                               list(prod=levels(db$prod)[4], dtn=x, dtn2=x^2))
                c0$X
              })
c2X <- t(c2X)[,-p]
rownames(c2X) <- dtin

t(apply(c2X, 1, do.contrast, b=b, V=V, df=df))

#------------------------------------------------------------------------------------------
# contraste 3. cloro vs mdx em 8 datas de interesse

c3X <- sapply(dtin,
              function(x){
                c0 <- contrast(m0, type="average",
                               list(prod=levels(db$prod)[2], dtn=x, dtn2=x^2),
                               list(prod=levels(db$prod)[3], dtn=x, dtn2=x^2))
                c0$X
              })
c3X <- t(c3X)[,-p]
rownames(c3X) <- dtin

t(apply(c3X, 1, do.contrast, b=b, V=V, df=df))

#------------------------------------------------------------------------------------------
# fim
#------------------------------------------------------------------------------------------
