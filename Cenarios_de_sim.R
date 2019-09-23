library(magrittr)
library(MASS)

#----------------------------------------------------------------------- Função Cpm

Cpm <- function(usl,lsl,sigma,mu){
  tao <- mean(usl,lsl)
  return((usl-lsl)/(6*sqrt((sigma^2)+(mu-tao)^2)))
}

m=c(3,5,10,15,20)
n=c(3,5,10)
cor=c(0.5,0.8,1)

Tao=10
USL=11
LSL=9
dist=Tao-LSL

mu=c(Tao,Tao+dist*0.5,Tao+dist)
sd=c(dist*0.5,dist,dist*1.5) # Verificar coeficiente de variacao
# Nove cenários de Cpm. tres valores de Cpm, com tres comb de sigma e mu p cada Cpm

Cpm(USL,LSL,sd,mu)

teste_cpm <- expand.grid(mu,sd)

Cpm(usl=USL,lsl=LSL,mu=teste_cpm[,1],sigma=teste_cpm[,2])
(cenarios <- expand.grid(m,n,cor,mu,sd))