
# bootstrap: Modarris (2007) >> reamostragem para amostragem por conjuntos ordenados
#                            >> n�o faz men��o a ordena��o perfeita

# Artigo escrito com o Angelo: Cpm com estima��o pontual >> replicar um cen�rio de simula��o utilizado no artigo.


#----------------------------------------------------------------------- Gerando as amostras. A AAS vai ser a ACO com corr = 0.

# tira o output da AAS. A AAS � o ACO com correla��o igual a zero.

# para intervalo bootstrap, temos que armazenar os valores da covari�vel tamb�m...
ACOeAAS=function(cor,sdx,sdy,tamanho,ciclos,mu_x,mu_y){
  if(cor==0){
    AAS <- rnorm(tamanho*ciclos,mean=mu_y,sd=sdy)
    return(AAS)
  }else{
    sigma <- matrix(c(sdx**2,cor*sdx*sdy,
                      cor*sdx*sdy,sdy**2),byrow=FALSE,ncol=2)
    obj <- MASS::mvrnorm(n=(tamanho**2)*ciclos,mu=c(mu_x,mu_y),
                         Sigma=sigma)
    library(magrittr)
    ACO <- split.data.frame(obj,
                            seq(1,ciclos*tamanho^2,
                                by=tamanho)) %>%
      lapply(function(x){x[order(x[,1]),2]})  %>%
      unlist() %>%
      split(rep(c(1:ciclos),each=tamanho*tamanho)) %>% 
      lapply(function(x){matrix(x,ncol=tamanho,byrow = FALSE)}) %>%
      lapply(diag) %>%
      unlist()
    return(ACO)
  }
}
# tirar a AAS a partir da ACO. Fazer com que a f� retorne dois outputs

#----------------------------------------------------------------------- Fun��o Cpm

Cpm <- function(usl,lsl,sigma,mu){
  tao <- mean(c(usl,lsl))
  return((usl-lsl)/(6*sqrt((sigma^2)+(mu-tao)^2)))
}

#----------------------------------------------------------------------- Gerando os IC's assint�ticos (come�a na p�gina 78)
# o "delta" foi definido que ser� o n�vel de signific�ncia

# Marcucci and Beazley (1988)
ic.mb <- function(cpm, delta, n){
  cpmL <- cpm*sqrt(qchisq(delta/2,df=n)/n)
  cpmU <- cpm*sqrt(qchisq(1-delta/2,df=n)/n)
  return(c("Inf" = cpmL,"Sup" = cpmU))
}

# Chan et.al (1990) 
ic.cxz <- function(cpm, delta, usl, lsl, s, x_bar, t, n){
  d <- (usl-lsl)/2
  sigma_m_chapeu <- ((d/3)^2)*((s^2)*((x_bar-t)^2)+((s^4)/2))/((s^2)+(x_bar-t)^2)^3
  cpmL <- cpm - qnorm(1-delta/2)* sqrt(sigma_m_chapeu/n)
  cpmU <- cpm + qnorm(1-delta/2)* sqrt(sigma_m_chapeu/n)
  return(c("Inf" = cpmL,"Sup" = cpmU))
} 

# Zimmer and Hubele (1997)
ic.zh <- function(cpm,n,s,delta,x_bar,t){
  lambda <- (n*(x_bar-t)^2)/(s^2)
  cpmL <- cpm*sqrt(qchisq(delta/2,n,ncp=lambda)/(n+lambda))
  cpmU <- cpm*sqrt(qchisq(1-delta/2,n,ncp=lambda)/(n+lambda))
  return(c("Inf" = cpmL,"Sup" = cpmU))
}

# Boyles (1991)
ic.bo <- function(cpm,delta,n,x_bar,s,t){ 
  epslon <- (x_bar-t)/(((n-1)/n)*s)
  v_chapeu=(n*(1+epslon^2)^2) / (1+2*epslon^2) 
  cpmL=cpm*sqrt(qchisq(delta/2, df= v_chapeu)/v_chapeu)
  cpmU=cpm*sqrt(qchisq(1-delta/2, df= v_chapeu)/v_chapeu)
  return(c("Inf" = cpmL,"Sup" = cpmU))
}

#----------------------------------------------------------------------- Simula��o
# Tive que separar em duas para armazenar as amostras tamb�m.

amostras <- function(parametros){
  
  cor <- parametros[[1]]
  sdx <- sqrt(parametros[[2]])
  sdy <- sqrt(parametros[[2]])     # Tome nota: para este experimento, n�o faz diferen�a se a dist de Y tiver par�metros diferentes da dist de X.
  tamanho <- parametros[[3]]
  ciclos <- parametros[[4]]
  mu_x <- parametros[[5]]
  mu_y <- parametros[[5]]    # Tome nota: para este experimento, n�o faz diferen�a se a dist de Y tiver par�metros diferentes da dist de X.
  usl <- parametros[[6]]
  lsl <- parametros[[7]]
  t <- parametros[[8]]
  
  # Primeiro: gerar a amostra. Tanto a AAS quanto a ACO.
  amostra <- ACOeAAS(cor,sdx,sdy,tamanho,ciclos,mu_x,mu_y)

  return(amostra)
}

calculo_ICs <- function(parametros,amostra){
  tamanho <- parametros[[3]]
  ciclos <- parametros[[4]]
  usl <- parametros[[6]]
  lsl <- parametros[[7]]
  t <- parametros[[8]]
  
  # Primeiro: calculando o Cpm estimado.
  cpm <- Cpm(usl,lsl,sd(amostra),mean(amostra)) 
  
  # Terceiro: calculando os IC's.
  Marcucci_e_Beazley_1988 <- c(ic.mb(cpm,.05,tamanho*ciclos),
                               ic.mb(cpm,.10,tamanho*ciclos))
  Chan_et_al_1990 <- c(ic.cxz(cpm,.05,usl,lsl,sd(amostra),mean(amostra),t,tamanho*ciclos),
                       ic.cxz(cpm,.10,usl,lsl,sd(amostra),mean(amostra),t,tamanho*ciclos))
  Zimmer_e_Hubele_1997 <- c(ic.zh(cpm,tamanho*ciclos,sd(amostra),.05,mean(amostra),t),
                            ic.zh(cpm,tamanho*ciclos,sd(amostra),.10,mean(amostra),t))
  Boyles_1991 <- c(ic.bo(cpm,.05,tamanho*ciclos,mean(amostra),sd(amostra),t),
                   ic.bo(cpm,.10,tamanho*ciclos,mean(amostra),sd(amostra),t))
  
  # Quarto: montando a sa�da. Escrevendo o t�tulo e juntando os resultados.
  Titulo <- c(apply(expand.grid(c("Inf","Sup"),
                                c("95%","90%"),
                                c("Marcucci and Beazley (1988)", 
                                  "Chan et.al (1990)", 
                                  "Zimmer and Hubele (1997)", 
                                  "Boyles (1991)")), 1, function(x){paste(x, collapse = " - ")}))
  Conteudo <- c(Marcucci_e_Beazley_1988,
                Chan_et_al_1990,
                Zimmer_e_Hubele_1997,
                Boyles_1991) #%>% matrix(nrow=1)
  #colnames(Conteudo) <- Titulo
  names(Conteudo) <- Titulo
  # Quinto: retornando o resultado
  return(Conteudo)
}


#----------------------------------------------------------------------- Envelopando a fun��o

# Definir par�metros de simula��o

cor=c(1,.5,.8,0) 
varx_e_y=c(1,1.777778,4.020075,15.841192)
tamanho=c(3,5)
ciclos=c(5,10)
# proposta: fixar o tamanho amostral final (n*m) em 15 e 60. S�o quatro cen�rios finais.
# fixando o tamanho em 3 e 5
mu_x_e_y=c(1000,1000.882,1001.738,1003.852)
usl=1008
lsl=992
t=1000


parametros=expand.grid(cor,varx_e_y,tamanho,ciclos,mu_x_e_y,usl,lsl,t) 
parametros=parametros[(parametros[,5]==1000 & parametros[,2]==1.777778) |
                      (parametros[,5]==1000 & parametros[,2]==4.020075) |
                      (parametros[,5]==1000 & parametros[,2]==15.841192) |
                      (parametros[,5]==1000.882 & parametros[,2]==1) |
                      (parametros[,5]==1001.738 & parametros[,2]==1) |
                      (parametros[,5]==1003.852 & parametros[,2]==1),]

# Fun��o replicadora 

replicadora <- function(linha_parametros,n_vezes){ 
  library(magrittr)
  ams <- vector(length = n_vezes, mode ="list") %>% 
         lapply(function(x){amostras(linha_parametros)})
  ICs <- ams %>% lapply(function(x){calculo_ICs(linha_parametros,x)})  %>%
         unlist() %>% matrix(ncol=16,byrow=TRUE)
  return(list(ams,ICs))
}

# rbenchmark::benchmark(replicadora(parametros[1,],10000))

# Agora, replicar para cada cen�rio de simula��o

# Vamos salvar a amostra

set.seed(12345)

n_vezes <- 10000 # reduzi p/ 10^3. Depois valido, a� eu vou para 10^4
resultado <- vector("list",length = nrow(parametros))
backup_ams <- vector("list",length = nrow(parametros))

for(i in 1:nrow(parametros)){
  cenario <- replicadora(parametros[i,],n_vezes)
  resultado[[i]] <- cenario[[2]]
  backup_ams[[i]] <- cenario[[1]]
}

save(resultado,parametros,file="resultado.RData")
save(backup_ams,file="amostras.RData")

# Fazer descritiva e enviar ao Cesar 12/04
# bootstrap param�trico >> estimativa dos par�metro

s (i.e: assumindo normalidade, estimativa p/ m�dia  e vari�ncia) 
# Definir qual intervalo bootstrap: percentil (1ra ordem),
#                                   t-bootstrap (estimativa do Erro, que gera um duplo-bootstrap) e BCA (2nda ordem)
#  Estimativa para a m�dia e dp (bootstrap): p/ cada amostra simulada (na AAS , usas a AAS na simula��o bootstrap. Na ACO, usa ACO na simula��o bootstrap com ordena��o perfeita).
# No caso da ACO : vai ter que usar a ACO 
# Analise de vari�ncia e metodo de compara�oes multiplas para testar o que afeta mais a taxa de cobertura e precis�o
# 

## 1) Pesquisa bilbiogr�fica do bootstrap: o professor compartilhou 3 artigos, favor verificar sobre a cita��o destes artigos feitos.
## 2) no artigo do SDSa, tente entender os tr�s algoritmos diferentes.
# To-Do list:
#
# Testar IC com outros casos da literatura
# Definir os par�metros para o estudo de simula��o (Cpm, m, n e cor)
# Rela��es de equival�ncias entre o Cpm e o Cp para verificar os parametros de simula��o

# Fa





















#----------------------------------------------------------------------- Estimando o intervalo percentil bootstrap


#####################################################################
# BOOTSTRAP PERCENTIL
#####################################################################
# 1) Gerar amostras bootstrap via ACO.
# Resp: feito. A fun��o utilizada � a mesma. Quando for AAS, a cor=0.
# 2) Armazen�-las
# Resp: feito. Ser� armazenado a m�dia de cada reamostra bootstrap. 
# Esse procedimento ser� repetido 1000 vezes (a princ�pio).
# 3) Extrarir os ICs
# Resp: a partir da dist. emp�rica, ser� extra�do o IC
# com base no percentil da mesma (e n�o do intervalo via Wald).
#####################################################################

percentil_boot_sample <- function(cor,tamanho,ciclos,amostra){ # dado que a dist teorica � a dist normal
  if(cor==0){
    AAS <- rnorm(tamanho*ciclos,mean=mean(amostra),sd=sd(amostra))
    return(AAS)
  }else{
    sdx <- sd(amostra) ; sdy <- sdx
    mu_x <- mean(amostra) ; mu_y <- mu_x
    
    sigma <- matrix(c(sdx**2,cor*sdx*sdy,
                      cor*sdx*sdy,sdy**2),byrow=FALSE,ncol=2)
    obj <- MASS::mvrnorm(n=(tamanho**2)*ciclos,mu=c(mu_x,mu_y),
                         Sigma=sigma)
    library(magrittr)
    ACO <- split.data.frame(obj,
                            seq(1,ciclos*tamanho^2,
                                by=tamanho)) %>%
      lapply(function(x){x[order(x[,1]),2]})  %>%
      unlist() %>%
      split(rep(c(1:ciclos),each=tamanho*tamanho)) %>% 
      lapply(function(x){matrix(x,ncol=tamanho,byrow = FALSE)}) %>%
      lapply(diag) %>%
      unlist()
    return(ACO)
  }
}

percentil_boot_IC <- function(boot, amostra, ciclos, tamanho, cor, alfa){ # boot � o n�mero de reamostragens bootstrap
  Lower_B_pos <- round(boot*(alfa/2),0)
  Upper_B_pos <- round(boot*(1-alfa/2),0)
  
  lsl <- 992
  usl <- 1008
  tao <- 1000
  
  # usar CPM AO INVES DA MEDIA. Concertar.
  amostra_boot <- vector(length = boot)
  amostra_boot[1] <- (usl-lsl)/(6*sqrt((sd(amostra)**2)+(mean(amostra)-tao)**2))
  for(i in 2:boot){
    iteracao <- percentil_boot_sample(cor, tamanho, ciclos, amostra)
    amostra_boot[i] <- (usl-lsl)/(6*sqrt((sd(iteracao)**2)+(mean(iteracao)-tao)**2))
  }
  amostra_boot <- sort(amostra_boot)
  return(c("LCB"=amostra_boot[Lower_B_pos],"UCB"=amostra_boot[Upper_B_pos]))
}

microbenchmark::microbenchmark({ # pRONTO. Vamos estimar o tempo com a Toy Sample
  percentil_boot_IC(99,as.numeric(unlist(backup_ams[[1]][1])),3,3,1,0.05)
})


percentil_boot_IC_Wald <- function(boot, amostra, ciclos, tamanho, cor, alfa){ # boot � o n�mero de reamostragens bootstrap
  
  lsl <- 992
  usl <- 1008
  tao <- 1000
  
  # usar CPM AO INVES DA MEDIA. Concertar.
  amostra_boot <- vector(length = boot)
  amostra_boot[1] <- (usl-lsl)/(6*sqrt((sd(amostra)**2)+(mean(amostra)-tao)**2))
  for(i in 2:boot){
    iteracao <- percentil_boot_sample(cor, tamanho, ciclos, amostra)
    amostra_boot[i] <- (usl-lsl)/(6*sqrt((sd(iteracao)**2)+(mean(iteracao)-tao)**2))
  }
  amostra_boot <- sort(amostra_boot)
  return(c("LCB"=mean(amostra_boot)-qnorm(1-alfa/2)*sd(amostra_boot),
           "UCB"=mean(amostra_boot)+qnorm(1-alfa/2)*sd(amostra_boot)))
}


microbenchmark::microbenchmark({ # Estimando o tempo com a Toy Sample
  percentil_boot_IC(99,as.numeric(unlist(backup_ams[[1]][1])),3,5,1,0.05)
})

# Ter feito a Toy Sample e mostrar o tempo estimado (se der para rodar um cen�rio, fazer)

# Unit: milliseconds
# expr      min       lq
# {     percentil_boot_IC(99, as.numeric(unlist(backup_ams[[1]][1])),          3, 5, 1, 0.05) } 114.0586 117.1912
# mean   median       uq      max neval
# 128.9917 119.7476 129.5701 344.8911   100

# Resp: usando n*m=15, vimos que cada reamostragem bootstrap leva em m�dia 0.001292929 segundos. Sendo 96 cen�rios,
# considerando 10000 reamostragens bootstraps, estimamos que o algoritimo leve 143.6588 dias para
# calcular os ICs bootstraps


# Ter feito a Toy Sample e mostrar o tempo estimado (se der para rodar um cen�rio, fazer)
# Definir os cenarios a serem simulados (d� prefer�ncia a aqueles mais contrastantes)
# Sugestao: fixar n e variar m
# Estrutura do relat�rio final
# Espectativa de subcobertura

# Gerar IC bootstrap pelo m�todo Wald (confirmar no material de EstCompII como fazer)
# Caso ache outra constru��o (e.g: BCA), implementar tamb�m.
# Fa�a outra fun��o para AAS, por conta do esfor�o computacional.
# 4 cen�rios de n*m que resultam {15,60} x 6 combina��es de Cpm que resultam {2; 1,33; 0,67} x Correla��es 1, 0.8 e 0 (AAS)

































# 29/06

# Os Intervalos bootstrap foram gerados (levou 5 dias na minha m�quina)
# e os c�digos foram transferidos para este script, e no "Codigos2.R"
# ser� analisado os dados (at� para n�o misturar uma coisa com a outra)

#####################################################################
# Estimar B e M. Identificar cen�rios
#####################################################################
# O objetivo desta etapa � identificar os 
# cen�rios ideais para fazer os intervalos bootstraps,
# e definir os tamanhos de B e M ideais
#####################################################################

#####################################################################
# Rodada 1
# B: {1000 e 2000} ; M: {3000 e 5000}
#####################################################################

# B <- c(999, 1999) # n�mero de reamostras bootstrap
# 
# M <- c(3000,5000) # n�mero entre 0 e 10000 de amostras a serem utilizadas.
# 
# rodada1 <- expand.grid(B,M)
# 
# toy_sample <- unlist(backup_ams[[1]])

### Temos que fazer um for, que
### teste cada linha daquela rodada1, e armazene o output
### da fun��o microbenchmark. Vide exemplo abaixo.

# teste_tempo1 <- vector(mode="list",length = nrow(rodada1))

### Beleza! Vamos fazer o teste: tamanho (n*m) � 15, corr=1 e alfa 0.05. Intervalo percentil,
### n�o o de Wald

# for(i in 1:nrow(rodada1)){
#   teste_tempo1[[i]] <- {
#     start_time <- Sys.time()
#     percentil_boot_IC(rodada1[i,1],
#                       toy_sample,
#                       rodada1[i,2],
#                       15,
#                       1,
#                       0.05)
#     end_time <- Sys.time()
#     end_time-start_time
#   }
# }

### �timo. Agora, vamos ver o de Wald

# teste_tempo2 <- vector(mode="list",length = nrow(rodada1))
# 
# for(i in 1:nrow(rodada1)){
#   teste_tempo2[[i]] <- {
#     start_time <- Sys.time()
#     percentil_boot_IC_Wald(rodada1[i,1],
#                            toy_sample,
#                            rodada1[i,2],
#                            15,
#                            1,
#                            0.05)
#     end_time <- Sys.time()
#     end_time-start_time
#   }
# }

# N�o deu certo o teste, vamos ter que pegar o cen�rio 
# mais econ�mico e seguir com ele.

### Microbenchmark est� levando muito tempo... Por isso alterei

# 14/06

## Parte escrita
### Baixei o livro do Brian Manly!
### Baixar o livro do Daveson e Hinkley - Baixou
### Baixar o livro do Efron e Tibshirani _ Baixou

### terminar o relat�rio final (antes do resultados e discuss�es)

### O que que tem de estima��o intervalar
### p/ outros valores no Boootstrap (suspeita: subcobertura)
### estima��o do (intraclass coeficient correlation)
###  

## Parte experimento
### Programa��o paralela!
### Fechar as simulacoes (mesmo com o B e o M pequenos)
### Come�ar a edi��o dos dados pois os mesmos v�o tomar muito tempo.
### Edi��o dos dados (ver trabalho do Vinicius como referencia)
### Resultados e discussoes: ver obj geral e especifico do projeto


### Vamos gerar agora os intervalos via bootstrap.

#### PROGRAMA��O PARALELA

require(doParallel)
ncl <- detectCores() # Checa quantos n�cleos existem na m�quina
cl <- makeCluster(ncl)
registerDoParallel(cl) # Registra os n�cleos a serem utilizados

#### PAR�METROS DE SIMULA��O E GERA��O DE RESULTADOS

IC_bootstrap_percentil <- vector(mode="list",length = nrow(parametros))
IC_bootstrap_Wald <- vector(mode="list",length = nrow(parametros))

B <- 999
M <- 3000

for(j in 1:nrow(parametros)){
  IC_bootstrap_percentil[[j]] <- foreach(i=1:M)%dopar%percentil_boot_IC(boot=B,
                                                                        amostra=backup_ams[[j]][[i]],
                                                                        ciclos=parametros[j,4],
                                                                        tamanho=parametros[j,3],
                                                                        cor=parametros[j,1],
                                                                        alfa=0.05)
  IC_bootstrap_Wald[[j]] <- foreach(i=1:M)%dopar%percentil_boot_IC_Wald(boot=B,
                                                                        amostra=backup_ams[[j]][[i]],
                                                                        ciclos=parametros[j,4],
                                                                        tamanho=parametros[j,3],
                                                                        cor=parametros[j,1],
                                                                        alfa=0.05)
} 


save(IC_bootstrap_percentil,parametros,
     file="IC_bootstrap_percentil.RData")
save(IC_bootstrap_Wald,
     file="IC_bootstrap_Wald.RData")

### ENCERRAR PROGRAMA��O PARALELA

stopCluster(cl)









































# 01/07
# Deu pau no IC Wald. Os limites est�o iguais, vamos corrigir.
# Corrigi a fun��o no come�o, vou fazer outra sess�o do R, e copiar e colar o "for" adaptado aqui
# Depois, preciso corrigir o "for" logo acima

### Vamos gerar agora os intervalos via bootstrap.

#### PROGRAMA��O PARALELA

require(doParallel)
ncl <- detectCores() # Checa quantos n�cleos existem na m�quina
cl <- makeCluster(ncl)
registerDoParallel(cl) # Registra os n�cleos a serem utilizados

#### PAR�METROS DE SIMULA��O E GERA��O DE RESULTADOS

IC_bootstrap_Wald <- vector(mode="list",length = nrow(parametros))

B <- 999
M <- 3000

for(j in 1:nrow(parametros)){
  IC_bootstrap_Wald[[j]] <- foreach(i=1:M)%dopar%percentil_boot_IC_Wald(boot=B,
                                                                        amostra=backup_ams[[j]][[i]],
                                                                        ciclos=parametros[j,4],
                                                                        tamanho=parametros[j,3],
                                                                        cor=parametros[j,1],
                                                                        alfa=0.05)
} 

save(IC_bootstrap_Wald,
     file="IC_bootstrap_Wald.RData")

### ENCERRAR PROGRAMA��O PARALELA

stopCluster(cl)

