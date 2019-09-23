# Esse script conterá a análise pós-simulação


# Resumo: mudar a parte (ICs da verossimilhança) para IC's construidas sob o delineamento AAS

# 1ra familia de ICs: não teve um grande ganho de precisão para as amostras.

# Teste de hipóteses para a taxa de cobertura: H0 = Tx de cobert = 0.95 >> teste para proporção

#----------------------------------------------------------------------- Função para gerar taxa de cobertura pós-simulação

load(file="amostras.RData")
load(file="resultado.RData")

Cpm <- function(vetor){ # primeiro, identificar o Cpm de cada cenário
  usl <- 1008
  lsl <- 992
  sigma <- sqrt(vetor[2])
  mu <- vetor[5]
  tao <- mean(c(usl,lsl))
  return((usl-lsl)/(6*sqrt((sigma**2)+(mu-tao)**2)))
}

plot(unlist(lapply(resultado,mean)),ylab="Média dos dados", xlab="Índice")

Cpm_s <- round(apply(parametros,1,Cpm),2)

taxa_de_cobertura <- vector(mode="list",length = length(resultado))
precisao <- vector(mode="list",length = length(resultado))

gabarito_para_o_for <- data.frame(seq(1,ncol(resultado[[1]]),by=2),1:8)

# Faremos da seguinte forma: para cada um dos 96 dataframes de 16x10000, vamos separá-los em 
# 8 dataframes 2x10000. E então, a função taxa_de_cobertura_fun e a precisao_fun irá calcular
# para cada um desses 8 dataframes as nossas medidas de interesse, e depois juntaremos os va-
# lores no final.

taxa_de_cobertura_fun <- function(data_frame,valor){
  resultado <- apply(data_frame,1,function(x){return((x[1] < valor)*(x[2] > valor))})
  return(mean(resultado))
}

teste <- vector(mode="list",length = 8)

for(i in seq(1,15,2)){
  teste[[gabarito_para_o_for[gabarito_para_o_for[,1]==i,2]]] <- resultado[[1]][,c(i,i+1)]
}

##################################################################################################################

for(i in 1:length(resultado)){
  Cpm_da_vez <- Cpm_s[i]
  cenario_da_vez <- as.data.frame(resultado[[i]])
  base_binaria <- vector(mode = "list" , length = ncol(cenario_da_vez)/2)
  base_amplitude <- vector(mode = "list" , length = ncol(cenario_da_vez)/2)
  for(j in seq(1,ncol(cenario_da_vez),by=2)){ # vamos fazer uma base binária onde 1 é que o intervalo conteve o Cpm e 0 é que o intervalo não conteve o Cpm.
    base_binaria[[gabarito_para_o_for[gabarito_para_o_for[,1]==j,2]]] = apply(dplyr::select(cenario_da_vez,c(j,j+1)),1,function(x){ifelse(x[1] < Cpm_da_vez && x[2] > Cpm_da_vez,1,0)})
    base_amplitude[[gabarito_para_o_for[gabarito_para_o_for[,1]==j,2]]] = apply(dplyr::select(cenario_da_vez,c(j,j+1)),1,function(x){return(abs(x[2]-x[1]))})
  }
  taxa_de_cobertura[[i]] <- unlist(lapply(base_binaria,mean))
  precisao[[i]] <- unlist(lapply(base_amplitude,mean))
}

taxas_de_cobertura <- data.frame(Cpm_s,parametros,matrix(unlist(taxa_de_cobertura),ncol=8,byrow=TRUE))
precisoes <- data.frame(Cpm_s,parametros,matrix(unlist(precisao),ncol=8,byrow = TRUE))

#----------------------------------------------------------------------- Uma vez gerada as taxas de cobertura, vamos proceder com as análises!

##---------------------------------------------------------------------- 1) Auditoria: gráfico "Resultado" deve estar próximo do "Esperado"

names(taxas_de_cobertura) <- c("Cpm_s",c("cor","varx_e_y","tamanho","ciclos","mu_x_e_y","usl","lsl","T"),apply(expand.grid(
                                                 c("95%","90%"),
                                                 c("Marcucci and Beazley (1988)", 
                                                   "Chan et.al (1990)", 
                                                   "Zimmer and Hubele (1997)", 
                                                   "Boyles (1991)")), 1, function(x){paste(x, collapse = " - ")}))

names(precisoes) <- c("Cpm_s",c("cor","varx_e_y","tamanho","ciclos","mu_x_e_y","usl","lsl","T"),apply(expand.grid(
                                                  c("95%","90%"),
                                                  c("Marcucci and Beazley (1988)", 
                                                    "Chan et.al (1990)", 
                                                    "Zimmer and Hubele (1997)", 
                                                    "Boyles (1991)")), 1, function(x){paste(x, collapse = " - ")}))

# save(list=c("taxa_de_cobertura","precisao"),file="resultados_tabuldados.RData")

# X11()
# par(mfrow=c(2,1))
# plot(Cpm_s, xlab="", ylab="Cpm estimado",main="Esperado",ylim = c(0,2.5))
# plot(unlist(lapply(resultado,mean)), xlab="",ylab="Média dos dados",main="Resultado",ylim = c(0,2.5))
# abline(h=2,col="red")
# abline(h=1.33,col="red")
# abline(h=0.67,col="red")

##---------------------------------------------------------------------- Feito

##---------------------------------------------------------------------- 2) Montar a tabela usando o pacote xtable

tamanho <- unique(parametros$Var3)

ciclos <- unique(parametros$Var4)

corr <- unique(parametros$Var1)

expand1 <- expand.grid(corr,ciclos,tamanho)

tabela <- data.frame("n, m" = paste(expand.grid(corr,ciclos,tamanho)[,3],
                                    expand.grid(corr,ciclos,tamanho)[,2],sep=", "),
                     "corr"= expand.grid(corr,ciclos,tamanho)[,1])

# Verificar ordem das linhas por conta dos ciclos! Obs: tive que reordenar as colunas
# das tabelas para aparecerem no formato que a gente quer

# Cpm 2
## mu=1000 e var=1.777778
## mu=1000.882 e var=1

Cpm2 <- data.frame(tabela,
                   precisoes[(precisoes$Cpm_s==2 & precisoes$mu_x_e_y==1000),
                             c(11,13,15,17,10,12,14,16)],
                   taxas_de_cobertura[(taxas_de_cobertura$Cpm_s==2 & taxas_de_cobertura$mu_x_e_y==1000),
                                      c(11,13,15,17,10,12,14,16)],
                   precisoes[(precisoes$Cpm_s==2 & precisoes$varx_e_y==1),
                             c(11,13,15,17,10,12,14,16)],
                   taxas_de_cobertura[(taxas_de_cobertura$Cpm_s==2 & taxas_de_cobertura$varx_e_y==1),
                                      c(11,13,15,17,10,12,14,16)])


# Arredondando para inserir a tabela no trabalho

# Cpm2[,-c(1,2)] <- round(Cpm2[,-c(1,2)],2) 

names(Cpm2)[-c(1,2)] <- apply(expand.grid(c("Marcucci and Beazley (1988)", 
                                            "Chan et.al (1990)", 
                                            "Zimmer and Hubele (1997)", 
                                            "Boyles (1991)"),
                                          c("90","95"),c("Precisoes","Taxas de cobertura"),
                                          c("mu=1000","sigma=1"))[c(4,3,2,1)],
                              1,function(x){paste(x, collapse = " - ")})

rownames(Cpm2) <- NULL

#tab2 <- xtable::xtable(Cpm2)

# Cpm 1.33
## mu=1000 e var=4.020075
## mu=1001.738 e var=1

Cpm1.33 <- data.frame(tabela,
                   precisoes[(precisoes$Cpm_s==1.33 & precisoes$mu_x_e_y==1000),
                             c(11,13,15,17,10,12,14,16)],
                   taxas_de_cobertura[(taxas_de_cobertura$Cpm_s==1.33 & taxas_de_cobertura$mu_x_e_y==1000),
                                      c(11,13,15,17,10,12,14,16)],
                   precisoes[(precisoes$Cpm_s==1.33 & precisoes$varx_e_y==1),
                             c(11,13,15,17,10,12,14,16)],
                   taxas_de_cobertura[(taxas_de_cobertura$Cpm_s==1.33 & taxas_de_cobertura$varx_e_y==1),
                                      c(11,13,15,17,10,12,14,16)])


# Arredondando para inserir a tabela no trabalho

# Cpm1.33[,-c(1,2)] <- round(Cpm1.33[,-c(1,2)],2) 


names(Cpm1.33)[-c(1,2)] <- apply(expand.grid(c("Marcucci and Beazley (1988)", 
                                            "Chan et.al (1990)", 
                                            "Zimmer and Hubele (1997)", 
                                            "Boyles (1991)"),
                                          c("90%","95%"),c("Precisoes","Taxas de cobertura"),
                                          c("mu=1000","sigma=1"))[c(4,3,2,1)],
                              1,function(x){paste(x, collapse = " - ")})

#tab1.33 <- xtable::xtable(Cpm1.33)

# Cpm 0.67
## mu=1000 e var=15.841192
## mu=1003.852 e var=1

Cpm0.67 <- data.frame(tabela,
                   precisoes[(precisoes$Cpm_s==0.67 & precisoes$mu_x_e_y==1000),
                             c(11,13,15,17,10,12,14,16)],
                   taxas_de_cobertura[(taxas_de_cobertura$Cpm_s==0.67 & taxas_de_cobertura$mu_x_e_y==1000),
                                      c(11,13,15,17,10,12,14,16)],
                   precisoes[(precisoes$Cpm_s==0.67 & precisoes$varx_e_y==1),
                             c(11,13,15,17,10,12,14,16)],
                   taxas_de_cobertura[(taxas_de_cobertura$Cpm_s==0.67 & taxas_de_cobertura$varx_e_y==1),
                                      c(11,13,15,17,10,12,14,16)])

# Arredondando para inserir a tabela no trabalho

# Cpm0.67[,-c(1,2)] <- round(Cpm0.67[,-c(1,2)],2) 


names(Cpm0.67)[-c(1,2)] <- apply(expand.grid(c("Marcucci and Beazley (1988)", 
                                            "Chan et.al (1990)", 
                                            "Zimmer and Hubele (1997)", 
                                            "Boyles (1991)"),
                                          c("90%","95%"),c("Precisoes","Taxas de cobertura"),
                                          c("mu=1000","sigma=1"))[c(4,3,2,1)],
                              1,function(x){paste(x, collapse = " - ")})

#tab0.67 <- print(xtable::xtable(Cpm0.67),include.rownames=FALSE)

#----------------------------------------------------------------------- Estimando o intervalo percentil bootstrap


#####################################################################
# BOOTSTRAP PERCENTIL
#####################################################################
# 1) Gerar amostras bootstrap via ACO.
# Resp: feito. A função utilizada é a mesma. Quando for AAS, a cor=0.
# 2) Armazená-las
# Resp: feito. Será armazenado a média de cada reamostra bootstrap. 
# Esse procedimento será repetido 1000 vezes (a princípio).
# 3) Extrarir os ICs
# Resp: a partir da dist. empírica, será extraído o IC
# com base no percentil da mesma (e não do intervalo via Wald).
#####################################################################

percentil_boot_sample <- function(cor,tamanho,ciclos,amostra){ # dado que a dist teorica é a dist normal
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

percentil_boot_IC <- function(boot, amostra, ciclos, tamanho, cor, alfa){ # boot é o número de reamostragens bootstrap
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

# microbenchmark::microbenchmark({ # pRONTO. Vamos estimar o tempo com a Toy Sample
#   percentil_boot_IC(99,as.numeric(unlist(backup_ams[[1]][1])),3,5,1,0.05)
# })

percentil_boot_IC_Wald <- function(boot, amostra, ciclos, tamanho, cor, alfa){ # boot é o número de reamostragens bootstrap
  
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


# microbenchmark::microbenchmark({ # Estimando o tempo com a Toy Sample
#   percentil_boot_IC(99,as.numeric(unlist(backup_ams[[1]][1])),3,5,1,0.05)
# })

# Ter feito a Toy Sample e mostrar o tempo estimado (se der para rodar um cenário, fazer)

# Unit: milliseconds
# expr      min       lq
# {     percentil_boot_IC(99, as.numeric(unlist(backup_ams[[1]][1])),          3, 5, 1, 0.05) } 114.0586 117.1912
# mean   median       uq      max neval
# 128.9917 119.7476 129.5701 344.8911   100

# Resp: usando n*m=15, vimos que cada reamostragem bootstrap leva em média 0.001292929 segundos. Sendo 96 cenários,
# considerando 10000 reamostragens bootstraps, estimamos que o algoritimo leve 143.6588 dias para
# calcular os ICs bootstraps


# Definir os cenarios a serem simulados (dê preferência a aqueles mais contrastantes)
# Sugestao: fixar n e variar m
# Estrutura do relatório final
# Espectativa de subcobertura

# Gerar IC bootstrap pelo método Wald (confirmar no material de EstCompII como fazer)
# Caso ache outra construção (e.g: BCA), implementar também.
# Faça outra função para AAS, por conta do esforço computacional.
# 4 cenários de n*m que resultam {15,60} x 6 combinações de Cpm que resultam {2; 1,33; 0,67} x Correlações 1, 0.8 e 0 (AAS)

# Avaliando ICs 95%

# View(precisoes[precisoes$Cpm_s==2 & precisoes$cor==1,-c(7,8,9,11,13,15,17)]) #Cpm2_corr1
# View(precisoes[precisoes$Cpm_s==2 & precisoes$cor==0,-c(7,8,9,11,13,15,17)]) #Cpm2_corr0

# 29/05


# Livro: Brian Manly para método bootstrap >> objetivo, 
# encontrar embasamento para encontrar o B e M.
# Observar os artigos que listamos na reunião (e.g: B=1000 e M=3000)

# 1) Estudar combinaçoes de B e M, verificar o tempo que leva para cada um.

# simulaçoes sequenciais de B para avaliar os quantis!
# 2) estudo de convergência: qual o tamanho ideal de B?

# 3) verificar o nro de processadores para prog paralela. E preparar um código p/ prog paralela.



















# 29/06

# Vamos analisar os intervalos bootstrap!
load(file="IC_bootstrap_percentil.RData")
load(file="IC_bootstrap_Wald.RData")

IC_bootstrap_percentil <- lapply(IC_bootstrap_percentil,
                                 function(x){matrix(unlist(x),byrow=TRUE,ncol=2)}) # não rode duas vezes! Se fizer, rode o "load" junto
IC_bootstrap_Wald <- lapply(IC_bootstrap_Wald,
                            function(x){matrix(unlist(x),byrow=TRUE,ncol=2)}) # não rode duas vezes! Se fizer, rode o "load" junto

taxa_de_cobertura_boot <- vector(mode="list",length = length(resultado))
precisao_boot <- vector(mode="list",length = length(resultado))

# Precisamos adaptar aqui...

## FIX THIS!!!!

taxa_de_cobertura_boot_fun <- function(lista1,lista2,valor){
  resultado1 <- mean(apply(lista1,1,function(x){(x[1] < valor)*(x[2] > valor)}))
  resultado2 <- mean(apply(lista2,1,function(x){(x[1] < valor)*(x[2] > valor)}))
  return(c(resultado1,resultado2))
}

precisao_boot_fun <- function(lista1,lista2){
  resultado1 <- mean(apply(lista1,1,function(x){x[2]-x[1]}))
  resultado2 <- mean(apply(lista2,1,function(x){x[2]-x[1]}))
  return(c(resultado1,resultado2))
}

##################################################################################################################

for(i in 1:nrow(parametros)){
  taxa_de_cobertura_boot[[i]] <- taxa_de_cobertura_boot_fun(IC_bootstrap_percentil[[i]],
                                                            IC_bootstrap_Wald[[i]],
                                                            as.numeric(Cpm_s[i]))
  precisao_boot[[i]] <- precisao_boot_fun(IC_bootstrap_percentil[[i]],
                                          IC_bootstrap_Wald[[i]])
}

taxas_de_cobertura_boot <- data.frame(Cpm_s,parametros,matrix(unlist(taxa_de_cobertura_boot),ncol=2,byrow=TRUE))
precisoes_boot <- data.frame(Cpm_s,parametros,matrix(unlist(precisao_boot),ncol=2,byrow = TRUE))

#----------------------------------------------------------------------- Uma vez gerada as taxas de cobertura, vamos proceder com as análises!

##---------------------------------------------------------------------- 1) Auditoria: gráfico "Resultado" deve estar próximo do "Esperado"

names(taxas_de_cobertura_boot) <- c("Cpm_s","cor","varx_e_y","tamanho","ciclos","mu_x_e_y","usl","lsl","T",
                                    "IC bootstrap percentil","IC bootstrap Wald")

names(precisoes_boot) <- c("Cpm_s","cor","varx_e_y","tamanho","ciclos","mu_x_e_y","usl","lsl","T",
                           "IC bootstrap percentil","IC bootstrap Wald")




###############################################################
# ANÁLISES
###############################################################

# Vamos acrescentar uma coluna com tamanho*ciclos
taxas_de_cobertura_boot <- dplyr::mutate(taxas_de_cobertura_boot,
                                         nm=taxas_de_cobertura_boot$tamanho*taxas_de_cobertura_boot$ciclos)
precisoes_boot <- dplyr::mutate(precisoes_boot,
                                nm=precisoes_boot$tamanho*precisoes_boot$ciclos)

# montar data frames para as análises

## Taxas de cobertura

bd_para_melt <- data.frame(taxas_de_cobertura_boot[,c(1,2,3,6,12,10,11)],precisoes_boot[,c(10,11)])
## Vamos transformar as col media e var em uma única, dizendo (shift média, shift var)

bd_para_melt$shift <- apply(bd_para_melt[,c(3,4)],1,function(x){ifelse(x[1]==1,"Media","Variância")})

bd_para_melt <- bd_para_melt[,c(1,2,10,5,6,7,8,9)]

gabarito <- paste(bd_para_melt[,1],bd_para_melt[,2],bd_para_melt[,4])
teste <- bd_para_melt
teste$gabarito <- gabarito

#1) Será que os "Shifts" causam algum impacto nas mensurações?
# fazer testes de hipóteses pode não ser o mais adequado por conta
# do tamanho das amostras. No entanto, testes de aleatorização podem ser adequados! Mas.. vamos fazer uma descritiva antes


X11()
par(mfrow=c(2,2))

shift_media <- teste[teste$shift=="Media",c(5,9)]
shift_var <- teste[teste$shift=="Variância",c(5,9)]
shift <- merge(shift_media,shift_var,gabarito,by.x="gabarito",by.y="gabarito")
hist(shift[,2]-shift[,3], main = "Tx Cobertra - Percentil")

shift_media <- teste[teste$shift=="Media",c(6,9)]
shift_var <- teste[teste$shift=="Variância",c(6,9)]
shift <- merge(shift_media,shift_var,gabarito,by.x="gabarito",by.y="gabarito")
hist(shift[,2]-shift[,3], main = "Tx Cobertra - Wald")

shift_media <- teste[teste$shift=="Media",c(7,9)]
shift_var <- teste[teste$shift=="Variância",c(7,9)]
shift <- merge(shift_media,shift_var,gabarito,by.x="gabarito",by.y="gabarito")
hist(shift[,2]-shift[,3], main = "Precisão - Percentil")

shift_media <- teste[teste$shift=="Media",c(8,9)]
shift_var <- teste[teste$shift=="Variância",c(8,9)]
shift <- merge(shift_media,shift_var,gabarito,by.x="gabarito",by.y="gabarito")
hist(shift[,2]-shift[,3], main = "Precisão - Wald")

par(mfrow=c(1,1))

# De fato, esse shift traz uma certa diferença, ainda que pequena
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Vamos tentar fazer os lolliplots!

#View(bd_para_melt)

names(bd_para_melt) <- c("Cpms","cor","shift","nm",
                         "Tx_Cobertura_Percentil",
                         "Tx_Cobertura_Wald",
                         "Precisao_Percentil",
                         "Precisao_Wald")
# Precisamos fazer um truque para parear os tipos diferentes de bootstrap no mesmo gráfico...
# E esse truque vai ser criar um factor para a variável "x"

X <- paste(c(unique(bd_para_melt$nm),
             unique(bd_para_melt$nm)),
           rep(c("percentil","wald"),each=4),sep=" ")
X <- X[c(1,5,2,6,3,7,4,8)]
niveis <- c(X[1],X[2],"e1","e2",
            X[3],X[4],"e3","e4",
            X[5],X[6],"e5","e6",
            X[7],X[8])
#x <- factor(x,levels = niveis)
# plot(x,1:8) aparentemente deu certo!!

x <- factor(paste(c(bd_para_melt$nm,
                    bd_para_melt$nm),
                  rep(c("percentil","wald"),
                      each=length(bd_para_melt$nm)),
                  sep=" "), levels=niveis)
bd_plot <- data.frame(x,rbind(bd_para_melt[,c(1,2,3,4)],
                              bd_para_melt[,c(1,2,3,4)]),
                      "Tx_Cobertura"=c(bd_para_melt[,5],bd_para_melt[,6]),
                      "Precisao"=c(bd_para_melt[,7],bd_para_melt[,8]))

###########################################################
# Gráfico Taxas de Cobertura

X11()
layout(matrix(c( 1, 1, 2, 2, 3, 3,13,
                 1, 1, 2, 2, 3, 3,13,
                 4, 4, 5, 5, 6, 6,13,
                 4, 4, 5, 5, 6, 6,13,
                 7, 7, 8, 8, 9, 9,14,
                 7, 7, 8, 8, 9, 9,14,
                 10,10,11,11,12,12,14,
                 10,10,11,11,12,12,14),ncol=8,byrow=FALSE))

bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==0,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=0.67, corr=0 (AAS)")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==0,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=1.33, corr=0 (AAS)")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==0,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=2, corr=0 (AAS)")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))

bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==0.5,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=0.67, corr=0.5")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==0.5,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=1.33, corr=0.5")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==0.5,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=2, corr=0.5")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))

bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==0.8,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=0.67, corr=0.8")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==0.8,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=1.33, corr=0.8")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==0.8,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=2, corr=0.8")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))
bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==1,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=0.67, corr=1")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==1,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=1.33, corr=1")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==1,c(1,4,5,6)]

x <- bd_teste$x
y <- bd_teste$Tx_Cobertura
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Tx_Cobertura),xaxt="n",pch=z,col=w,cex=2,ylab="Taxa de cobertura",xlab="nm",main="Cpm=2, corr=1")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))

par(mai=c(0,0,0,0))
plot.new()
legend("center",legend=c("Bootstrap percentil",
                         "Bootstrap Wald"),ncol=2,fill=c("blue","red"),cex=1.25)
par(mai=c(0,0,0,0))
plot.new()
legend("center",legend=c("V: shift na variância",
                         "M: shift na média"),ncol=2,cex=1.25)

par(mfrow=c(1,1))

###########################################################
# Gráfico precisões

X11()
layout(matrix(c( 1, 1, 2, 2, 3, 3,13,
                 1, 1, 2, 2, 3, 3,13,
                 4, 4, 5, 5, 6, 6,13,
                 4, 4, 5, 5, 6, 6,13,
                 7, 7, 8, 8, 9, 9,14,
                 7, 7, 8, 8, 9, 9,14,
                10,10,11,11,12,12,14,
                10,10,11,11,12,12,14),ncol=8,byrow=FALSE))

bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==0,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=0.67, corr=0 (AAS)")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==0,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=1.33, corr=0 (AAS)")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==0,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=2, corr=0 (AAS)")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))

bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==0.5,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=0.67, corr=0.5")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==0.5,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=1.33, corr=0.5")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==0.5,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=2, corr=0.5")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))

bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==0.8,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=0.67, corr=0.8")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==0.8,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=1.33, corr=0.8")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==0.8,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=2, corr=0.8")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))
bd_teste <- bd_plot[bd_plot$Cpms==0.67 & bd_plot$cor==1,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=0.67, corr=1")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==1.33 & bd_plot$cor==1,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=1.33, corr=1")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))


bd_teste <- bd_plot[bd_plot$Cpms==2 & bd_plot$cor==1,c(1,4,5,7)]

x <- bd_teste$x
y <- bd_teste$Precisao
z <- bd_teste$shift# tem que ser 17 e 16 p/ dar boa!
w <- rep(c("red","blue"),each=8)
#
plot(as.numeric(x),y,cex.main=2,cex.lab=1.5,ylim=range(bd_plot$Precisao),xaxt="n",pch=z,col=w,cex=2,ylab="Precisão",xlab="nm",main="Cpm=2, corr=1")

segments(as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][1]}),
         as.numeric(x[1:length(x)]),
         sapply(1:length(x),function(i){y[sapply(1:length(x),function(j){x[j]==x})[,i]][2]}),lwd=3,lty=3)
axis(1,at=as.numeric(x)[1:length(x)],labels=c("  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50",
                                              "  15","  25","  30","  50"))
par(mai=c(0,0,0,0))
plot.new()
legend("center",legend=c("Bootstrap percentil",
                "Bootstrap Wald"),ncol=2,fill=c("blue","red"),cex=1.25)
par(mai=c(0,0,0,0))
plot.new()
legend("center",legend=c("V: shift na variância",
                         "M: shift na média"),ncol=2,cex=1.25)


#https://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r





























###################################################
# Tabela
# Bootstrap >> Shift Média
print(xtable::xtable({bd_para_melt$nm <- as.integer(bd_para_melt$nm);bd_para_melt[bd_para_melt$shift=="Media",c(1,4,2,5,6,7,8)]}), include.rownames = FALSE,digits=c(2,0,1,4,4,4,4))
# Bootstrap >> Shift Variancia
print(xtable::xtable({bd_para_melt$nm <- as.integer(bd_para_melt$nm);bd_para_melt[bd_para_melt$shift=="Variância",c(1,4,2,5,6,7,8)]},digits=c(2,0,1,4,4,4,4)), include.rownames = FALSE)












########################
# aNÁLISE DE REGRESSÃO : TX DE PRECISÃO E COBERTURA EXPLICADA PELOS PARÂMETROS
# Vamos fazer uma análise de regressão para os dados.

# Independente da lista dos parâmetros vamos ver a dist por tipo de intervalo, da precisão e taxa de cobertura  
# Taxas de cobertura
boxplot_das_tx <- data.frame("IC_MB_1988"=taxas_de_cobertura$`95% - Marcucci and Beazley (1988)`,
                             "IC_C_1990"=taxas_de_cobertura$`95% - Chan et.al (1990)`,
                             "IC_ZH_1997"=taxas_de_cobertura$`95% - Zimmer and Hubele (1997)`,
                             "IC_BO_1991"=taxas_de_cobertura$`95% - Boyles (1991)`,
                             "IC_boot_perc"=taxas_de_cobertura_boot$`IC bootstrap percentil`,
                             "IC_boot_Wald"=taxas_de_cobertura_boot$`IC bootstrap Wald`)
X11()
par(cex.axis=1,mar=c(3,7,2.5,1)) 
boxplot(boxplot_das_tx, main="Taxa de cobertura",
        horizontal=TRUE,
        las=2,
        names=c("Marcucci e\nBeazley (1988)",
                               "Chan, Xiong e\nZhang (1990)",
                               "Boyles (1991)",
                               "Zimmer e\nHubele (1997)",
                               "Bootstrap\nPercentil",
                               "Bootstrap\nWald")) # Intervalos bootstrap mais consistentes. E uma super cobertura do MB
abline(v=0.95,col="red",lwd=2)

# Precisões
boxplot_das_prec <- data.frame("IC_MB_1988"=precisoes$`95% - Marcucci and Beazley (1988)`,
                             "IC_C_1990"=precisoes$`95% - Chan et.al (1990)`,
                             "IC_ZH_1997"=precisoes$`95% - Zimmer and Hubele (1997)`,
                             "IC_BO_1991"=precisoes$`95% - Boyles (1991)`,
                             "IC_boot_perc"=precisoes_boot$`IC bootstrap percentil`,
                             "IC_boot_Wald"=precisoes_boot$`IC bootstrap Wald`)
X11()
par(cex.axis=1,mar=c(3,7,2.5,1)) 
boxplot(boxplot_das_prec, main="Precisão",
        horizontal=TRUE,
        las=2,
        names=c("Marcucci e\nBeazley (1988)",
                "Chan, Xiong e\nZhang (1990)",
                "Boyles (1991)",
                "Zimmer e\nHubele (1997)",
                "Bootstrap\nPercentil",
                "Bootstrap\nWald")) # Intervalos bootstrap mais consistentes. E uma super cobertura do MB
#abline(h=0.95,col="red")


# Fazendo as regressoes

bd_reg_prec <- read.csv("precisoes_regressao.csv",header=TRUE,sep=";")[,-c(4,5)]
lm_prec1 <- lm(Precisoes~.,data=bd_reg_prec)
summary(lm_prec1)
X11()
par(mfrow=c(2,2))
plot(lm_prec1) # talvez... passar a raiz na variável resposta
par(mfrow=c(1,1))

lm_prec2 <- lm(sqrt(Precisoes)~.,data=bd_reg_prec)
summary(lm_prec2)
X11()
par(mfrow=c(2,2))
plot(lm_prec2) # beleza... só que piorou o qqplot
par(mfrow=c(1,1))

lm_prec3 <- lm(sqrt(Precisoes)~.,data={bd_reg_prec$nm_2 <- bd_reg_prec$nm**2 ; bd_reg_prec[,-4]})
summary(lm_prec3)
X11()
par(mfrow=c(2,2))
plot(lm_prec3) 
par(mfrow=c(1,1))








##############n###################
# Tabelas Nao-Proibitivas

# Precisão
prec_tabela <- read.csv("precisao_tabela.csv",header=TRUE,sep=";")
print(xtable::xtable(prec_tabela,digits=c(2,2,1,4,3,4,4,4,4,4,4)),include.rownames = FALSE)


# Taxa Cobertura
tx_cobert_tabela <- read.csv("taxa_de_cobertura_tabela.csv",header=TRUE,sep=";")
print(xtable::xtable(tx_cobert_tabela,digits=c(2,2,1,4,3,4,4,4,4,4,4)),include.rownames = FALSE)

X11()
par(mfrow=c(2,1))

boxplot(bd_reg_prec$Precisoes~bd_reg_prec$cor, ylab="Precisão", xlab="Correlação")


bd_reg_tx_cober <- read.csv("tx_cobert_regressao.csv",header=TRUE,sep=";")[,-c(4,5)]
boxplot(bd_reg_tx_cober$Taxa_de_cobertura~bd_reg_tx_cober$cor, ylab="Taxa de cobertura", xlab="Correlação")


par(mai=c(0,0,0,0))
plot.new()
legend("center",legend=c("Bootstrap percentil",
                         "Bootstrap Wald",
                         "V: shift na variância",
                         "M: shift na média"),ncol=3,fill=c("blue","red","white","white"))












##################################################################

# Essas tabelas vao ser para o maior e menor conjunto de mn. vao juntar info 
# da amplitude média e taxa de cobertura. E vai ter asteriscos que vão dizer se
# rejeita ou não hipótese de ser .95 a tx de cobertura.

# Os "p"s críticos serão feitos para os níveis de 1%, 5% e 10% de significância.
# E também, para os tamanhos de amostra 10000 e 3000. Dado que os testes de hipótese serão bilaterais, teremos que calcular 3x2x2=12 p-críticos

calc_p_crit <- function(signif, p0=0.95, nm){
  inf <- qnorm(signif/2)*sqrt((p0*(1-p0))/nm)+p0
  sup <- qnorm(1-signif/2)*sqrt((p0*(1-p0))/nm)+p0
  return(c(inf,sup))
}

p_criticos <- t(apply(expand.grid(c(.01,.05,.1),c(10000,3000)),1,function(x){calc_p_crit(signif=x[1],nm=x[2])}))

gabarito_para_asterisco <- data.frame("ICs"=names(prec_tabela)[c(5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10)],
                                      "niveis_signif"=c(rep(c(.01,.05,.1),times=6)),
                                      "inf"=p_criticos[c(rep(c(1,2,3),times=4),rep(c(4,5,6),times=2)),1],
                                      "sup"=p_criticos[c(rep(c(1,2,3),times=4),rep(c(4,5,6),times=2)),2])
  



tabela_top_nm50 <- read.csv("tabela_top_nm50.csv", header=TRUE, sep=";",stringsAsFactors = FALSE)
tabela_top_nm15 <- read.csv("tabela_top_nm15.csv", header=TRUE, sep=";",stringsAsFactors = FALSE)
nm15_side <- xtable::xtable(tabela_top_nm15)
print(nm15_side,include.rownames=FALSE, floating=TRUE, floating.environment="sidewaystable")
nm50_side <- xtable::xtable(tabela_top_nm50)
print(nm50_side,include.rownames=FALSE, floating=TRUE, floating.environment="sidewaystable")
