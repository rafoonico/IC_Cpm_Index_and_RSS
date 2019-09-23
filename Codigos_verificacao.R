teste <- list(`Amostragem aleatória simples` = c(1003.24998691387, 1012.78197394474, 
                                                 983.612320663312, 1006.76697217938, 988.466415154886, 995.02600151051, 
                                                 983.081260142846, 1000.40576075835, 1013.44394230378, 1006.26621154823, 
                                                 1010.21041191412, 1002.37180472286, 997.138095577785, 1003.46998899343, 
                                                 1000.45215757801), `Amostragem por conjuntos ordenados` = c(999.185764088234, 
                                                                                                             995.607867129665, 991.81097344686, 1010.21041191412, 990.082597450663, 
                                                                                                             997.040466236032, 1008.36701303413, 989.96049970435, 991.428901431509, 
                                                                                                             1000.45215757801, 998.90254450408, 988.466415154886, 1000.40576075835, 
                                                                                                             992.146423816159, 1012.69087823268))

Cpm <- function(usl,lsl,sigma,mu){
  tao <- mean(c(usl,lsl))
  return((usl-lsl)/(6*sqrt((sigma^2)+(mu-tao)^2)))
}

(Cpm_AAS <- Cpm(1050,950,sd(teste[[1]]),mean(teste[[1]])))
(Cpm_ACO <- Cpm(1050,950,sd(teste[[2]]),mean(teste[[2]])))

delta <- 0.05

n <- length(teste[[1]])

# 1ro - Ok!

(MB_L1 <- Cpm_AAS*sqrt(qchisq(delta/2,df=n)/n))
(MB_S1 <- Cpm_AAS*sqrt(qchisq(1-delta/2,df=n)/n)) ; MB_L1 < Cpm_AAS && MB_S1 > Cpm_AAS 

(MB_L2 <- Cpm_ACO*sqrt(qchisq(delta/2,df=n)/n))
(MB_S2 <- Cpm_ACO*sqrt(qchisq(1-delta/2,df=n)/n)) ; MB_L2 < Cpm_ACO && MB_S2 > Cpm_ACO 

# 2ro - Ok!

d <- (1050-950)/2

t <- (1050+950)/2

sigma_m_chapeu_AAS <- ((d/3)^2)*((sd(teste[[1]])^2)*((mean(teste[[1]])-t)^2)+((sd(teste[[1]])^4)/2))/((sd(teste[[1]])^2)+(mean(teste[[1]])-t)^2)^3
sigma_m_chapeu_ACO <- ((d/3)^2)*((sd(teste[[2]])^2)*((mean(teste[[2]])-t)^2)+((sd(teste[[2]])^4)/2))/((sd(teste[[2]])^2)+(mean(teste[[2]])-t)^2)^3


(CXZ_L1 <- Cpm_AAS-qnorm(1-delta/2)*sigma_m_chapeu_AAS)
(CXZ_S1 <- Cpm_AAS+qnorm(1-delta/2)*sigma_m_chapeu_AAS) ; CXZ_L1 < Cpm_AAS && CXZ_S1 > Cpm_AAS 

(CXZ_L2 <- Cpm_ACO-qnorm(1-delta/2)*sigma_m_chapeu_ACO)
(CXZ_S2 <- Cpm_ACO+qnorm(1-delta/2)*sigma_m_chapeu_ACO) ; CXZ_L2 < Cpm_ACO && CXZ_S2 > Cpm_ACO

# 3ro - Ok!

lambda_AAS <- (n*(mean(teste[[1]])-t)^2)/(sd(teste[[1]])^2)

lambda_ACO <- (n*(mean(teste[[2]])-t)^2)/(sd(teste[[2]])^2)

(ZH_L1 <- Cpm_AAS*sqrt(qchisq(delta/2,n)/(n+lambda_AAS)))
(ZH_S1 <- Cpm_AAS*sqrt(qchisq(1-delta/2,n)/(n+lambda_AAS)))  ; ZH_L1 < Cpm_AAS && ZH_S1 > Cpm_AAS 

(ZH_L2 <- Cpm_ACO*sqrt(qchisq(delta/2,n)/(n+lambda_ACO)))
(ZH_S2 <- Cpm_ACO*sqrt(qchisq(1-delta/2,n)/(n+lambda_ACO)))  ; ZH_L1 < Cpm_ACO && ZH_S1 > Cpm_ACO 

# 4to - Ok!

epslon_AAS <- (mean(teste[[1]])-t)/sd(teste[[1]])
epslon_ACO <- (mean(teste[[2]])-t)/sd(teste[[2]])

v_chapeu_AAS=(n*(1+epslon_AAS^2)^2) / (1+2*epslon_AAS^2) 
v_chapeu_ACO=(n*(1+epslon_ACO^2)^2) / (1+2*epslon_ACO^2)

(BO_L1 <- Cpm_AAS*sqrt(qchisq(delta/2,df=v_chapeu_AAS)/v_chapeu_AAS))
(BO_S1 <- Cpm_AAS*sqrt(qchisq(1-delta/2,df=v_chapeu_AAS)/v_chapeu_AAS))  ; BO_L1 < Cpm_AAS && BO_S1 > Cpm_AAS 

(BO_L2 <- Cpm_ACO*sqrt(qchisq(delta/2,df=v_chapeu_ACO)/v_chapeu_ACO))
(BO_S2 <- Cpm_ACO*sqrt(qchisq(1-delta/2,df=v_chapeu_ACO)/v_chapeu_ACO))  ; BO_L2 < Cpm_ACO && BO_S2 > Cpm_ACO 

# 5to - Ok!

v_chapeu_AAS=(n*(1+epslon_AAS^2)^2)/(1+2*epslon_AAS^2)
v_chapeu_ACO=(n*(1+epslon_ACO^2)^2)/(1+2*epslon_ACO^2)

(LIVRO_L1 <- Cpm_AAS*(qnorm(delta/2)*sqrt(2/(9*v_chapeu_AAS))+1-2/(9*v_chapeu_AAS))^(3/2))
(LIVRO_S1 <- Cpm_AAS*(qnorm(1-delta/2)*sqrt(2/(9*v_chapeu_AAS))+1-2/(9*v_chapeu_AAS))^(3/2))  ; LIVRO_L1 < Cpm_AAS && LIVRO_S1 > Cpm_AAS

(LIVRO_L2 <- Cpm_ACO*(qnorm(delta/2)*sqrt(2/(9*v_chapeu_ACO))+1-2/(9*v_chapeu_ACO))^(3/2))
(LIVRO_S2 <- Cpm_ACO*(qnorm(1-delta/2)*sqrt(2/(9*v_chapeu_ACO))+1-2/(9*v_chapeu_ACO))^(3/2))  ; LIVRO_L2 < Cpm_ACO && LIVRO_S2 > Cpm_ACO

#------------------------------------------------------------------------- COMPARAR COM AS FÓRMULAS

# 1ro - Ok!

# Marcucci and Beazley (1988)

ic.mb <- function(cpm,delta,n){
  cpmL <- cpm*sqrt(qchisq(delta/2,df=n)/n)
  cpmU <- cpm*sqrt(qchisq(1-delta/2,df=n)/n)
  return(c("Inf" = cpmL,"Sup" = cpmU)) # IC exato quando o processo está sob controle. Caso contrario, o IC é conservador.
}

MB1 <- ic.mb(Cpm_AAS,delta,n) ; as.numeric(MB1[1]) == MB_L1 && as.numeric(MB1[2]) == MB_S1
MB2 <- ic.mb(Cpm_ACO,delta,n) ; as.numeric(MB2[1]) == MB_L2 && as.numeric(MB2[2]) == MB_S2

# 2do - Ok

# Chan et.al (1990)

ic.cxz <- function(cpm,alfa,usl,lsl,s,x_bar,t){ # esse alfa é nivel de confiança
  d <- (usl-lsl)/2
  sigma_m_chapeu <- ((d/3)^2)*((s^2)*((x_bar-t)^2)+((s^4)/2))/((s^2)+(x_bar-t)^2)^3
  cpmL <- cpm-qnorm(1-alfa/2)*sigma_m_chapeu
  cpmU <- cpm+qnorm(1-alfa/2)*sigma_m_chapeu
  return(c("Inf" = cpmL,"Sup" = cpmU))
}

CXZ1 <- ic.cxz(Cpm_AAS,delta,1050,950,sd(teste[[1]]),mean(teste[[1]]),t) ; as.numeric(CXZ1[1]) == CXZ_L1 && as.numeric(CXZ1[2]) == CXZ_S1
CXZ2 <- ic.cxz(Cpm_ACO,delta,1050,950,sd(teste[[2]]),mean(teste[[2]]),t) ; as.numeric(CXZ2[1]) == CXZ_L2 && as.numeric(CXZ2[2]) == CXZ_S2

# 3ro - Ok!

# Zimmer and Hubele (1997)

ic.zh <- function(cpm,n,sigma,delta,mu,t){ # lambda é parâmetro de não centralidade
  lambda <- (n*(mu-t)^2)/(sigma^2)
  cpmL <- cpm*sqrt(qchisq(delta/2,n)/(n+lambda))
  cpmU <- cpm*sqrt(qchisq(1-delta/2,n)/(n+lambda))
  return(c("Inf" = cpmL,"Sup" = cpmU))
}

ZH1 <- ic.zh(Cpm_AAS,n,sd(teste[[1]]),delta,mean(teste[[1]]),t) ; as.numeric(ZH1[1]) == ZH_L1 && as.numeric(ZH1[2]) == ZH_S1
ZH2 <- ic.zh(Cpm_ACO,n,sd(teste[[2]]),delta,mean(teste[[2]]),t) ; as.numeric(ZH2[1]) == ZH_L2 && as.numeric(ZH2[2]) == ZH_S2

# 4to - Ok!

# Boyles (1991)

ic.bo <- function(cpm,delta,epslon,n){ 
  v_chapeu=(n*(1+epslon^2)^2) / (1+2*epslon^2) 
  cpmL=cpm*sqrt(qchisq(delta/2,df=v_chapeu)/v_chapeu)
  cpmU=cpm*sqrt(qchisq(1-delta/2,df=v_chapeu)/v_chapeu)
  return(c("Inf" = cpmL,"Sup" = cpmU))
}

BO1 <- ic.bo(Cpm_AAS,delta,epslon_AAS,n) ; as.numeric(BO1[1]) == BO_L1 && as.numeric(BO1[2]) == BO_S1  
BO2 <- ic.bo(Cpm_ACO,delta,epslon_ACO,n) ; as.numeric(BO2[1]) == BO_L2 && as.numeric(BO2[2]) == BO_S2  

# 5to - Ok!

# Sugerido pelo livro

ic <- function(cpm,delta,epslon,n){ 
  v_chapeu=(n*(1+epslon^2)^2)/(1+2*epslon^2)
  cpmL=cpm*(qnorm(delta/2)*sqrt(2/(9*v_chapeu))+1-2/(9*v_chapeu))^(3/2)
  cpmU=cpm*(qnorm(1-delta/2)*sqrt(2/(9*v_chapeu))+1-2/(9*v_chapeu))^(3/2)
  return(c("Inf" = cpmL,"Sup" = cpmU))
}

LIVRO1 <- ic(Cpm_AAS,delta,epslon_AAS,n) ; as.numeric(LIVRO1[1]) == LIVRO_L1 && as.numeric(LIVRO1[2]) == LIVRO_S1  
LIVRO2 <- ic(Cpm_ACO,delta,epslon_ACO,n) ; as.numeric(LIVRO2[1]) == LIVRO_L2 && as.numeric(LIVRO2[2]) == LIVRO_S2 

#------------------------------------------------------------------------- COMPARAR COM O VERDADEIRO PARÂMETRO

CPM <- Cpm(1050,950,10,1000)

MB1[[1]] < CPM && MB1[[2]] > CPM
MB2[[1]] < CPM && MB2[[2]] > CPM

CXZ1[[1]] < CPM && CXZ1[[2]] > CPM
CXZ2[[1]] < CPM && CXZ2[[2]] > CPM

ZH1[[1]] < CPM && ZH1[[2]] > CPM
ZH2[[1]] < CPM && ZH2[[2]] > CPM

BO1[[1]] < CPM && BO1[[2]] > CPM
BO2[[1]] < CPM && BO2[[2]] > CPM

LIVRO1[[1]] < CPM && LIVRO1[[2]] > CPM
LIVRO2[[1]] < CPM && LIVRO2[[2]] > CPM

#------------------------------------------------------------------------- COMPARAR FUNCAO SIMULADORA

epslon <- function(mu,t,sigma){
  return((mu-t)/sigma)
}
  
teste <- list(`Amostragem aleatória simples` = c(1003.24998691387, 1012.78197394474, 
                                                   983.612320663312, 1006.76697217938, 988.466415154886, 995.02600151051, 
                                                   983.081260142846, 1000.40576075835, 1013.44394230378, 1006.26621154823, 
                                                   1010.21041191412, 1002.37180472286, 997.138095577785, 1003.46998899343, 
                                                   1000.45215757801), `Amostragem por conjuntos ordenados` = c(999.185764088234, 
                                                                                                               995.607867129665, 991.81097344686, 1010.21041191412, 990.082597450663, 
                                                                                                               997.040466236032, 1008.36701303413, 989.96049970435, 991.428901431509, 
                                                                                                               1000.45215757801, 998.90254450408, 988.466415154886, 1000.40576075835, 
                                                                                                               992.146423816159, 1012.69087823268))
  
  
  # cor <- parametros[[1]]
  # sdx <- sqrt(parametros[[2]])
  # sdy <- sqrt(parametros[[2]])     # Tome nota: para este experimento, não faz diferença se a dist de Y tiver parâmetros diferentes da dist de X.
  # tamanho <- parametros[[3]]
  # ciclos <- parametros[[4]]
  # mu_x <- parametros[[5]]
  # mu_y <- parametros[[5]]    # Tome nota: para este experimento, não faz diferença se a dist de Y tiver parâmetros diferentes da dist de X.
  # usl <- parametros[[6]]
  # lsl <- parametros[[7]]
  # t <- parametros[[8]]
  # 
  # # Primeiro: gerar a amostra. Tanto a AAS quanto a ACO.
  # amostra <- ACOeAAS(cor,sdx,sdy,tamanho,ciclos,mu_x,mu_y)
  AAS <- teste[[1]]
  ACO <- teste[[2]]
  
  usl <- 1050
  lsl <- 950
  tamanho <- 5
  ciclos <- 3
  t <- 1000
  
  # Segundo: calculando o Cpm estimado.
  cpm_aas <- Cpm(usl,lsl,sd(AAS),mean(AAS)) 
  cpm_aco <- Cpm(usl,lsl,sd(ACO),mean(ACO))
  
  # Terceiro: calculando os IC's.
  Marcucci_e_Beazley_1988 <- c(ic.mb(cpm_aas,.05,tamanho*ciclos),
                               ic.mb(cpm_aco,.05,tamanho*ciclos),
                               ic.mb(cpm_aas,.10,tamanho*ciclos),
                               ic.mb(cpm_aco,.10,tamanho*ciclos))
  Chan_et_al_1990 <- c(ic.cxz(cpm_aas,.05,usl,lsl,sd(AAS),mean(AAS),t),
                       ic.cxz(cpm_aco,.05,usl,lsl,sd(ACO),mean(ACO),t),
                       ic.cxz(cpm_aas,.10,usl,lsl,sd(AAS),mean(AAS),t),
                       ic.cxz(cpm_aco,.10,usl,lsl,sd(ACO),mean(ACO),t))
  Zimmer_e_Hubele_1997 <- c(ic.zh(cpm_aas,tamanho*ciclos,sd(AAS),.05,mean(AAS),t),
                            ic.zh(cpm_aco,tamanho*ciclos,sd(ACO),.05,mean(ACO),t),
                            ic.zh(cpm_aas,tamanho*ciclos,sd(AAS),.10,mean(AAS),t),
                            ic.zh(cpm_aco,tamanho*ciclos,sd(ACO),.10,mean(ACO),t))
  Boyles_1991 <- c(ic.bo(cpm_aas,.05,epslon(mean(AAS),t,sd(AAS)),tamanho*ciclos),
                   ic.bo(cpm_aco,.05,epslon(mean(ACO),t,sd(ACO)),tamanho*ciclos),
                   ic.bo(cpm_aas,.10,epslon(mean(AAS),t,sd(AAS)),tamanho*ciclos),
                   ic.bo(cpm_aco,.10,epslon(mean(ACO),t,sd(ACO)),tamanho*ciclos))
  Sugestao_livro <- c(ic(cpm_aas,.05,epslon(mean(AAS),t,sd(AAS)),tamanho*ciclos),
                      ic(cpm_aco,.05,epslon(mean(ACO),t,sd(ACO)),tamanho*ciclos),
                      ic(cpm_aas,.10,epslon(mean(AAS),t,sd(AAS)),tamanho*ciclos),
                      ic(cpm_aco,.10,epslon(mean(ACO),t,sd(ACO)),tamanho*ciclos))
  
  # Quarto: montando a saída. Escrevendo o título e juntando os resultados.
  Titulo <- c(apply(expand.grid(c("Inf","Sup"),
                                c("AAS","ACO"),
                                c("95%","90%"),
                                c("Marcucci and Beazley (1988)", 
                                  "Chan et.al (1990)", 
                                  "Zimmer and Hubele (1997)", 
                                  "Boyles (1991)", 
                                  "Sugerido pelo livro")), 1, function(x){paste(x, collapse = " - ")}))
  Conteudo <- c(Marcucci_e_Beazley_1988,
                Chan_et_al_1990,
                Zimmer_e_Hubele_1997,
                Boyles_1991,
                Sugestao_livro) #%>% matrix(nrow=1)
  #colnames(Conteudo) <- Titulo
  names(Conteudo) <- Titulo
  # Quinto: retornando o resultado
Conteudo

