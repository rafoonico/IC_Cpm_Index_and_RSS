load(file="resultado.RData")

load(file="amostras.RData")

Cpm <- function(vetor){ # primeiro, identificar o Cpm de cada cenário
  usl <- 1008
  lsl <- 992
  sigma <- sd(vetor)
  mu <- mean(vetor)
  tao <- mean(c(usl,lsl))
  return((usl-lsl)/(6*sqrt((sigma**2)+(mu-tao)**2)))
}

names(parametros) <- c("corr","sd","n","m","mu","lse","lie","T")

parametros$nm <- parametros$n*parametros$m

rownames(parametros) <- 1:nrow(parametros)

#------------------------------------------------------------------ Precisamos fazer as bases para plotar os histogramas

tamanhos <-  unique(parametros$nm)

correlacoes <- unique(parametros$corr)

bases_para_plot <- vector(mode="list",length = length(tamanhos)*length(correlacoes))

# amostras utilizadas (teste)

cenario <- backup_ams[[as.numeric(rownames(parametros[parametros$corr==1 & parametros$sd==1 & parametros$mu==1003.852 & parametros$nm==50,]))]]

teste_Cpms <- unlist(lapply(cenario,Cpm))

dist_teor <- dnorm(seq(min(teste_Cpms),max(teste_Cpms),by=0.00001),mean=mean(teste_Cpms),sd=sd(teste_Cpms))

# Fazendo a função

# Exemplo para exibir título
# par(mfrow = c(2, 2))
# plot(iris$Petal.Length, iris$Petal.Width)
# plot(iris$Sepal.Length, iris$Petal.Width)
# plot(iris$Sepal.Width, iris$Petal.Width)
# plot(iris$Sepal.Length, iris$Petal.Width)
# mtext("My 'Title' in a strange place", side = 3, line = -21, outer = TRUE)

plot_diag <- function(vetor, titulo, cor, min, max){ # Agora, fazer incluindo grid de plots para diferentes corr
  return({
  par(mar = c(5, 4, 4, 4) + 0.3) 
  hist(vetor, breaks=seq(min(vetor),max(vetor),
                         by=diff(range(vetor))/nclass.Sturges(vetor)),
       xlab="Cpm",ylab="Freq.",main=titulo,xlim=c(min,max),xaxt="n")
  par(new=TRUE)
  dist_teor <- dnorm(seq(min(vetor),max(vetor),by=0.00001),mean=mean(vetor),sd=sd(vetor))
  plot(dist_teor~seq(min(vetor),max(vetor),by=0.00001),
       type="l",col=cor,xlab="",ylab = "",yaxt='n',xlim=c(min,max))
  axis(side=4)
  })
}

frame_gif <- function(nm,sd,media){
  library(magrittr)
  
  # corr=0 
  corr0 <- backup_ams[[as.numeric(rownames(parametros[parametros$corr==0 & 
                                           parametros$sd==sd & 
                                           parametros$mu==media & 
                                           parametros$nm==nm,]))]] %>% lapply(Cpm) %>% unlist()
  # corr=0.5 
  corr0.5 <- backup_ams[[as.numeric(rownames(parametros[parametros$corr==0.5 & 
                                                        parametros$sd==sd & 
                                                        parametros$mu==media & 
                                                        parametros$nm==nm,]))]] %>% lapply(Cpm) %>% unlist()
  # corr=0.8 
  corr0.8 <- backup_ams[[as.numeric(rownames(parametros[parametros$corr==0.8 & 
                                                        parametros$sd==sd & 
                                                        parametros$mu==media & 
                                                        parametros$nm==nm,]))]] %>% lapply(Cpm) %>% unlist()
  # corr=1 
  corr1<- backup_ams[[as.numeric(rownames(parametros[parametros$corr==1 & 
                                                     parametros$sd==sd & 
                                                     parametros$mu==media & 
                                                     parametros$nm==nm,]))]] %>% lapply(Cpm) %>% unlist()
  
  min <- min(c(corr0,corr0.5,corr0.8,corr1))
  max <- max(c(corr0,corr0.5,corr0.8,corr1))
  
  return({
    
    par(mfrow=c(2,2))
    plot_diag(corr0,expression(paste(rho," = 0",sep="")),cor="red",min,max)
    plot_diag(corr0.5,expression(paste(rho," = 0.5",sep="")),cor="blue",min,max)
    plot_diag(corr0.8,expression(paste(rho," = 0.8",sep="")),cor="green",min,max)
    plot_diag(corr1,expression(paste(rho," = 1",sep="")),cor="orange",min,max)
    
    title(bquote(mu~" = "~.(media)~", "~sigma~" = "~.(sd)~", e nm ="~.(nm)),
          line=-1.5, 
          outer=TRUE, cex.main=3)
    
    list("eqm0" =mean((2-corr0)^2),
    "eqm0.5" = mean((2-corr0.5)^2),
    "eqm0.8" = mean((2-corr0.8)^2),
    "eqm1" = mean((2-corr1)^2))
    
  })
}




frame_gif(nm=15, sd=1.777778, media=1000)
frame_gif(nm=25, sd=1.777778, media=1000)
frame_gif(nm=30, sd=1.777778, media=1000)
frame_gif(nm=50, sd=1.777778, media=1000)
frame_gif(nm=15, sd=4.020075, media=1000)
frame_gif(nm=25, sd=4.020075, media=1000)
frame_gif(nm=30, sd=4.020075, media=1000)
frame_gif(nm=50, sd=4.020075, media=1000)
frame_gif(nm=15, sd=15.841192, media=1000)
frame_gif(nm=25, sd=15.841192, media=1000)
frame_gif(nm=30, sd=15.841192, media=1000)
frame_gif(nm=50, sd=15.841192, media=1000)
frame_gif(nm=15, sd=1, media=1000.882)
frame_gif(nm=25, sd=1, media=1000.882)
frame_gif(nm=30, sd=1, media=1000.882)
frame_gif(nm=50, sd=1, media=1000.882)
frame_gif(nm=15, sd=1, media=1001.738)
frame_gif(nm=25, sd=1, media=1001.738)
frame_gif(nm=30, sd=1, media=1001.738)
frame_gif(nm=50, sd=1, media=1001.738)
frame_gif(nm=15, sd=1, media=1003.852)
frame_gif(nm=25, sd=1, media=1003.852)
frame_gif(nm=30, sd=1, media=1003.852)
frame_gif(nm=50, sd=1, media=1003.852)







