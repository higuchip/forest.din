#============================================================================================================================
# forest.din                            
# Função para determinação da dinâmica de comunidades de espécies arbóreas                        
#---------------------------------------------------------------------------------------------------------------------------- 
#Determinação das taxas de mortalidade, recrutamento, perda e ganho em área basal, mudanças líquidas e rotatividade,
#baseado em:
#
#
#KORNING, J.; BALSLEV, H. Growth and mortality of trees in Amazonian tropical rain forest in Ecuador. Journal of Vegetation Science,
#v.5, n.1, p.77-86, 1994.
#OLIVEIRA FILHO, A. T. et a. Dinâmica da comunidade e populações arbóreas da borda e interior de um remanescente 
#florestal na Serra da Mantiqueira, Minas Gerais, em um intervalo de cinco anos (1999-2004). 
#Revista Brasileira de Botânica, v.30, n.1, p.149-161, 2007.
#SALAMI, B. et al. Influência de variáveis ambientais na dinâmica do componente arbóreo em um fragmento de Floresta
#Ombrófila Mista em Lages, SC. Scientia Forestalis, v.42, n.102, p.197-207, 2014.
#SHEIL, D.; DAVID, BURSLEM, D. F. R. P.; ALDER, D. The interpretation and misinterpretation of mortality rate measures. 
#Journal of Ecology, v.83, n.2, p.331-333, 1995.
#SHEIL, D.; JENNINGS, S.; SAVILL, P. Long-term permanent plot observations of vegetation dynamics in Budongo, a Ugandan 
#rain forest. Journal of Tropical Ecology, v.16, n.6, p.865-882, 2000.
#
#Estimativa de Biomassa Acima do Solo baseado em Chave et al. (2014), considerando fórmula com as seguintes variáveis: 
# D = Diâmetro (cm)
# WD = Densidade da madeira (g.cm-3)
# E = Estimativa de Stress Ambiental, baseado na coordenada geografica (coord)
# Referência:
#CHAVE et al. (2014) Improved allometric models to estimate the aboveground biomass of tropical trees, Global Change Biology, 20 (10), 3177-3190
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#Autor:  Pedro Higuchi                                   
# 01/04/2017							
#Como citar:
#
#HIGUCHI, P. forest.dyn: Função em linguagem de programação estatística R para a determinação da dinâmica de comunidades de espécies arbóreas 2018. DOI: 10.5281/zenodo.1297702 Disponvel em https://github.com/higuchip/forest.din
#
#REJOU-MECHAIN, M.; TANGUY, A.; PIPONIOT, C.; CHAVE, J.; HERAULT, B. 	BIOMASS: Estimating Aboveground Biomass and Its Uncertainty in Tropical Forests. R package version 1.2.
# https://CRAN.R-project.org/package=BIOMASS										                              
#														                              
# Observações:											                      
# a) O argumento x (planilha de dados) terá que conteras colunas Parcela (identificação das parcelas),	Especie (id. espécies), DAP1 (DAP no ano 1) e           
# DAP2 (DAP no ano 2)   
# b) arquivo exemplo de entrada, disponível em https://raw.githubusercontent.com/higuchip/forest.din/master/dados_exemplo.csv
# c) O argumento t, representa o tempo entre inventários  
# d) Argumento coord deve ser do tipo c(long, lat), com valores de graus decimais
# e) Argumento add_wd representa um data.frame com valores de densidade da madeira (g.cm-3) formato com três colunas ("genus", "species", "wd"). Caso argumento add_wd não seja fornecido, a densidade da madeira será estimada com a função getWoodDensity do pacote BIOMASS, baseado em Zanne et al. Global wood density database. Dryad. Identifier: http://datadryad.org/handle/10255/dryad.235 (2009).
# f) Requer pacote BIOMASS
#
# Modificações:
#
# Data: 04/04/2017
# *Add: Determinação da riqueza para os diferentes anos
# 
# Data: 14/05/2017
# *Add: Correcao na determinacao da riqueza para os diferentes anos
# 
# Data: 29/05/2017
# *Add: Correcao na determinacao da taxa de ganho em área basal
#
# Data: 30/05/2017
# *Add: Inclusao de calculo de taxas para a comunidade como um todo
#
#Data: 25/06/2018
#*Add: Inclusao estimativa de biomassa acima do solo
#----------------------------------------------------------------------------------------------------------------------



forest.din<-function(x,t,coord,add_wd = NULL)
{
  
  require(BIOMASS)
  x[is.na(x)]<-0
  
  spp<-x$Especie
  as.character(spp)
  words <- strsplit(as.character(spp), " ")
  genero<-sapply(words, "[", 1)
  especie<-sapply(words, "[", 2)
  
  
  # Obtenção densidade da madeira
  #Chave, Jerome, et al. Towards a worldwide wood economics spectrum. Ecology letters 12.4 (2009): 351-366.
  #Zanne, A. E., et al. Global wood density database. Dryad. Identifier: http://hdl. handle. net/10255/dryad 235 (2009).
  
  
  wd<-getWoodDensity(genus=genero, species = especie, addWoodDensityData = add_wd)
  wd$meanWD
  x$wd<-wd$meanWD
  
  #Determinacao Biomassa Acima do Solo (tonelada)
  #Chave et al. (2014) Improved allometric models to estimate the aboveground biomass of tropical trees, Global Change Biology, 20 (10), 3177-3190
  
  
  #Ano1 
  AGB1 <- computeAGB(x$DAP1, x$wd, coord = coord)
  
  #Ano2 
  AGB2 <- computeAGB(x$DAP2, x$wd, coord = coord)
  
  x$AGB1<-AGB1
  x$AGB2<-AGB2
  

  #calculo área seccional
  x$AS1<-(pi*x$DAP1^2)/40000
  x$AS2<-(pi*x$DAP2^2)/40000
  x$ASdif<-x$AS2-x$AS1
  
  
  
  
  #formula para determinação do número de sobreviventes
  sob<- function(x) {
    sob.p <- nrow(x[x$DAP1>0 & x$DAP2>0,])
    return(sob.p)
  }
  
  #formula para determinação do número de mortas
  mort<- function(x) {
    mort.p <- nrow(x[x$DAP1>0 & x$DAP2==0,])
    return(mort.p)
  }
  
  #formula para determinação do número de recrutas 
  recr<- function(x) {
    recr.p <- nrow(x[x$DAP1==0 & x$DAP2>0,])
    return(recr.p)
  }
  
  #numero de sobreviventes, mortos, recrutas, n0 e n1 por parcela
  n.sob.parc<-as.matrix((by(x, x$Parcela, sob)))
  n.mort.parc<-as.matrix((by(x, x$Parcela, mort)))
  n.recr.parc<-as.matrix((by(x, x$Parcela, recr)))
  n.n0.parc<-n.sob.parc+n.mort.parc
  n.n1.parc<-n.sob.parc+n.recr.parc
  
  #taxas por parcela  
  tx.mort.parc<-(1-(((n.n0.parc-n.mort.parc)/n.n0.parc)^(1/t)))*100
  tx.recr.parc<-(1-(1-(n.recr.parc/n.n1.parc))^(1/t))*100
  tx.nc.parc<-(((n.n1.parc/n.n0.parc)^(1/t))-1)*100
  turn.parc<-(tx.mort.parc+tx.recr.parc)/2
  
  
  #numero de sobreviventes, mortos, recrutas, n0 e n1 por especies
  n.sob.spp<-as.matrix((by(x, x$Especie, sob)))
  n.mort.spp<-as.matrix((by(x, x$Especie, mort)))
  n.recr.spp<-as.matrix((by(x, x$Especie, recr)))
  n.n0.spp<-n.sob.spp+n.mort.spp
  n.n1.spp<-n.sob.spp+n.recr.spp
  
  #taxas por spp
  tx.mort.spp<-(1-(((n.n0.spp-n.mort.spp)/n.n0.spp)^(1/t)))*100
  tx.recr.spp<-(1-(1-(n.recr.spp/n.n1.spp))^(1/t))*100
  tx.nc.spp<-(((n.n1.spp/n.n0.spp)^(1/t))-1)*100
  turn.spp<-(tx.mort.spp+tx.recr.spp)/2
  
  din.parc.ind<-cbind(n.n0.parc, n.sob.parc, n.mort.parc, n.recr.parc, n.n1.parc,
                      round(tx.mort.parc,2),  round(tx.recr.parc,2),
                      round(tx.nc.parc,2), round(turn.parc,2))
  
  din.spp.ind<-cbind(n.n0.spp, n.sob.spp, n.mort.spp, n.recr.spp, n.n1.spp,
                     round(tx.mort.spp,2),  round(tx.recr.spp,2),
                     round(tx.nc.spp,2), round(turn.spp,2))
  colnames(din.parc.ind)<-c("N0", "sob", "mort", "recr", "N1", 
                            "TX.MORT", "TX.RECR", "TX.NC", "TURN")
  colnames(din.spp.ind)<-c("N0", "sob", "mort", "recr", "N1", 
                           "TX.MORT", "TX.RECR", "TX.NC", "TURN")
  
  #AREA BASAL parc
  
  
  sob.ganho.ab.subset<- function(x) {
    sob.subset<-(x[x$DAP1>0 & x$DAP2>0,])
    sob.ganho.ab.subset<-sob.subset[sob.subset$ASdif >0,]
    return(sum(sob.ganho.ab.subset$ASdif))
  }
  sob.ganho.ab.parc<-as.matrix((by(x, x$Parcela, sob.ganho.ab.subset)))
  
  sob.perda.ab.subset<- function(x) {
    sob.subset<-(x[x$DAP1>0 & x$DAP2>0,])
    sob.perda.ab.subset<-sob.subset[sob.subset$ASdif <0,]
    return(sum(sob.perda.ab.subset$ASdif))
  }
  sob.perda.ab.parc<-as.matrix((by(x, x$Parcela, sob.perda.ab.subset)))
  
  
  recr.ab.subset<- function(x) {
    recr.subset<-(x[x$DAP1==0 & x$DAP2>0,])
    return(sum(recr.subset$AS2))
  }
  recr.ab.parc<-as.matrix((by(x, x$Parcela, recr.ab.subset)))
  
  mort.ab.subset<- function(x) {
    mort.subset<-(x[x$DAP1>0 & x$DAP2==0,])
    return(sum(mort.subset$AS1))
  }
  mort.ab.parc<-as.matrix((by(x, x$Parcela, mort.ab.subset)))
  
  ABganho.parc<-sob.ganho.ab.parc+recr.ab.parc
  ABperda.parc<-sob.perda.ab.parc-mort.ab.parc
  
  ab.n0.parc<-aggregate(x$AS1, by=list(Parc=x$Parcela), FUN=sum)
  
  ab.n1.parc<-aggregate(x$AS2, by=list(Parc=x$Parcela), FUN=sum)
  
  tx.perda.ab.parc<-(1-(((ab.n0.parc$x+ABperda.parc)/ab.n0.parc$x)^(1/t)))*100
  tx.ganho.ab.parc<-(1-(1-(ABganho.parc/ab.n1.parc$x))^(1/t))*100
  tx.nc.ab.parc<-(((ab.n1.parc$x/ab.n0.parc$x)^(1/t))-1)*100
  turn.ab.parc<-(tx.perda.ab.parc+tx.ganho.ab.parc)/2
  
  #Biomassa parcela
  biomassa.n0.parc<-aggregate(x$AGB1, by=list(Parc=x$Parcela), FUN=sum)
  biomassa.n1.parc<-aggregate(x$AGB2, by=list(Parc=x$Parcela), FUN=sum)
  
  
  din.parc.ab<- cbind(round(ab.n0.parc$x,4), round(biomassa.n0.parc$x,4), round(sob.ganho.ab.parc,4),
                      round(sob.perda.ab.parc,4),round(mort.ab.parc,4), 
                      round(recr.ab.parc,4),round(ab.n1.parc$x,4),round(biomassa.n1.parc$x,4),
                      round(tx.perda.ab.parc,4),  round(tx.ganho.ab.parc,4),
                      round(tx.nc.ab.parc,4), round(turn.ab.parc,4))
  colnames(din.parc.ab)<-c("AB0", "BAS0", "G.sob", "P.sob", "AB.m",
                           "AB.r", "AB1", "BAS1", "Tx.perda.AB",
                           "Tx.ganho.AB", "Tx.nc.AB", "Turn.AB")
  
  #AREA BASAL SPP
  
  sob.ganho.ab.spp<-as.matrix((by(x, x$Especie, sob.ganho.ab.subset)))
  sob.perda.ab.spp<-as.matrix((by(x, x$Especie, sob.perda.ab.subset)))
  recr.ab.spp<-as.matrix((by(x, x$Especie, recr.ab.subset)))
  mort.ab.spp<-as.matrix((by(x, x$Especie, mort.ab.subset)))
  ABganho.spp<-sob.ganho.ab.spp+recr.ab.spp
  ABperda.spp<-sob.perda.ab.spp-mort.ab.spp
  ab.n0.spp<-aggregate(x$AS1, by=list(Parc=x$Especie), FUN=sum)
  ab.n1.spp<-aggregate(x$AS2, by=list(Parc=x$Especie), FUN=sum)
  
  #Biomassa spp
  
  biomassa.n0.spp<-aggregate(x$AGB1, by=list(Parc=x$Especie), FUN=sum)
  biomassa.n1.spp<-aggregate(x$AGB2, by=list(Parc=x$Especie), FUN=sum)
  
  
  tx.perda.ab.spp<-(1-(((ab.n0.spp$x+ABperda.spp)/ab.n0.spp$x)^(1/t)))*100
  tx.ganho.ab.spp<-(1-(1-(ABganho.spp/ab.n1.spp$x))^(1/t))*100
  tx.nc.ab.spp<-(((ab.n1.spp$x/ab.n0.spp$x)^(1/t))-1)*100
  turn.ab.spp<-(tx.perda.ab.spp+tx.ganho.ab.spp)/2
  
  din.spp.ab<- cbind(round(ab.n0.spp$x,4), round(biomassa.n0.spp$x,4),round(sob.ganho.ab.spp,4),
                     round(sob.perda.ab.spp,4),round(mort.ab.spp,4), 
                     round(recr.ab.spp,4),round(ab.n1.spp$x,4),round(biomassa.n1.spp$x,4),
                     round(tx.perda.ab.spp,4),  round(tx.ganho.ab.spp,4),
                     round(tx.nc.ab.spp,4), round(turn.ab.spp,4))
  colnames(din.spp.ab)<-c("AB0", "BAS0", "G.sob", "P.sob", "AB.m",
                          "AB.r", "AB1","BAS1", "Tx.perda.AB",
                          "Tx.ganho.AB", "Tx.nc.AB", "Turn.AB")
  
  dinamica<-  list(n.parc=din.parc.ind, 
                   n.spp=din.spp.ind, 
                   ab.parc=din.parc.ab, 
                   ab.spp=din.spp.ab)
  
  #Abundância
  n0 <- sum(din.parc.ind[,1])
  n.mort <- sum(din.parc.ind[,3])
  n.recr <- sum(din.parc.ind[,4])
  n1 <- sum(din.parc.ind[,5])
  n0.desv <- sd(din.parc.ind[,1])
  n1.desv <- sd(din.parc.ind[,5])
  
  #Área basal
  ab.ganho.sob<-sum(din.parc.ab[,3])
  ab.perda.sob<- sum(din.parc.ab[,4])
  ab.recr<-sum(din.parc.ab[,6])
  ab.mort<-sum(din.parc.ab[,5])
  ab0<-sum(din.parc.ab[,1])
  ab0.desv <- sd(din.parc.ab[,1])
  ab1<-sum(din.parc.ab[,7])
  ab1.desv <- sd(din.parc.ab[,7])

  #Biomassa
  BAS1<-sum(din.parc.ab[,2])
  BAS1.desv <- sd(din.parc.ab[,2])
  BAS2<-sum(din.parc.ab[,8])
  BAS2.desv <- sd(din.parc.ab[,8])
                          
  
  #Riqueza
  subset.ano1<-x[x$DAP1>0,]
  matriz.spp.ano1<-table(subset.ano1$Parcela,subset.ano1$Especie)
  s.ano1<-ncol(matriz.spp.ano1[,apply(matriz.spp.ano1, 2, sum)>0])
  
  subset.ano2<-x[x$DAP2>0,]
  matriz.spp.ano2<-table(subset.ano2$Parcela,subset.ano2$Especie)
  s.ano2<-ncol(matriz.spp.ano2[,apply(matriz.spp.ano2, 2, sum)>0])
  
  
  # N Total
  dinamica_total.n<-apply(dinamica[[1]], 2,sum)
  tx.mort.total.n<-(1-(((dinamica_total.n[1]-dinamica_total.n[3])/dinamica_total.n[1])^(1/t)))*100
  tx.recr.total.n<-(1-(1-(dinamica_total.n[4]/dinamica_total.n[5]))^(1/t))*100
  tx.nc.total.n<-(((dinamica_total.n[5]/dinamica_total.n[1])^(1/t))-1)*100
  tx.turn.total.n<-(tx.mort.total.n+tx.recr.total.n)/2
  
  
  #AB Total  
  dinamica_total_ab1<-apply(dinamica[[3]], 2,sum)
  ABganho.total<-dinamica_total_ab1[3]+dinamica_total_ab1[6]
  ABperda.total<-dinamica_total_ab1[4]-dinamica_total_ab1[5]
  tx.perda.ab.total<-(1-(((dinamica_total_ab1[1]+ABperda.total)/dinamica_total_ab1[1])^(1/t)))*100
  tx.ganho.ab.total<-(1-(1-(ABganho.total/dinamica_total_ab1[7]))^(1/t))*100
  tx.nc.ab.total<-(((dinamica_total_ab1[7]/dinamica_total_ab1[1])^(1/t))-1)*100
  turn.ab.total<-(tx.perda.ab.total+tx.ganho.ab.total)/2
  
  
  print(dinamica)
  
  cat("DINAMICA DA COMUNIDADE TOTAL", fill=TRUE)
  
  cat("Riqueza ano 1 = ",s.ano1,"especies", fill=TRUE)
  cat("Riqueza ano 2 = ",s.ano2,"especies", fill=TRUE)
  cat("Abundancia ano 1 = ",round(n0,digits=2),"+/-",round(n0.desv,digits=2),"ind", fill=TRUE)
  cat("Abundancia ano 2 = ",round(n1,digits=2),"+/-",round(n1.desv,digits=2),"ind", fill=TRUE)
  cat("Taxa de Mortalidade = ", round(tx.mort.total.n,digits=2), "%.ano-1", fill=TRUE)
  cat("Taxa de Recrutamento = ",round(tx.recr.total.n,digits=2), "%.ano-1", fill=TRUE)
  cat("Taxa de Mudança Líquida em n = ",round( tx.nc.total.n,digits=2), "%.ano-1", fill=TRUE)
  cat("Taxa de Rotatividade Líquida em n = ",round(tx.turn.total.n,digits=2), "%.ano-1", fill=TRUE)
  cat("Area basal ano 1 = ",round(ab0,digits=2),"+/-",round(ab0.desv,digits=2),"m2", fill=TRUE)
  cat("Area basal ano 2 = ",round(ab1,digits=2),"+/-",round(ab1.desv,digits=2),"m2", fill=TRUE)
  cat("Taxa de Perda em AB = ", round(tx.perda.ab.total,digits=2), "%.ano-1", fill=TRUE)
  cat("Taxa de Ganho em AB = ",round(tx.ganho.ab.total,digits=2), "%.ano-1", fill=TRUE)
  cat("Taxa de Mudança Líquida em AB = ",round(tx.nc.ab.total,digits=2), "%.ano-1", fill=TRUE)
  cat("Taxa de Rotatividade Líquida em AB = ",round(turn.ab.total,digits=2), "%.ano-1", fill=TRUE)
  cat("Biomassa ano 1 = ",round(BAS1,digits=2),"+/-",round(BAS1.desv,digits=2),"ton.", fill=TRUE)
  cat("Biomassa ano 2 = ",round(BAS2,digits=2),"+/-",round(BAS2.desv,digits=2),"ton.", fill=TRUE)
    


  write.table(dinamica[[1]], file = "dinamica_n_parcelas.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  write.table(dinamica[[2]], file = "dinamica_n_especies.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  write.table(dinamica[[4]], file = "dinamica_ab_especies.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  write.table(dinamica[[3]], file = "dinamica_ab_parcelas.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  
 
}

#FIM
#Aplicação
#
#dados_exemplo <- read.table("https://raw.githubusercontent.com/higuchip/forest.din/master/dados_exemplo.csv",
#                            header=T, sep = ";", dec=",")
#source("https://raw.githubusercontent.com/higuchip/forest.din/master/forest.din.R")
#forest.din(dados_exemplo, 5, c(-50.17,-27.71)) #onde, 5 representa o tempo entre intervalos e c(long,lat) representa as coordenadas do local

