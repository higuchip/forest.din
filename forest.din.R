#--------------------------------------------------------------------------------------------------------------------
# forest.din versão beta                                   
# Função para cálculo de taxas demográficas para comunidades de espécies arbóreas                        
# Determinação das taxas de mortalidade, recrutamento, perda e ganho em área basal, mudanças líquidas e rotatividade,
# baseado em:
#
#
#KORNING, J.; BALSLEV, H. Growth and mortality of trees in Amazonian tropical rain forest in Ecuador. Journal of Vegetation Science,
#v.5, n.1, p.77-86, 1994.
#OLIVEIRA FILHO, A. T. et a. Dinâmica da comunidade e populações arbóreas da borda e interior de um remanescente 
#florestal na Serra da Mantiqueira, Minas Gerais, em um intervalo de cinco anos (1999-2004). 
#Revista Brasileira de Botânica, v.30, n.1, p.149-161, 2007.
#SALAMI, B. et al. Influência de variáveis ambientais na dinâmica do componente arbóreo em um fragmento de Floresta
#Ombrófila Mista em Lages, SC. Scientia Forestalis, v.42, n.102, p.197-207, 2014.
#SHEIL, D.; DAVID, BURSLEM, D. F. R. P.; ALDER, D. The interpretation and misinterpretation of mortality rate measures. 
#Journal of Ecology, v.83, n.2, p.331-333, 1995.
#SHEIL, D.; JENNINGS, S.; SAVILL, P. Long-term permanent plot observations of vegetation dynamics in Budongo, a Ugandan 
#rain forest. Journal of Tropical Ecology, v.16, n.6, p.865-882, 2000.
#
#
# Autor:  Pedro Higuchi                                   
# 01/04/2017											                        
#														                              
#														                              
# Observações:											                      
# a) O argumento din (planilha de dados) terá que conter	
# as colunas Parcelas (identificação das parcelas),		    
# Especie (id. espécies), DAP1 (DAP no ano 1) e           
# DAP2 (DAP no ano 2)   
# b) arquivo exemplo de entrada, disponível em https://www.dropbox.com/s/snnco7rkpl22mp2/dados_exemplo.xlsx?dl=0
# c) O argumento t, representa o tempo entre inventários  		     
#----------------------------------------------------------------------------------------------------------------------



forest.din<-function(x,t)
{
  
  x[is.na(x)]<-0
  
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
  tx.ganho.ab.parc<-(1-(1-(sob.ganho.ab.parc+recr.ab.parc/ab.n1.parc$x))^(1/t))*100
  tx.nc.ab.parc<-(((ab.n1.parc$x/ab.n0.parc$x)^(1/t))-1)*100
  turn.ab.parc<-(tx.perda.ab.parc+tx.ganho.ab.parc)/2
  
  din.parc.ab<- cbind(round(ab.n0.parc$x,4), round(sob.ganho.ab.parc,4),
                      round(sob.perda.ab.parc,4),round(mort.ab.parc,4), 
                      round(recr.ab.parc,4),round(ab.n1.parc$x,4),
                      round(tx.perda.ab.parc,4),  round(tx.ganho.ab.parc,4),
                      round(tx.nc.ab.parc,4), round(turn.ab.parc,4))
  colnames(din.parc.ab)<-c("AB0", "G.sob", "P.sob", "AB.m",
                           "AB.r", "AB1", "Tx.perda.AB",
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
  
  tx.perda.ab.spp<-(1-(((ab.n0.spp$x+ABperda.spp)/ab.n0.spp$x)^(1/t)))*100
  tx.ganho.ab.spp<-(1-(1-(sob.ganho.ab.spp+recr.ab.spp/ab.n1.spp$x))^(1/t))*100
  tx.nc.ab.spp<-(((ab.n1.spp$x/ab.n0.spp$x)^(1/t))-1)*100
  turn.ab.spp<-(tx.perda.ab.spp+tx.ganho.ab.spp)/2
  
  din.spp.ab<- cbind(round(ab.n0.spp$x,4), round(sob.ganho.ab.spp,4),
                     round(sob.perda.ab.spp,4),round(mort.ab.spp,4), 
                     round(recr.ab.spp,4),round(ab.n1.spp$x,4),
                     round(tx.perda.ab.spp,4),  round(tx.ganho.ab.spp,4),
                     round(tx.nc.ab.spp,4), round(turn.ab.spp,4))
  colnames(din.spp.ab)<-c("AB0", "G.sob", "P.sob", "AB.m",
                          "AB.r", "AB1", "Tx.perda.AB",
                          "Tx.ganho.AB", "Tx.nc.AB", "Turn.AB")
  
  dinamica<-  list(n.parc=din.parc.ind, 
                   n.spp=din.spp.ind, 
                   ab.parc=din.parc.ab, 
                   ab.spp=din.spp.ab)
  
  n0 <- sum(din.parc.ind[,1])
  n.mort <- sum(din.parc.ind[,3])
  n.recr <- sum(din.parc.ind[,4])
  n1 <- sum(din.parc.ind[,5])
  n0.desv = sd(din.parc.ind[,1])
  n1.desv = sd(din.parc.ind[,5])
  
  ab.ganho.sob<-sum(din.parc.ab[,2])
  ab.perda.sob<- sum(din.parc.ab[,3])
  ab.recr<-sum(din.parc.ab[,5])
  ab.mort<-sum(din.parc.ab[,4])
  ab0<-sum(din.parc.ab[,1])
  ab0.desv = sd(din.parc.ab[,1])
  ab1<-sum(din.parc.ab[,6])
  ab1.desv = sd(din.parc.ab[,6])
  print(dinamica)
  
  
  cat("Abundância ano 1 = ",round(n0,digits=2),"±",round(n0.desv,digits=2),"ind", fill=TRUE)
  cat("Abundância ano 2 = ",round(n1,digits=2),"±",round(n1.desv,digits=2),"ind", fill=TRUE)
  cat("Área basal ano 1 = ",round(ab0,digits=2),"±",round(ab0.desv,digits=2),"m2", fill=TRUE)
  cat("Área basal ano 2 = ",round(ab1,digits=2),"±",round(ab1.desv,digits=2),"m2", fill=TRUE)
  
  write.table(dinamica[[1]], file = "dinamica_n_parcelas.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  write.table(dinamica[[2]], file = "dinamica_n_especies.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  write.table(dinamica[[4]], file = "dinamica_ab_especies.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  write.table(dinamica[[3]], file = "dinamica_ab_parcelas.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  
}

#FIM
#Aplicação
#din<-read.table("dados_exemplo.csv", header = T, sep=";", dec=",")
#forest.din(din,5), onde o primeiro argumento representa o arquivo com os dados de dinâmica e o segundo argumento o tempo entre inventários
