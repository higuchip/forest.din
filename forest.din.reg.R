#============================================================================================================================
# forest.din.reg                            
# Função para determinação da dinâmica de comunidades regenerantes de espécies arbóreas                        
#---------------------------------------------------------------------------------------------------------------------------- 
#Determinação das taxas de mortalidade, recrutamento, mudanças de classes de tamanho e mudanças líquidas e rotatividade,
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


forest.din.reg <-function(x,t)
{
 
  x[is.na(x)]<-0
   
  
  #formula para determinação do número de sobreviventes
  sob<- function(x) {
    sob.p <- nrow(x[x$C1>0 & x$C2>0,])
    return(sob.p)
  }

  
  #formula para determinação do número de mortas
  mort<- function(x) {
    mort.p <- nrow(x[x$C1>0 & x$C2==0,])
    return(mort.p)
  }
  

  #formula para determinação do número de recrutas 
  recr<- function(x) {
    recr.p <- nrow(x[x$C1==0 & x$C2>0,])
    return(recr.p)
  }
  


  #formula para determinação mudanca para classes superiores 
  up<- function(x) {
    recr.p <- nrow(x[x$C1>0 & x$C2 > x$C1,])
    return(recr.p)
  }
  
  

  #formula para determinação mudanca para classes inferiores 
  
  down<- function(x) {
    recr.p <- nrow(x[x$C2>0 & x$C2 < x$C1,])
    return(recr.p)
  }
  
  #numero de sobreviventes, mortos, recrutas, n0 e n1 por parcela
  n.sob.parc<-as.matrix((by(x, x$Parcela, sob)))
  n.mort.parc<-as.matrix((by(x, x$Parcela, mort)))
  n.recr.parc<-as.matrix((by(x, x$Parcela, recr)))
  n.up.parc<-as.matrix((by(x, x$Parcela, up)))
  n.down.parc<-as.matrix((by(x, x$Parcela, down)))
  n.n0.parc<-n.sob.parc+n.mort.parc
  n.n1.parc<-n.sob.parc+n.recr.parc
  
  #taxas por parcela  
  tx.mort.parc<-(1-(((n.n0.parc-n.mort.parc)/n.n0.parc)^(1/2)))*100
  tx.recr.parc<-(1-(1-(n.recr.parc/n.n1.parc))^(1/2))*100
  tx.nc.parc<-(((n.n1.parc/n.n0.parc)^(1/2))-1)*100
  turn.parc<-(tx.mort.parc+tx.recr.parc)/2
  tx.up.parc<-(1-((n.n0.parc-n.up.parc)/n.n0.parc)^(1/2))*100
  tx.down.parc<-(1-((n.n0.parc-n.down.parc)/n.n0.parc)^(1/2))*100
 
  #numero de sobreviventes, mortos, recrutas, n0 e n1 por especies
  n.sob.spp<-as.matrix((by(x, x$Especie, sob)))
  n.mort.spp<-as.matrix((by(x, x$Especie, mort)))
  n.recr.spp<-as.matrix((by(x, x$Especie, recr)))
  n.up.spp<-as.matrix((by(x, x$Especie, up)))
  n.down.spp<-as.matrix((by(x, x$Especie, down)))
  n.n0.spp<-n.sob.spp+n.mort.spp
  n.n1.spp<-n.sob.spp+n.recr.spp
  
  #taxas por spp
  tx.mort.spp<-(1-(((n.n0.spp-n.mort.spp)/n.n0.spp)^(1/2)))*100
  tx.recr.spp<-(1-(1-(n.recr.spp/n.n1.spp))^(1/2))*100
  tx.nc.spp<-(((n.n1.spp/n.n0.spp)^(1/2))-1)*100
  turn.spp<-(tx.mort.spp+tx.recr.spp)/2
  tx.up.spp<-(1-((n.n0.spp-n.up.spp)/n.n0.spp)^(1/2))*100
  tx.down.spp<-(1-((n.n0.spp-n.down.spp)/n.n0.spp)^(1/2))*100
  
  din.parc.ind<-cbind(n.n0.parc, n.sob.parc, n.mort.parc, n.recr.parc, n.n1.parc,
                      round(tx.mort.parc,2),  round(tx.recr.parc,2), round(tx.up.parc,2),
                      round(tx.down.parc,2),round(tx.nc.parc,2), round(turn.parc,2))
  
  din.spp.ind<-cbind(n.n0.spp, n.sob.spp, n.mort.spp, n.recr.spp, n.n1.spp,
                     round(tx.mort.spp,2),  round(tx.recr.spp,2),round(tx.up.spp,2),
                     round(tx.down.spp,2),round(tx.nc.spp,2), round(turn.spp,2))
  colnames(din.parc.ind)<-c("N0", "sob", "mort", "recr", "N1", 
                            "TX.MORT", "TX.RECR", "TX.UP", "TX.DOWN", "TX.NC", "TURN")
  colnames(din.spp.ind)<-c("N0", "sob", "mort", "recr", "N1", 
                           "TX.MORT", "TX.RECR", "TX.UP", "TX.DOWN","TX.NC", "TURN")
  
  
  dinamica<-  list(n.parc=din.parc.ind, 
                  n.spp=din.spp.ind)
  
  #Abundância
  n0 <- sum(din.parc.ind[,1])
  n.mort <- sum(din.parc.ind[,3])
  n.recr <- sum(din.parc.ind[,4])
  n1 <- sum(din.parc.ind[,5])
  n0.desv <- sd(din.parc.ind[,1])
  n1.desv <- sd(din.parc.ind[,5])
  n.up<- sum(n.up.parc)
  n.down<-sum(n.down.parc)
  
  
  #Riqueza
  subset.ano1<-x[x$C1>0,]
  matriz.spp.ano1<-table(subset.ano1$Parcela,subset.ano1$Especie)
  s.ano1<-ncol(matriz.spp.ano1[,apply(matriz.spp.ano1, 2, sum)>0])
  
  subset.ano2<-x[x$C2>0,]
  matriz.spp.ano2<-table(subset.ano2$Parcela,subset.ano2$Especie)
  s.ano2<-ncol(matriz.spp.ano2[,apply(matriz.spp.ano2, 2, sum)>0])
  
  
  # N Total
  dinamica_total.n<-apply(dinamica[[1]], 2,sum)
  tx.mort.total.n<-(1-(((dinamica_total.n[1]-dinamica_total.n[3])/dinamica_total.n[1])^(1/2)))*100
  tx.recr.total.n<-(1-(1-(dinamica_total.n[4]/dinamica_total.n[5]))^(1/2))*100
  tx.nc.total.n<-(((dinamica_total.n[5]/dinamica_total.n[1])^(1/2))-1)*100
  tx.turn.total.n<-(tx.mort.total.n+tx.recr.total.n)/2
  tx.up.n<-(1-((dinamica_total.n[1]-n.up)/dinamica_total.n[1])^(1/2))*100
  tx.down.n<-(1-((dinamica_total.n[1]-n.down)/dinamica_total.n[1])^(1/2))*100
  
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
  cat("Taxa de Up = ", round(tx.up.n,digits=2), "%.ano-1", fill=TRUE)
  cat("Taxa de Down = ", round(tx.down.n,digits=2), "%.ano-1", fill=TRUE)
  
  
  write.table(dinamica[[1]], file = "dinamica_reg_parcelas.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  write.table(dinamica[[2]], file = "dinamica_reg_especies.csv", row.names = TRUE, dec=",", sep=";", quote=FALSE)
  
}
