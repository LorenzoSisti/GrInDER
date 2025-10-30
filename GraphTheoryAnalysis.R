################################################################################
################################################################################
#libraries
################################################################################
################################################################################

library(bio3d)
library(igraph)
library(networkD3)
library(reshape2)
library(plyr)

################################################################################
################################################################################
#Setting work directories
################################################################################
################################################################################

#It is important to set work directories since if we set them, we could later recall
#the generic work directory (and we can substitute the path just here, at the start of the code)
#Otherwise, we would have to change the path in all the parts of the code that present the code


work_dir <- "/Users/lorenzosisti/Documents/Milanetti_lab_rot/"
InterCont_dir <- paste(work_dir, "InterContactMatrix/", sep="")
GT_output_dir <- paste(work_dir, "GT_output/", sep="")



################################################################################
################################################################################
#Setting functions
################################################################################
################################################################################
#la g di sotto è il nome della variabile
#quando richiamo la funzione nel testo posso anche mettere net_aus
Local_Network_Des_Bin <- function(g){
  BetweennessCentrality <- betweenness(g, v= V(g), directed = FALSE, normalized = TRUE)
  ClosenessCentrality <- closeness(g, v= V(g), normalized = TRUE)
  #Queste due righe di sopra mi ritornano entrambe un vettore con un valore di BC e di CC per ogni residuo cioè per ogni nodo
  #Shortest path <- distances(net_aus, algorithm = "Dijkstra") 
  Degree <- degree(g)
  #ClusteringCoefficient <- transitivity(g, type = "barrat")
  
  df_des <- as.data.frame(cbind(BetweennessCentrality, ClosenessCentrality, Degree)) #così creo un df
  return((df_des))
}

#Main

#Fino ad oggi noi abbiamo calcolato i grafi solamente con una thd di 15 Angstrom
#Ciò è stato fatto con la riga di codice presente sotto che ad ora è commentata
#matL <- list.files(InterCont_dir, pattern = "_15\\.csv$")

#Adesso ci calcoliamo i grafi e poi anche il file .csv per tutte le thd nel vettore delle treshold
#Lo facciamo con la riga di codice di sotto
matL <- list.files(InterCont_dir)


#noi siamo ora in grado di leggere questo file in questa nuova cartella e di crearci la rete

df_NetDes_general <- data.frame() #dataframe di descrittori di rete generale.  

for(i in 1:length(matL)){
  path_aus <- paste(InterCont_dir, matL[i], sep="")
  mat_aus <- as.matrix(read.csv(path_aus, row.names = 1))
  #fai > image(as.matrix(mat_aus)) Vediamo che è la matrice simmetrica di cui prima
  
  #Volevamo solamente un grafo più piccolo da guardare con la p, altrimenti sarebbe venuta un casino.
  #Adesso possiamo anche eliminare questa riga
  #Nel caso in cui volessimo stampare il grafo, dobbiamo ridurlo
  
  #mat_aus <- mat_aus[300:470,300:470]

  
  #NETWORK GENERATION
  
  net_aus <- graph_from_adjacency_matrix(mat_aus)
  #fai plot(net_aus)
  
  Isolated = which(degree(net_aus)==0)
  net_aus_2 = delete.vertices(net_aus, Isolated)
  #plot(net_aus_2, vertex.size=5)
  
  df_pairs <- melt(mat_aus)
  colnames(df_pairs) <- c("Res1", "Res2", "Contacts")
  df_pairs <- df_pairs[df_pairs$Contacts !=0, ]
  p <- simpleNetwork(df_pairs, height = "200px", width = "100px")
  #p
  #La p di sopra è commentata in quanto altrimenti ci vorrebbe tantissimo acommentare il grafo grande ( e illeggibile)
  #Nel caso in cui volessimo stampare il grafo dei contatti intermolecolari, dovremmo imporre la condizione a mat-aus di cui sopra
  
  
  #Graph theory based descriptors. Dal grafo (net_aus_2) ci calcoliamo i descrittori
  #degree.des <- degree(net_aus_2)
  #barplot(degree.des, las=2, cex.names = 0.3)
  
  df_Net <- Local_Network_Des_Bin(net_aus)
  df_Net[df_Net$ClosenessCentrality == "NaN", "ClosenessCentrality"] <- 0
  #Potremmo fare un patto. Potremmo dire che se il degree è nullo (no contatto) eliminiamo la riga
  
  df_Net_BS <- df_Net[df_Net$Degree != 0, ]
  col_aus <- colnames(df_Net_BS)
  
  #la condizione alla fine è per selezionare la terza colonna
  #qui di sotto vado a inventare una nuova colonna che non ho dopo il dollaro
  
  df_Net_BS$Chain <- do.call(rbind, strsplit(rownames(df_Net_BS), "_"))[,3]
  
  head(df_Net_BS)
  
  mean_NetDes <- apply(df_Net_BS[,col_aus], 2, mean)
  
  #Questo di sotto per calcolare la media di A e a
  
  
  df_NetDes <- data.frame()
  for (j in 1:length(col_aus)) {
  
    df_Net_BS_aus <- df_Net_BS[,c(col_aus[j], "Chain")]
  
    colnames(df_Net_BS_aus) <- c("NetDes", "Chain")
    df_NetDes_aus <- ddply(df_Net_BS_aus, .(Chain), summarize, mean = mean(NetDes))
    df_NetDes_aus$NetDes <- col_aus[j]
    
    mean_NetDes_j <- mean_NetDes[names(mean_NetDes) == col_aus[j]]
    TwoChains <- paste(unique(df_NetDes_aus$Chain)[1], unique(df_NetDes_aus$Chain)[2], sep = "-")

    
    vet_aus <- c(TwoChains, mean_NetDes_j, names(mean_NetDes_j)) 
    
    df_NetDes_aus <- rbind(df_NetDes_aus, vet_aus)
    
    df_NetDes <- rbind(df_NetDes, df_NetDes_aus)
  
  }
  
  pdb_aus <- strsplit(matL[i], "_")
  
  
  unlist(pdb_aus)
  unlist(pdb_aus)[1]
  df_NetDes$pdb <- unlist(pdb_aus)[1]
  df_NetDes$DistCutoff <- unlist(strsplit(unlist(pdb_aus)[length(unlist(pdb_aus))], ".csv"))
  
  #la funzione ddply vuole in primo luogo il df su cui operare
  #successivamente vuole i blocchi in cui dividere il df
  
  df_NetDes_general <- rbind(df_NetDes_general, df_NetDes)
  
}

#Fino ad oggi abbiamo salvato i risultati di questo codice usando le righe di sotto
#commenterò adesso queste due righe in modo da poter salvare un nuovo file in cui ci sono i risultati per tutte le thd
#file_name <- "df_NetDes_general"
#write.csv(df_NetDes_general,   paste(GT_output_dir, file_name, ".csv", sep = "")) 

#Adesso proveremo a creare un nuovo file .csv che tenga in considerazione tutte le thd e non solamente 15 Angstrom
file_name <- "df_NetDes_general_allthd"
write.csv(df_NetDes_general,   paste(GT_output_dir, file_name, ".csv", sep = "")) 





#simaRNAweb
#hdock server



  
  
  
  
  