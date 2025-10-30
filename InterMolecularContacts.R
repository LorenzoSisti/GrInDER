################################################################################
################################################################################
#libraries
################################################################################
################################################################################

library(bio3d)

library(igraph)

################################################################################
################################################################################
#Setting work directories
################################################################################
################################################################################

#It is important to set work directories since if we set them, we could later recall
#the generic work directory (and we can substitute the path just here, at the start of the code)
#Otherwise, we would have to change the path in all the parts of the code that present the code


work_dir <- "/Users/lorenzosisti/Documents/Milanetti_lab_rot/"
pdb_dir <- paste(work_dir, "pdb/", sep="")


################################################################################
################################################################################
#Setting variables
################################################################################
################################################################################

name_mat_dir <- "distance_matrix"
name_InterCont_dir <- "InterContactMatrix"
DistCutoff <- 25

vet_DistCutoff <- c(15,20,25,30,35)

list.files(pdb_dir) 

#Questa lista di sopra è un vettore. Per allocarla dobbiamo darle un nome.


################################################################################
################################################################################
#Functions. Here we are defining functions.
################################################################################
################################################################################


#il df a cui ci riferiamo qui di sotto è un pdb qualsiasi 
#ricordiamoci che i file pdb sono dataframe


CA_selection_fn <- function(df){
  df_aus <- df[df$elety == "CA", ]
  return(df_aus) 
}

ContactMatrix_fn <- function(matdist, cutoff){
  matdist[matdist <= cutoff] <- 1
  matdist[matdist > cutoff] <- 0
  return(matdist) 
}


#La funzione CA_selection_fn è una funzione che seleziona Calpha

#La funzione ContactMatrix_fn è una funzione che a partire da una matrice delle 
#distanze e dalla treshold ci fornisce una matrice di contatto binarizzata.

################################################################################
################################################################################
#Main
################################################################################
################################################################################

#pdbList è la stringa che si occupa dicreare una lista dei files che si trovano 
#nella work directory precedentemente definita. Adesso ne abbiamo solo, ma presto ne avremo centinaia 

#Queste righe di codice qui sotto servono essenzialmente per prendere l'eemento i-esimo
#che abbiamo nella cartella (directory) pdb nel nostro computer

pdbL <- list.files(pdb_dir)  
for(i in 1:length(pdbL) ){
  #di sotto c vale è il cutoff
  for (c in 1:length(vet_DistCutoff)) {
    
    DistCutoff <- vet_DistCutoff[c]
  
    path_aus <- paste(pdb_dir, pdbL[i], sep="")
    pdb_aus <- read.pdb(path_aus)
    df_pdb_aus <- pdb_aus$atom 
    df_pdb_aus <- df_pdb_aus[df_pdb_aus$type == "ATOM", ]
    unique(df_pdb_aus$type)
    
    
    
    chain_1 <- unique(df_pdb_aus$chain)[1]
    chain_2 <- unique(df_pdb_aus$chain)[2]
    #Se in PyMOL apriamo 6ikm ci accorgiamo che è composto da moltissime catene. Devo trovare il dimero.
    #Adesso mi serve solamente la catena A e a, perché mi serve solo il df relativo al dimero. Come prendiamo solo le catene del dimero?
    df_dimer_aus <- df_pdb_aus[df_pdb_aus$chain == chain_1 | df_pdb_aus$chain == chain_2, ] 
    #Vedremo che c'è anche un modo per prendere sue catene generiche che non si chiamano A e a, ma lo vedremo un altro giorno.
    df_dimer_CA_aus <- CA_selection_fn(df_dimer_aus)
    #Noi arriveremo a calcolare le matrici delle distanze e quelle dei contatti.
    
    #I pdb di ingresso danno matrici. Le matrici le vogliamo salvare? 
    #Se si devo mettere il file "distance_matrix" nella nostra cartella
  
    
    
    #CA Selection (la matrice delle distanze e dei contatti va calcolata sulla base dei carboni alpha)
    
    
    df_dimer_CA_xyz <- df_dimer_CA_aus[, c("x", "y", "z")]
    #La c di sopra sta per "colonne". Diamo ora un nome ai rownames, perché altrimenti sono brutti.
    row_names <- paste(df_dimer_CA_aus$resid, df_dimer_CA_aus$resno, df_dimer_CA_aus$chain, sep="_")
    rownames(df_dimer_CA_xyz) <- row_names
    
    DistMat <- as.matrix(dist(df_dimer_CA_xyz))
    #as.matrix attempts to turn its argument into a matrix
    #Come salviamo questa matrice in un file? Dobbiamo cambiarle l'estensione.
    
    #Questa riga di sotto serve per cambiare l'estensione da pdb a csv.
    #Prima dobbiamo strippare l'estensione pdb dai file i-esimi presenti nella work directory.
    #Il risultato dello stripping è un tipo di file che è una lista, ma a noi non serve una lista ma un character (che è un vettore?).
    #IMPORTANTE! A che cosa serve utilizzare "unlist"? Lo utilizziamo poiché
    file_name <- unlist( strsplit(pdbL[i], ".pdb"))
    paste(work_dir, name_mat_dir, "/", file_name, ".csv", sep = "")
    #cosa = distmat mentre dove = quello dopo
    
    #saving the distance matrix
    
    write.csv(DistMat,   paste(work_dir, name_mat_dir, "/", file_name, ".csv", sep = "")) 
    #questa matrice è simmetrica
  
    
    #Contact matrix
    
    #la cosa migliore è creare una funzione che crei dei contatti
    
    ContMat <- ContactMatrix_fn(DistMat, DistCutoff) 
    #La matrice dei contatti viene chiamata matrice di adiacenza in Graph theory
    #prova a rifartelo sia con 15 che con 25
    #con 25 possiamo vedere che nel fold più piccolo c'è una periodicità (struttura altamente ripetuta)
    #infatti su pyMOL si vede bene la struttura beta-sheets - loop
    #dalla matrice dei contatti intermolecolari possiamo crearci una rete. 
    #Come si crea una rete? 
    
    image(ContMat)
    
    #INTERMOLECULAR CONTACTS. Questa sotto è una matrice rettangolare
    
    sele_x <- grep(paste("_", chain_1, sep = ""), rownames(DistMat))
    sele_y <- grep(paste("_", chain_2, sep = ""), rownames(DistMat))
    
    #vedo le dimensioni della ContMat
    
    InterContMat <- ContMat[sele_x,sele_y]
    image(InterContMat)
    
    #INTER-MOLECULAR CONTACTS CONSIDERING A SYMMETRIC MATRIX. Questo mi permette di avere matrici quadrate che poi mi servono per i grafi 
    
    InterContMat_Sym <- ContMat
    InterContMat_Sym[sele_x,sele_x] <- 0
    InterContMat_Sym[sele_y,sele_y] <- 0
    #Fai image
    
    #se cambi i valori di treshold devi runnare tutto quello sotto al for senza far correre il for
    
    #Saving the intermolecular contacts
    
    file_name2 <- paste(file_name, "_CutOff_", DistCutoff, sep = "")
    write.csv(InterContMat_Sym,   paste(work_dir, name_InterCont_dir, "/", file_name2, ".csv", sep = "")) 
    
  
  } 
}
