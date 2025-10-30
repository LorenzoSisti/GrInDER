################################################################################
#libraries
################################################################################

library(bio3d)
library(igraph)
library(networkD3)
library(reshape2)
library(plyr)
library(MASS)
library(pROC)

################################################################################
#Setting work directories
################################################################################

work_dir <- "/Users/lorenzosisti/Documents/Milanetti_lab_rot/"
files_dir <- paste(work_dir, "files/", sep="")
GT_output_dir <- paste(work_dir, "GT_output/", sep="")
script_dir <- paste(work_dir, "script/", sep = "")
function_script <- paste(script_dir, "Function_file.R", sep = "")

#Calling the function file

source(function_script)

#Se a questo pezzo di codice scritto sopra attacco il main di GT_results_ottimizzato.R il codice funziona come GT_results_ottimizzato.R originale
#Procedo a estendere quello pensato la scorsa volta in IIT insieme

################################################################################
#Main
################################################################################

df_kd <- read.csv(paste(files_dir, "Dias_complexes.csv", sep = ""), sep = ";")
df_GT <- read.csv(paste(GT_output_dir, "df_NetDes_general.csv", sep = ""), sep = ",")

df_GT$kd <- 0 #Con la riga di codice sopra creiamo df_GT$kd e la inizializziamo inizialmente a zero

#Il seguente for ha lo scopo di popolare la colonna kd di df_GT 
for(i in 1:nrow(df_GT)){
  pdb_aus <- df_GT[i,"pdb"] 
  
  kd_aus <- df_kd[df_kd$PDB_ID == pdb_aus, "Value.M."] #ma questo potrebbe anche essere zero
  
  if(length(kd_aus) != 0){
    df_GT[i,"kd"] <- kd_aus
  }
}

file_name <- "df_NetDes_general_kd_ottimizzato"
write.csv(df_GT,   paste(GT_output_dir, file_name, ".csv", sep = "")) 

################################################################################
#GT_results
################################################################################

#ESTRAGGO I DESCRITTORI DI RETE UNICI
#df$NetDes estrae la colonna "NetDes" dal dataframe.unique(): Rimuove i duplicati, quindi ottieni solo i tre descrittori unici.


vet_descrittori <- unique(df_GT$NetDes) 

#Sezione fancy color
################################################################################
#Per il monomero
col_trasp <- rgb(red = 0.5, green = 0, blue = 0.5, alpha = 0.5)
#Per il dimero
col_trasp_dim <- rgb(red = 0.7, green = 0.3, blue = 0, alpha = 0.5)
################################################################################
  
for (j in vet_descrittori) {
  
  #Creiamo e inizializziamo i dataframe per A, B e A-B
  df_subset_A <- df_GT[df_GT$NetDes == j & df_GT$Chain == "A", ]
  df_subset_B <- df_GT[df_GT$NetDes == j & df_GT$Chain == "B", ]
  df_subset_AB <- df_GT[df_GT$NetDes == j & df_GT$Chain == "A-B", ]
  
  #Chiamiamo ora la funzione di correlazione e plot definita nel blocco delle funzioni
  perform_correlation_plot(df_subset_A, df_subset_B, df_subset_AB, j, "Chain A", "Chain B", "Chain A-B")
  
  
}

################################################################################
#ROC
################################################################################

#Creiamo e inizializziamo i dataframe per A, B e A-B per ogni descrittore di nostro interesse (ottimizzeremo in seguito)

df_subset_A_degree <- df_GT[df_GT$NetDes == "Degree" & df_GT$Chain == "A", ]
df_subset_B_degree <- df_GT[df_GT$NetDes == "Degree" & df_GT$Chain == "B", ]

#Occhio alla parte in cui sopprimo gli assi con xaxt e yaxt, vedi bene se puoi migliorarla.
df_subset_monomer_degree <- rbind(df_subset_A_degree, df_subset_B_degree)
vet_bin_monomer_degree <- vet_roc_fun(df_subset_monomer_degree, -15)
df_roc_monomer_degree <- as.data.frame(cbind(df_subset_monomer_degree$mean, vet_bin_monomer_degree))
colnames(df_roc_monomer_degree) <- c("value","status")
plot(density(df_roc_monomer_degree$value[df_roc_monomer_degree$status == 1]), lwd=2, col="red", ylim=c(0,0.15), xlim=c(0,45))
par(new=T)
plot(density(df_roc_monomer_degree$value[df_roc_monomer_degree$status == 0]), lwd=2, col="blue", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0,0.15), xlim=c(0,45))
roc_obj_monomer_degree <- roc(df_roc_monomer_degree$status, df_roc_monomer_degree$value)
auc(roc_obj_monomer_degree)
plot(roc_obj_monomer_degree, main = paste("ROC_monomer_Degree", "AUC = ", round(auc(roc_obj_monomer_degree), digits = 4)))

#creazione e inizializzazione df dimer
df_subset_dimer_degree <- df_GT[df_GT$NetDes == "Degree" & df_GT$Chain == "A-B", ]

vet_bin_dimer_degree <- vet_roc_fun(df_subset_dimer_degree, -15)
df_roc_dimer_degree <- as.data.frame(cbind(df_subset_dimer_degree$mean, vet_bin_dimer_degree))
colnames(df_roc_dimer_degree) <- c("value","status")
plot(density(df_roc_dimer_degree$value[df_roc_dimer_degree$status == 1]), lwd=2, col="red", ylim=c(0,0.15), xlim=c(0,45))
par(new=T)
plot(density(df_roc_dimer_degree$value[df_roc_dimer_degree$status == 0]), lwd=2, col="blue", ylim=c(0,0.15), xlim=c(0,45))
roc_obj_dimer_degree <- roc(df_roc_dimer_degree$status, df_roc_dimer_degree$value)
auc(roc_obj_dimer_degree)
plot(roc_obj_dimer_degree, main = paste("ROC_dimer_Degree", "AUC = ", round(auc(roc_obj_dimer_degree), digits = 4)))

df_subset_A_clos <- df_GT[df_GT$NetDes == "ClosenessCentrality" & df_GT$Chain == "A", ]
df_subset_B_clos <- df_GT[df_GT$NetDes == "ClosenessCentrality" & df_GT$Chain == "B", ]

df_subset_monomer_clos <- rbind(df_subset_A_clos, df_subset_B_clos)
vet_bin_monomer_clos <- vet_roc_fun(df_subset_monomer_clos, -15)
df_roc_monomer_clos <- as.data.frame(cbind(df_subset_monomer_clos$mean, vet_bin_monomer_clos))
colnames(df_roc_monomer_clos) <- c("value","status")
plot(density(df_roc_monomer_clos$value[df_roc_monomer_clos$status == 1]), lwd=2, col="red", ylim=c(0,10), xlim=c(0,01))
par(new=T)
plot(density(df_roc_monomer_clos$value[df_roc_monomer_clos$status == 0]), lwd=2, col="blue", ylim=c(0,10), xlim=c(0,01))
roc_obj_monomer_clos <- roc(df_roc_monomer_clos$status, df_roc_monomer_clos$value)
auc(roc_obj_monomer_clos)
plot(roc_obj_monomer_clos, main = paste("ROC_monomer_Clos", "AUC = ", round(auc(roc_obj_monomer_clos), digits = 4)))

#creazione e inizializzazione df dimer
df_subset_dimer_clos <- df_GT[df_GT$NetDes == "ClosenessCentrality" & df_GT$Chain == "A-B", ]

vet_bin_dimer_clos <- vet_roc_fun(df_subset_dimer_clos, -15)
df_roc_dimer_clos <- as.data.frame(cbind(df_subset_dimer_clos$mean, vet_bin_dimer_clos))
colnames(df_roc_dimer_clos) <- c("value","status")
plot(density(df_roc_dimer_clos$value[df_roc_dimer_clos$status == 1]), lwd=2, col="red", ylim=c(0,10), xlim=c(0,01))
par(new=T)
plot(density(df_roc_dimer_clos$value[df_roc_dimer_clos$status == 0]), lwd=2, col="blue", ylim=c(0,10), xlim=c(0,01))
roc_obj_dimer_clos <- roc(df_roc_dimer_clos$status, df_roc_dimer_clos$value)
auc(roc_obj_dimer_clos)
plot(roc_obj_dimer_clos, main = paste("ROC_dimer_Clos", "AUC = ", round(auc(roc_obj_dimer_clos), digits = 4)))


