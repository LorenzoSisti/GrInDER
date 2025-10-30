################################################################################
################################################################################
#Functions file - GT results
################################################################################
################################################################################

#FUNZIONE PER PLOTTARE PUNTI log(kd)-log(mean) E CALCOLARE CORRELAZIONE

perform_correlation_plot <- function(df_subset_A, df_subset_B, df_subset_AB, descriptor, chain_type_A, chain_type_B, chain_type_AB) {
  
  # Controlla e rimuovi NA e valori <= 0 per la catena A
  df_subset_A <- df_subset_A[df_subset_A$kd > 0 & df_subset_A$mean > 0, ]
  # Controlla e rimuovi NA e valori <= 0 per la catena B
  df_subset_B <- df_subset_B[df_subset_B$kd > 0 & df_subset_B$mean > 0, ]
  # Controlla e rimuovi NA e valori <= 0 per il complesso A-B
  df_subset_AB <- df_subset_AB[df_subset_AB$kd > 0 & df_subset_AB$mean > 0, ]
  
  #Utilizzando due condizionali vado ora a creare due scenari:
  #scenario uno in cui abbiamola combinazione della catena A e della catena B (da ora in poi chiamata "MONOMERO")
  #scenario due in cui abbiamo il complesso A-B (da ora in poi chiamata "DIMERO")
  
  # Combinazione dei dati di kd e mean delle catene A e B per la regressione
  combined_kd <- c(log(df_subset_A$kd), log(df_subset_B$kd))
  combined_mean <- c(log(df_subset_A$mean), log(df_subset_B$mean))
  
  #Porzione della funzione che crea il plot per la combinazione di catene A e B
  if (nrow(df_subset_A) > 0 & nrow(df_subset_B) > 0) {
    cor_value_combined <- cor(combined_kd, combined_mean, use = "complete.obs")
    
    plot(log(df_subset_A$kd), log(df_subset_A$mean), 
         main = paste("Corr. for", descriptor, " (monomer)", "\nCor_monomer:", round(cor_value_combined, 4)), 
         xlab = "log(Kd)", ylab = "log(Mean)", pch = 20, col = "purple")
    
    points(log(df_subset_B$kd), log(df_subset_B$mean), pch = 20, col = "pink")
    #La funz. points() serve per aggiungere punti a un grafico già esistente in R (appunto per me)
    
    # Aggiungiamo la legenda
    legend("topright", legend = c(chain_type_A, chain_type_B), 
           col = c("purple", "pink"), pch = 20)
    
    
    # Aggiungiamo la linea di regressione
    abline(lm(combined_mean ~ combined_kd), col = "firebrick")
    
    
    cat("Correlation for A and B,", descriptor, ":", cor_value_combined, "\n")
  } else {
    warning("No valid data for correlation and plotting for A and B.")
  }
  
  #Setto le coordinate dei quadranti sulla base dei valori discussi in videochiamata
  #Tutto ciò lo faccio solo su closeness centrality e degree in quanto avevamo notato che la distribuzione dei punti di betweenness era alquanto sferica
  
  if (descriptor == "Degree"){
    abline(h = 2.9, lwd=1, col = "black")
    abline(v = -15, lwd=1, col = "black")
  }
  else if (descriptor == "ClosenessCentrality"){
    abline(h = -0.8, lwd=1, col = "black")
    abline(v = -15, lwd=1, col = "black")
  }
  
  
  
  
  # Qui di seguito metto righe di codice per creare il plot per il dimero A-B
  if (nrow(df_subset_AB) > 0) {
    cor_value_AB <- cor(log(df_subset_AB$kd), log(df_subset_AB$mean), use = "complete.obs")
    
    plot(log(df_subset_AB$kd), log(df_subset_AB$mean), 
         main = paste("Correlation for", descriptor, " (dimer)", "\nCor_dimer:", round(cor_value_AB, 4)), 
         xlab = "log(Kd)", ylab = "log(Mean)", pch = 20, col = "goldenrod1")
    
    # Aggiungiamo la linea di regressione
    abline(lm(log(df_subset_AB$mean) ~ log(df_subset_AB$kd)), col = "firebrick")
    
    cat("Correlation for A-B,", descriptor, ":", cor_value_AB, "\n")
  } else {
    warning("No valid data for correlation and plotting for A-B.")
  }
  
  #Setto le coordinate dei quadranti
  if (descriptor == "Degree"){
    abline(h = 2.9, lwd=1, col = "black")
    abline(v = -15, lwd=1, col = "black")
  }
  else if (descriptor == "ClosenessCentrality"){
    abline(h = -0.8, lwd=1, col = "black")
    abline(v = -15, lwd=1, col = "black")
  }
  
}

################################################################################

#FUNZIONE PER SEPARARE I DATI SULLA BASE DI UNA TRESHOLD (PRELUDIO ALLO STUDIO DI DISTIBUZIONE NORMALE)

separate_by_treshold <- function(df, descriptor, treshold) {
  
  #Le righe di sotto separano i punti sopra e sotto la soglia basata sul log(mean)
  
  df_above_treshold <- df[log(df$mean) > treshold, ]
  df_below_treshold <- df[log(df$mean) <= treshold, ]
  return(list(above = df_above_treshold, below = df_below_treshold))
}

################################################################################


replace_invalid_values <- function(vec) {
  n <- length(vec)
  
  for (i in 1:n) {
    # Se il valore è -Inf o 0
    if (vec[i] == -Inf || vec[i] == 0) {
      # Controlla se ci sono valori precedenti e successivi validi
      previous <- if (i > 1) vec[i - 1] else NA
      following <- if (i < n) vec[i + 1] else NA
      
      # Sostituisci con la media dei valori precedente e successivo, se validi
      if (!is.na(previous) && !is.na(following)) {
        vec[i] <- mean(c(previous, following))
      } else if (!is.na(previous)) {
        vec[i] <- previous
      } else if (!is.na(following)) {
        vec[i] <- following
      }
    }
  }
  
  return(vec)
}

################################################################################

#La funzione di sotto dovrebbe fornire una gaussiana 
gaussian <- function(x, mu, sigma, A) {
  A * exp(-(x - mu)^2 / (2 * sigma^2))
}

################################################################################

#Funzione per fitting e plotting delle gaussiane

fit_and_plot_gaussian <- function(log_kd, log_mean, descriptor, threshold, color) {
  # Pulizia dei dati (eliminazione dei -Inf)
  log_kd_cleaned <- replace_invalid_values(log_kd)
  log_mean_cleaned <- replace_invalid_values(log_mean)
  
  #Esegui il test di Shapiro-Wilk su log(Kd) e log(Mean)
  #In realtà mi sto interessando solo alla normalità dei punti log(kd) e non a log(mean)
  #Magari domani lo commento
  shapiro_kd <- shapiro.test(log_kd_cleaned)
  #shapiro_mean <- shapiro.test(log_mean_cleaned)
  
  # Stampa i risultati del test di Shapiro-Wilk
  cat("Shapiro-Wilk test per log(Kd): W =", shapiro_kd$statistic, ", p-value =", shapiro_kd$p.value, "\n")
  #cat("Shapiro-Wilk test per log(Mean): W =", shapiro_mean$statistic, ", p-value =", shapiro_mean$p.value, "\n")
  
  
  #Aggiunta degli istogrammi per visualizzare la distribuzione
  #Lo faccio poiché sono insoddisfatto di alcuni fit gaussiani
  
  hist(log_kd_cleaned, main = paste(threshold, descriptor, "\n log(Kd) distribution"), 
       xlab = "log(Kd)", col = "goldenrod1", border = "black")
  #hist(log_mean_cleaned, main = paste(threshold, descriptor, "\n log(Mean) distribution"), 
  #     xlab = "log(Mean)", col = "lightgreen", border = "black")
  
  # Inizializzazione dei parametri
  init_params <- list(mu = mean(log_kd_cleaned, na.rm = TRUE), 
                      sigma = sd(log_kd_cleaned, na.rm = TRUE), 
                      A = max(log_mean_cleaned, na.rm = TRUE))
  
  # Prova a fare il fitting, cattura eventuali errori
  fit <- tryCatch({
    nls(log_mean_cleaned ~ gaussian(log_kd_cleaned, mu, sigma, A), start = init_params)
  }, error = function(e) {
    warning(paste("Fit fallito per il descrittore:", descriptor, "con threshold:", threshold))
    return(NULL)
  })
  
  # Se il fit ha successo, continua
  if (!is.null(fit)) {
    # Estrai i parametri stimati
    params <- summary(fit)$coefficients
    cat("Descrittore:", descriptor, " (", threshold, ")\n")
    cat("Valore centrale (mu):", params["mu", 1], "\n")
    cat("Deviazione standard (sigma):", params["sigma", 1], "\n")
    cat("Ampiezza (A):", params["A", 1], "\n")
    
    # Genera il Q-Q plot per log_kd_cleaned
    qqnorm(log_kd_cleaned, main = paste("Q-Q Plot for log(Kd) \n", threshold, descriptor))
    qqline(log_kd_cleaned, col = "firebrick")
    
    # Plot dei dati con la gaussiana fittata
    plot(log_kd, log_mean, pch = 20, col = color, 
         main = paste("Fit gaussiano \n (", threshold, ")", descriptor), xlab = "log(Kd)", ylab = "log(Mean)",)
    curve(gaussian(x, params["mu", 1], params["sigma", 1], params["A", 1]), 
          add = TRUE, col = "firebrick", lwd = 2)
    
  }
}

################################################################################
#Funzione per creare il df per fare la ROC
#Separiamo i punti sulla base di una threshold che sta sulle ascisse a "-15", poi plottiamo delle density
################################################################################

vet_roc_fun <- function(df,BaCutoff){
  vet_Kd <- df[,7]
  vet_bin_Kd <- rep(0,length(vet_Kd))
  vet_bin_Kd[log(vet_Kd) < BaCutoff] <- 1
  return(vet_bin_Kd)
}

