# Installa i pacchetti necessari se non sono già installati
pacchetti <- c("ggplot2", "dplyr", "kableExtra", "webshot")

for (pacchetto in pacchetti) {
  if (!requireNamespace(pacchetto, quietly = TRUE)) {
    install.packages(pacchetto)
  }
}

# Installa il pacchetto webshot
if (!requireNamespace("webshot", quietly = TRUE)) {
  install.packages("webshot")
}

# Specifica la versione di PhantomJS che desideri installare
version <- "2.1.1"

# Installa PhantomJS con la versione specificata
webshot::install_phantomjs(version = version, force = TRUE)

options(scipen = 999)

# Carica i pacchetti
library(ggplot2)
library(dplyr)
library(kableExtra)

getwd()

# Funzione per leggere e formattare i dati
leggi_e_formatta_dati <- function(input_path) {
  dati <- read.csv(input_path, header = FALSE, sep = ";", skip=1, stringsAsFactors = FALSE)
  colnames(dati) <- c("Modality", "OMP", "BlockSize", "tree_creation_time","kernel_search_time", "program_execution_time", "speedup", "efficiency")
  write.csv(dati,file=input_path,row.names = FALSE, sep = ";")
  return(dati)
}

# Funzione per creare e salvare le tabelle
crea_e_salva_tabella <- function(input_path, output_folder) {
  dati <- leggi_e_formatta_dati(input_path)
  
  # Formatta la tabella
  tabella_formattata <- dati %>%
    kable("html") %>%
    kable_styling(full_width = FALSE)
  
  # Salva la tabella in un file HTML nella cartella specificata
  writeLines(as.character(tabella_formattata), file.path(output_folder, "dati.html"))
  
  # Salva la tabella in un file PNG nella cartella specificata
  save_kable(tabella_formattata, file = file.path(output_folder, "dati.png"))
}

# Percorsi delle cartelle e sottocartelle desiderate
cartelle <- c("OMP_CUDA")
sottocartelle <- c("opt0","opt1","opt2","opt3")
num_iterazioni <- c("10000", "250000", "3500000")

# Loop per creare tabelle per tutte le combinazioni di cartelle e sottocartelle
for (cartella in cartelle) {
  for (sottocartella in sottocartelle) {
    for (iterazione in num_iterazioni) {
      file_name <- paste0(iterazione, ".csv")
      path <- file.path("Results", cartella, sottocartella, file_name)
      output_folder <- file.path("Tables", cartella, sottocartella)
      
      # Assicurati che la cartella di output esista, altrimenti creala
      if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
      }
      
      cat("Creazione della tabella per:", path, "\n")
      crea_e_salva_tabella(path, output_folder)
      file.rename(file.path(output_folder, "dati.png"), file.path(output_folder, paste0(iterazione, ".png")))
      file.remove(file.path(output_folder, "dati.html"))
    }
  }
}
