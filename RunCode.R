library(doMPI)
# Define the simulation function
source("01_helper.R")
source("02_simulation.R")

args <- commandArgs(trailingOnly = TRUE)
ii <- as.numeric(args[1])

num_replications <- 1000 # 1000 Wiederholungen

# Suppress warnings and messages
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(foreach))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(MASS))))

temp.result <- list()

# Nicht-parallele Schleifen
  for (p in 1:num_replications) {
    
    # Seed ist schon in der Funktion gesetzt, daher nicht erneut setzen
    # seed <- (ii - 1) * num_replications + p
    # set.seed(seed)
    
    # Führe die single_case Funktion aus und fange eventuelle Fehler ab
    s <- try(single_case(ii, p))
    
    if (!class(s) == "try-error") {
      temp.result <- s
    } else {
      temp.result <- warning("This replication failed for some reason")
    }
    # Speichere das Ergebnis für die jeweilige Iteration
    saveRDS(temp.result, paste0("Data/Simulation_Des", ii, "_NumRep", p, ".rds"))
    
  }


# Speichere alle Ergebnisse
saveRDS(temp.result, "SimulationResults.rds")

# Kein Cluster zum Schließen in der nicht-parallelen Version
