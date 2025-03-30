rm(list = ls())
library(xtable)
library(foreach)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyr)

setwd("C:/Users/ytid13aw/Documents/Forschung/Projekt 3 Multilevel Imputation/multilevel_imputation_simluation/multilevel_imputation_simluation")

path <- getwd()
path <- paste0(path, "/Data/RepML_neu")

Design <- expand.grid(
  J = c(25, 50),
  N = c(1000),
  misrate = c(0.1, 0.2, 0.3 ,0.4, 0.5),
  model = c("random intercept", "random intercept and random slope"),
  mechanism = c("MAR","MNAR", "MCAR")
) 

sim.results <- vector("list", 60)

for (ii in 1:60) {
  # Liste für die aktuelle "ii"-Element erstellen
  sim.results[[ii]] <- vector("list", 1000)
  
  # Schleife über die 1000 "k"-Elemente
  for (k in 1:1000) {
    print(c(ii,k))
    # Dateiname entsprechend der angegebenen Struktur erstellen
    datei_name <- sprintf("Simulation_Des%d_NumRep%d.rds", ii, k)
    
    # Kompletten Pfad zur Datei erstellen
    datei_pfad <- file.path(path, datei_name)
    
    # Prüfen, ob die Datei existiert
    if (file.exists(datei_pfad)) {
      # Datei einlesen und in die verschachtelte Liste speichern
      sim.results[[ii]][[k]] <- readRDS(datei_pfad)
    } else {
      message(sprintf("Datei fehlt: %s", datei_pfad))
    }
  }
}


sim_counts <- sapply(sim.results, FUN = function(xx) length(xx[sapply(xx, is.list)]))
Design$N <- sim_counts


# Leere Liste zum Speichern der Ergebnisse
significance_counts <- list()

# Signifikanzniveau
alpha <- 0.05

# Anzahl der Koeffizienten, die wir erwarten (hier z.B. 12)
num_coefficients <- 13

# Schleife über jedes Design (1 bis 24)
for (i in seq_along(sim.results)) {
    design_result <- sim.results[[i]]
    
    # Initialisiere eine Liste für jede Methode im aktuellen Design mit einem Zähler für jeden Koeffizienten
    method_counts <- list(
        model = integer(num_coefficients),          # Anzahl der Koeffizienten für die Methode 'model'
        Mice = integer(num_coefficients),
        ranger = integer(num_coefficients),
        levelRanger = integer(num_coefficients),
        ranger5 = integer(num_coefficients),
        levelRanger5 = integer(num_coefficients),
        boost =  integer(num_coefficients),
        boost.dummies =  integer(num_coefficients)
    )
    
    # Schleife über die Simulationsdurchläufe (1 bis 1000)
    for (j in seq_along(design_result)) {
        # Hole die einzelnen Summary-Tabellen
        summaries <- design_result[[j]]
        
        # Überprüfe für jede Methode in summaries die Signifikanz der Koeffizienten
        for (method in names(method_counts)) {
            summary_name <- paste0("Summary.", method)
            
            # Prüfe, ob die Summary-Tabelle existiert
            if (summary_name %in% names(summaries)) {
                # Extrahiere die p-Werte der Methode
                summary_table <- summaries[[summary_name]]
                
                # Prüfe, welche p-Wert-Spalte existiert
                if ("Pr(>|t|)" %in% names(summary_table)) {
                    p_values <- summary_table$`Pr(>|t|)`
                } else if ("p.value" %in% names(summary_table)) {
                    p_values <- summary_table$p.value
                } else {
                    warning(paste("No p-value column found in", summary_name, "for design", i))
                    next  # Springe zur nächsten Methode, falls keine p-Wert-Spalte vorhanden ist
                }
                
                # Überprüfe, ob die Anzahl der p-Werte der erwarteten Anzahl der Koeffizienten entspricht
                if (length(p_values) == num_coefficients) {
                    # Zähle für jeden Koeffizienten separat, ob der p-Wert < alpha ist
                    method_counts[[method]] <- method_counts[[method]] + as.integer(p_values < alpha)
                } else {
                    # Falls die Anzahl der Koeffizienten nicht übereinstimmt, gib eine Warnung aus
                  warning(paste("Unexpected number of coefficients in", summary_name, "for design", i))
                    #q <- q+1
                }
            }
        }
    }
    
    # Speichere die Signifikanzzählungen für das aktuelle Design
    significance_counts[[i]] <- method_counts
}

# Ausgabe der Ergebnisse für jede Methode und jedes Design
significance_proportions <- lapply(seq_along(significance_counts), function(i) {
    # Extrahiere die Anzahl der Durchläufe N für das entsprechende Setting
    N_i <- Design$N[i]
    
    # Für jedes Verfahren (z.B. model, Mice, ranger) teile die Werte durch N_i
    lapply(significance_counts[[i]], function(method_counts) {
        method_counts / N_i
    })
})


zero_coeffs <- c(5, 6, 7, 11, 12, 13)

# Non-zero Koeffizienten (alles außer denjenigen, die 0 sind)
non_zero_coeffs <- setdiff(1:13, zero_coeffs)
method_order <- c("model", "Mice", "ranger", "levelRanger", "ranger5", "levelRanger5", "boost", "boost.dummies")

#################################### unter Nullhypo
plot_data <- do.call(rbind, lapply(seq_along(significance_proportions), function(i) {
    # Hole das aktuelle Setting und die Anteile signifikanter Tests
    setting <- Design[i, ]
    props <- significance_proportions[[i]]
    
    # Erstelle eine Zeile für jede Methode und jeden Koeffizienten
    data.frame(
        J = setting$J,
        misrate = setting$misrate,
        model = setting$model,
        mechanism = setting$mechanism,
        method = factor(rep(names(props), each = length(props[[1]])), levels = method_order),
        coefficient = rep(1:length(props[[1]]), times = length(props)),
        significance_rate = unlist(props)
    )
}))

# Filter für die Koeffizienten, die in Wahrheit 0 sind und entferne "MNAR"
plot_data_zero <- plot_data %>% 
    filter(coefficient %in% zero_coeffs & mechanism != "MNAR" & misrate %in% c(0.1, 0.3, 0.5))

create_model_plot <- function(model_name, mechanism_val) {
    ggplot(plot_data_zero %>% 
               filter(mechanism == mechanism_val & model == model_name), 
           aes(x = factor(coefficient), y = significance_rate, fill = method)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
        facet_grid(J ~ misrate, labeller = labeller(misrate = c(
            `0.1` = "Missing Rate: 0.1", 
            `0.3` = "Missing Rate: 0.3", 
            `0.5` = "Missing Rate: 0.5"
        ))) + 
        labs(
            title = paste("Rejection Rates (mechanism =", mechanism_val, ",", model_name, ")"),
            x = "Coefficient (Zero Coefficients)",
            y = "Rejection Rates",
            fill = NULL
        ) +
        ylim(0, 0.25) +
        scale_x_discrete(labels = c(
            "5" = expression(X[ij5]),
            "6" = expression(X[ij6]),
            "7" = expression(X[ij7]),
            "11" = expression(L[ij1]),
            "12" = expression(L[ij2]),
            "13" = expression(L[ij3])
        )) +  # Neue Beschriftung
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
        scale_fill_discrete(
            labels = c(
                "model" = "model",
                "Mice" = "mice",
                "ranger" = "ranger",
                "levelRanger" = "ranger.dummies",
                "ranger5" = "ranger5",
                "levelRanger5" = "ranger.dummies5",
                "boost" = "boost",
                "boost.dummies" = "boost.dummies"
            )
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.5),
            axis.text.y = element_text(size = 20, face = "bold"), 
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22),
            strip.text = element_text(size = 28),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 26),
            plot.title = element_text(size = 26, face = "bold"),
            legend.position = "bottom",
            panel.spacing = unit(1.5, "lines")
        ) +
        guides(fill = guide_legend(keywidth = 1.2, keyheight = 1.2))
}










# Erstelle für jedes Modell Plots für misrate 0.1 und 0.5

# Für "random intercept"
create_model_plot("random intercept", "MAR")
create_model_plot("random intercept", "MCAR")


# Für "random intercept and random slope"
create_model_plot("random intercept and random slope", "MAR")
create_model_plot("random intercept and random slope", "MCAR")




save_null_model_plots <- function() {
    conditions <- c("random intercept" = "ri", "random intercept and random slope" = "rirs")
    mechanisms <- c("MAR" = "MAR", "MCAR" = "MCAR")  # Großschreibung beibehalten
    
    for (cond in names(conditions)) {
        for (mech in names(mechanisms)) {
            # Erstellen und Speichern des Null Model Plots
            plot <- create_model_plot(cond, mech)
            filename <- paste0("Plot_", conditions[[cond]], "_", mechanisms[[mech]], "_Null.png")
            ggsave(filename, plot, width = 18, height = 7, bg = "white")
        }
    }
}

# Funktion aufrufen
save_null_model_plots()




####################################################################################################
############### coefficient bias ###################################################################
####################################################################################################

# Leere Liste zum Speichern der Ergebnisse
estimate_sum <- list()


# Anzahl der Koeffizienten, die wir erwarten (hier z.B. 12)
num_coefficients <- 13

# Schleife über jedes Design (1 bis 24)
for (i in seq_along(sim.results)) {
    design_result <- sim.results[[i]]
    
    # Initialisiere eine Liste für jede Methode im aktuellen Design mit einem Zähler für jeden Koeffizienten
    method_counts <- list(
        model = integer(num_coefficients),          # Anzahl der Koeffizienten für die Methode 'model'
        Mice = integer(num_coefficients),
        ranger = integer(num_coefficients),
        levelRanger = integer(num_coefficients),
        ranger5 = integer(num_coefficients),
        levelRanger5 = integer(num_coefficients),
        boost =  integer(num_coefficients),
        boost.dummies =  integer(num_coefficients)
    )
    
    # Schleife über die Simulationsdurchläufe (1 bis 1000)
    for (j in seq_along(design_result)) {
        # Hole die einzelnen Summary-Tabellen
        summaries <- design_result[[j]]
        
        # Überprüfe für jede Methode in summaries die Signifikanz der Koeffizienten
        for (method in names(method_counts)) {
            summary_name <- paste0("Summary.", method)
            
            # Prüfe, ob die Summary-Tabelle existiert
            if (summary_name %in% names(summaries)) {
                # Extrahiere die p-Werte der Methode
                summary_table <- summaries[[summary_name]]
                
                # Prüfe, welche p-Wert-Spalte existiert
                if ("Estimate" %in% names(summary_table)) {
                    estimates <- summary_table$Estimate
                } else if ("estimate" %in% names(summary_table)) {
                    estimates <- summary_table$estimate
                } else {
                    warning(paste("No estimatecolumn found in", summary_name, "for design", i))
                    next  # Springe zur nächsten Methode, falls keine p-Wert-Spalte vorhanden ist
                }
                
                # Überprüfe, ob die Anzahl der p-Werte der erwarteten Anzahl der Koeffizienten entspricht
                if (length(estimates) == num_coefficients) {
                    # Zähle für jeden Koeffizienten separat, ob der p-Wert < alpha ist
                    method_counts[[method]] <- method_counts[[method]] + estimates
                } else {
                    # Falls die Anzahl der Koeffizienten nicht übereinstimmt, gib eine Warnung aus
                    warning(paste("Unexpected number of coefficients in", summary_name, "for design", i))
                    #q <- q+1
                }
            }
        }
    }
    
    # Speichere die Signifikanzzählungen für das aktuelle Design
    estimate_sum[[i]] <- method_counts
}

# Ausgabe der Ergebnisse für jede Methode und jedes Design
estimate_proportions <- lapply(seq_along(estimate_sum), function(i) {
    # Extrahiere die Anzahl der Durchläufe N für das entsprechende Setting
    N_i <- Design$N[i]
    
    # Für jedes Verfahren (z.B. model, Mice, ranger) teile die Werte durch N_i
    lapply(estimate_sum[[i]], function(method_counts) {
        method_counts / N_i
    })
})

estimate_proportions 


# Ziel-Vektor zum Abziehen
target_vector <- c(0.3, 0.5, 1, 1.5, 0, 0, 0, 0.5, 1, 1.5, 0, 0, 0)

# Berechnung des Bias für jedes Design und jede Methode
bias_list <- lapply(seq_along(estimate_proportions), function(i) {
    lapply(estimate_proportions[[i]], function(method_proportions) {
        method_proportions - target_vector
    })
})

# Erstellen des `plot_bias`-Datenrahmens
plot_bias <- do.call(rbind, lapply(seq_along(bias_list), function(i) {
    # Extrahiere die Designs
    setting <- Design[i, ]
    biases <- bias_list[[i]]
    
    # Erstelle eine Zeile für jede Methode und jeden Koeffizienten
    data.frame(
        J = setting$J,
        misrate = setting$misrate,
        model = setting$model,
        mechanism = setting$mechanism,
        method = factor(rep(names(biases), each = length(biases[[1]])), levels = method_order),
        coefficient = factor(rep(1:length(biases[[1]]), times = length(biases)), levels = 1:13, labels = c("1", "2", "3", "4", "5", "6", "7", "L1", "L2", "L3", "L4", "L5", "L6")),
        bias = unlist(biases)
    )
}))

# Filter für Null-Koeffizienten (zero_coeffs) und non-zero Koeffizienten (non_zero_coeffs)
plot_bias_zero <- plot_bias %>% filter(as.numeric(coefficient) %in% zero_coeffs & mechanism != "MNAR")
plot_bias_non_zero <- plot_bias %>% filter(as.numeric(coefficient) %in% non_zero_coeffs & mechanism != "MNAR")


# Funktion zur Erstellung von Bias-Plots für Nullhypothese-Koeffizienten mit gedrehter Ansicht
# Funktion zur Erstellung von Bias-Plots für Nullhypothese-Koeffizienten mit gedrehter Ansicht
create_bias_plot <- function(model_name, mechanism_val) {
    ggplot(plot_bias_zero %>% 
               filter(mechanism == mechanism_val & model == model_name & misrate %in% c(0.1, 0.3, 0.5)),  # Filter nur für misrate 0.1, 0.3, 0.5
           aes(x = coefficient, y = bias, fill = method)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
        coord_flip() + # Plot drehen
        facet_grid(J ~ misrate, labeller = labeller(
            J = function(x) paste0(x),  # Zeigt nur die Clustergröße ohne "J:"
            misrate = label_both
        )) +
        labs(
            title = paste("Bias (mechanism =", mechanism_val, ",", model_name, ")"),
            x = " ",
            y = "Bias",
            fill = NULL # Entfernt den Titel der Legende
        ) +
        ylim(-1.2, 0.3) + # Einheitlicher y-Achsenbereich
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"), # Achsentexte fett
            axis.text.y = element_text(size = 20, face = "bold"), # Achsentexte fett
            axis.title.x = element_text(size = 28),
            axis.title.y = element_text(size = 28),
            strip.text = element_text(size = 28), # Größe der Facet-Titel
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 26),
            plot.title = element_text(size = 26, face = "bold"), # Überschrift vergrößert
            legend.position = "bottom", # Legende unterhalb des Plots
            panel.spacing = unit(1.5, "lines") # Erhöht den Abstand zwischen den Facet-Panels
        ) +
        scale_x_discrete( # Da coord_flip() aktiv ist, betrifft dies die Y-Achse
            labels = c(
                "5" = expression(X[ij4]),
                "6" = expression(X[ij5]),
                "7" = expression(X[ij6]),
                "L4" = expression(L[ij4]),
                "L5" = expression(L[ij5]),
                "L6" = expression(L[ij6])
            )
        ) +
        scale_fill_discrete(
            labels = c(
                "model" = "model",
                "Mice" = "mice",
                "ranger" = "ranger",
                "levelRanger" = "ranger.dummies",
                "ranger5" = "ranger5",
                "levelRanger5" = "ranger.dummies5",
                "boost" = "boost",
                "boost.dummies" = "boost.dummies"
            ) 
        ) +
        guides(fill = guide_legend(keywidth = 1.2, keyheight = 1.2)) # Legende anpassen
}






# Funktion zur Erstellung von Bias-Plots für non-zero Koeffizienten mit gedrehter Ansicht
create_non_zero_bias_plot <- function(model_name, mechanism_val, target_vector) {
    # Daten vorbereiten: Relativen Bias berechnen
    plot_data <- plot_bias_non_zero %>%
        filter(mechanism == mechanism_val & model == model_name & misrate %in% c(0.1, 0.3, 0.5)) %>% # Filter nur für misrate 0.1, 0.3, 0.5
        mutate(
            true_coefficient = target_vector[as.numeric(coefficient)],  # Wahre Koeffizienten zuordnen
            relative_bias = ifelse(true_coefficient != 0, bias / true_coefficient, NA) # Relativer Bias
        )
    
    ggplot(plot_data, aes(x = coefficient, y = relative_bias, fill = method)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
        coord_flip() + # Horizontaler Plot
        facet_grid(J ~ misrate, labeller = labeller(
            J = function(x) paste0(x),  # Zeigt nur die Clustergröße ohne "J:"
            misrate = label_both
        )) +
        labs(
            title = paste("Relative Bias (mechanism =", mechanism_val, ",", model_name, ")"),
            x = " ",
            y = "Relative Bias",
            fill = NULL # Entfernt den Titel der Legende
        ) +
        ylim(-1.2, 0.3) + # Einheitlicher y-Achsenbereich
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"), # Achsentexte fett
            axis.text.y = element_text(size = 20, face = "bold"), # Achsentexte fett
            axis.title.x = element_text(size = 28),
            axis.title.y = element_text(size = 28),
            strip.text = element_text(size = 28), # Größe der Facet-Titel
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 26),
            plot.title = element_text(size = 26, face = "bold"), # Überschrift vergrößert
            legend.position = "bottom", # Legende unterhalb des Plots
            panel.spacing = unit(1.5, "lines") # Erhöht den Abstand zwischen den Facet-Panels
        ) +
        scale_x_discrete( # Da coord_flip() aktiv ist, betrifft dies die Y-Achse
            labels = c(
                "1" = "intercept",
                "2" = expression(X[ij1]),
                "3" = expression(X[ij2]),
                "4" = expression(X[ij3]),
                "L1" = expression(L[ij8]),
                "L2" = expression(L[ij9]),
                "L3" = expression(L[ij10])
            )
        ) +
        scale_fill_discrete(
            labels = c(
                "model" = "model",
                "Mice" = "mice",
                "ranger" = "ranger",
                "levelRanger" = "ranger.dummies",
                "ranger5" = "ranger5",
                "levelRanger5" = "ranger.dummies5",
                "boost" = "boost",
                "boost.dummies" = "boost.dummies"
            )
        ) +
        guides(fill = guide_legend(keywidth = 1.2, keyheight = 1.2)) # Legende anpassen
}









# Plots erstellen
create_bias_plot("random intercept", "MAR")
create_bias_plot("random intercept", "MCAR")
create_bias_plot("random intercept and random slope", "MAR")
create_bias_plot("random intercept and random slope", "MCAR")

create_non_zero_bias_plot("random intercept", "MAR", target_vector)
create_non_zero_bias_plot("random intercept", "MCAR", target_vector)
create_non_zero_bias_plot("random intercept and random slope", "MAR", target_vector)
create_non_zero_bias_plot("random intercept and random slope", "MCAR", target_vector)


save_bias_plots <- function(target_vector) {
    conditions <- c("random intercept" = "ri", "random intercept and random slope" = "rirs")
    mechanisms <- c("MAR" = "mar", "MCAR" = "mcar")
    
    for (cond in names(conditions)) {
        for (mech in names(mechanisms)) {
            # Erstellen und Speichern des Bias-Plots
            plot1 <- create_bias_plot(cond, mech)
            filename1 <- paste0("biasplot_", conditions[[cond]], "_", mechanisms[[mech]], "_null.png")
            ggsave(filename1, plot1, width = 18, height = 9, bg = "white")
            
            # Erstellen und Speichern des Non-Zero Bias-Plots
            plot2 <- create_non_zero_bias_plot(cond, mech, target_vector)
            filename2 <- paste0("biasplot_", conditions[[cond]], "_", mechanisms[[mech]], "_H1.png")
            ggsave(filename2, plot2, width = 18, height = 9, bg = "white")
        }
    }
}

# Funktion aufrufen
save_bias_plots(target_vector)












############################################# alterhypo


colors <- c(
    "model" = "#F8766D",        # Rotes Rot
    "Mice" = "#CD9600" ,         # Blau
    "ranger" = "#7CAE00" ,       # Grünes Grün
    "levelRanger" =  "#00BE67" ,  # Orange
    "ranger5" = "#00BFC4",      # Helles Gelb
    "levelRanger5" = "#00A9FF" , # Rosa
    "boost" = "#C77CFF" ,        # Braun
    "boost.dummies" = "#FF61CC"  # Magenta
)
colors_non_zero <- colors[c("model", "Mice", "ranger5", "levelRanger5", "boost", "boost.dummies")]


plot_data <- do.call(rbind, lapply(seq_along(significance_proportions), function(i) {
    # Hole das aktuelle Setting und die Anteile signifikanter Tests
    setting <- Design[i, ]
    props <- significance_proportions[[i]]
    
    # Erstelle eine Zeile für jede Methode und jeden Koeffizienten
    data.frame(
        J = setting$J,
        misrate = setting$misrate,
        model = setting$model,
        mechanism = setting$mechanism,
        method = factor(rep(names(props), each = length(props[[1]])), levels = method_order),
        coefficient = rep(1:length(props[[1]]), times = length(props)),
        significance_rate = unlist(props)
    )
}))

# Filter für die non-zero Koeffizienten und entferne "MNAR"
plot_data_non_zero <- plot_data %>% 
    filter(coefficient %in% non_zero_coeffs & mechanism != "MNAR" & misrate %in% c(0.1, 0.3, 0.5))


# Funktion, um die Plots für jedes Modell und die non-zero Koeffizienten zu erstellen
create_non_zero_model_plot <- function(model_name, mechanism_val) {
    plot_data_filtered <- plot_data_non_zero %>% 
        filter(mechanism == mechanism_val & model == model_name & method %in% names(colors_non_zero))
    
    ggplot(plot_data_filtered, 
           aes(x = factor(coefficient), y = significance_rate, fill = method)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
        facet_grid(J ~ misrate, labeller = labeller(misrate = c(
            `0.1` = "Missing Rate: 0.1", 
            `0.3` = "Missing Rate: 0.3", 
            `0.5` = "Missing Rate: 0.5"
        ))) + # Facet-Titel anpassen
        labs(
            title = paste("Rejection Rates (mechanism =", mechanism_val, ",", model_name, ", non-zero coefficients)"),
            x = "Coefficient (Non-Zero Coefficients)",
            y = "Rejection Rates ",
            fill = NULL # Entfernt den Titel der Legende
        ) +
        ylim(0, 1) +
        scale_fill_manual(values = colors_non_zero) + # Verwende konsistente Farben
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"), # Schrift kleiner, Winkel auf 0
            axis.text.y = element_text(size = 20, face = "bold"), # Achsenzahlen fett
            axis.title.x = element_text(size = 28),
            axis.title.y = element_text(size = 28),
            strip.text = element_text(size = 28),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 26),
            plot.title = element_text(size = 26, face = "bold"), # Überschrift vergrößert
            legend.position = "bottom", # Legende unterhalb des Plots
            panel.spacing = unit(1.5, "lines") # Erhöht den Abstand zwischen den Facet-Panels
        ) +
        scale_x_discrete(
            labels = c(
                "1" = "intercept",
                "2" = expression(X[ij1]),
                "3" = expression(X[ij2]),
                "4" = expression(X[ij3]),
                "8" = expression(L[ij1]),
                "9" = expression(L[ij2]),
                "10" = expression(L[ij3])
            )
        ) + # X-Achsen-Beschriftungen im mathematischen Format
        scale_fill_discrete(
            labels = c(
                "model" = "model",
                "Mice" = "mice",
                "ranger" = "ranger",
                "levelRanger" = "ranger.dummies",
                "ranger5" = "ranger5",
                "levelRanger5" = "ranger.dummies5",
                "boost" = "boost",
                "boost.dummies" = "boost.dummies"
            ) 
        ) +
        guides(fill = guide_legend(keywidth = 1.2, keyheight = 1.2)) # Legende anpassen
}






# Erstelle für jedes Modell Plots für misrate 0.1 und 0.5 (non-zero Koeffizienten)

# Für "random intercept"
create_non_zero_model_plot("random intercept", "MAR")
create_non_zero_model_plot("random intercept", "MCAR")



# Für "random intercept and random slope"
create_non_zero_model_plot("random intercept and random slope", "MAR")
create_non_zero_model_plot("random intercept and random slope", "MCAR")


save_non_zero_model_plots <- function() {
    conditions <- c("random intercept" = "ri", "random intercept and random slope" = "rirs")
    mechanisms <- c("MAR" = "MAR", "MCAR" = "MCAR")  # Großschreibung beibehalten
    
    for (cond in names(conditions)) {
        for (mech in names(mechanisms)) {
            # Erstellen und Speichern des Non-Zero Model Plots
            plot <- create_non_zero_model_plot(cond, mech)
            filename <- paste0("Plot_", conditions[[cond]], "_", mechanisms[[mech]], "_H1.png")
            ggsave(filename, plot, width = 18, height = 8, bg = "white")
        }
    }
}

# Funktion aufrufen
save_non_zero_model_plots()







