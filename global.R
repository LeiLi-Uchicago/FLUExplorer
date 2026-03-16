library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(readr)
library(tidyr)
library(openxlsx)
library(plotly)
library(msaR) 
library(Biostrings)
library(msa)
library(waiter)
library(lubridate)
library(tidyverse)
library(leaflet)             # Fixes: could not find function "leafletOutput"
library(leaflet.minicharts)  # Required for the pie charts on the map
library(shinyWidgets)
library(shinyjs)

# ==========================================
# 1. GLOBAL DATA LOADING & SETUP
# ==========================================
# Version 2: Fixed Neuraminidase (NA) protein loading
RDS_CACHE <- "data/app_cache_flu_v2.rds"

# Subtypes to load
SUBTYPES <- c("H1N1", "H3N2")

if (file.exists(RDS_CACHE)) {
  # ---- FAST PATH: load everything from the pre-built cache ----
  message("Loading data from RDS cache: ", RDS_CACHE)
  cache <- readRDS(RDS_CACHE)
  
  # Check if all required objects are present
  required_objects <- c("metadata_global", "aa_usage_by_clade", "nt_usage_by_clade")
  if (all(required_objects %in% names(cache))) {
    metadata_global    <- cache$metadata_global
    total_raw          <- cache$total_raw
    total_parsed       <- cache$total_parsed
    aa_usage_by_clade  <- cache$aa_usage_by_clade
    aa_usage_by_year   <- cache$aa_usage_by_year
    aa_usage_by_year_month <- cache$aa_usage_by_year_month
    nt_usage_by_clade  <- cache$nt_usage_by_clade
    nt_usage_by_year   <- cache$nt_usage_by_year
    nt_usage_by_year_month <- cache$nt_usage_by_year_month
    important_pos_df   <- cache$important_pos_df
    
    rm(cache)   # free the wrapper list from memory
    message("RDS cache loaded successfully.")
    cache_loaded <- TRUE
  } else {
    message("RDS cache is outdated. Rebuilding...")
    cache_loaded <- FALSE
  }
} else {
  cache_loaded <- FALSE
}

if (!cache_loaded) {
  # ---- SLOW PATH: read from source files and build cache ----
  message("RDS cache not found. Loading from source files and building cache...")
  
  # --- Metadata ---
  all_metadata <- list()
  for (subtype in SUBTYPES) {
    meta_path <- paste0("data/", subtype, "/metadata_merged_annotated.csv")
    if (file.exists(meta_path)) {
      message("Loading metadata for ", subtype)
      # CRITICAL: na = character() ensures 'NA' (Neuraminidase) is NOT treated as a missing value
      meta <- read_csv(meta_path, show_col_types = FALSE, na = character(),
                       col_select = c("Isolate_Id", "Subtype", "Collection_Date", "Location", 
                                      "HA_clade", "NA_clade", "HA_subclade"))
      
      # Rename columns to match app expectations
      # Group = Subtype (A / H1N1)
      meta <- meta %>% 
        dplyr::rename(Group = Subtype, 
                      date = Collection_Date,
                      clade = HA_clade, 
                      G_clade = NA_clade) 
      
      # Standardize Group name to just H1N1 or H3N2 if it's "A / H1N1"
      meta$Group <- gsub("A / ", "", meta$Group)
      
      # Parse Location into region and country
      loc_split <- stringr::str_split(meta$Location, " / ", simplify = TRUE)
      meta$region <- loc_split[, 1]
      meta$country <- loc_split[, 2]
      
      all_metadata[[subtype]] <- meta
    }
  }
  metadata_global <- bind_rows(all_metadata)
  
  # Robust Date Handling
  metadata_global$Year <- stringr::str_extract(metadata_global$date, "^\\d{4}")
  metadata_global$YM <- stringr::str_extract(metadata_global$date, "^\\d{4}-\\d{2}")
  
  total_raw    <- scales::comma(nrow(metadata_global))
  metadata_global <- metadata_global %>% filter(!is.na(Year))
  total_parsed <- scales::comma(nrow(metadata_global))
  
  # --- Usage tables ---
  aa_clade_list <- list()
  aa_year_list  <- list()
  aa_ym_list    <- list()
  nt_clade_list <- list()
  nt_year_list  <- list()
  nt_ym_list    <- list()
  
  for (subtype in SUBTYPES) {
    message("Loading usage tables for ", subtype)
    
    # AA tables
    # CRITICAL: na = character() ensures 'NA' protein is not interpreted as missing
    aa_clade <- read_csv(paste0("data/", subtype, "/AA/aa_usage_by_HA_clade.csv"), show_col_types = FALSE, na = character()) %>%
      dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
      dplyr::rename_with(~ gsub("^HA_clade$", "Clade", .x), any_of("HA_clade")) %>%
      mutate(Group = subtype)
    
    aa_year <- read_csv(paste0("data/", subtype, "/AA/aa_usage_by_Year.csv"), show_col_types = FALSE, na = character()) %>%
      dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
      mutate(Group = subtype)
      
    aa_ym <- read_csv(paste0("data/", subtype, "/AA/aa_usage_by_Year_Month.csv"), show_col_types = FALSE, na = character()) %>%
      dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
      mutate(Group = subtype)
    
    aa_clade_list[[subtype]] <- aa_clade
    aa_year_list[[subtype]]  <- aa_year
    aa_ym_list[[subtype]]    <- aa_ym
    
    # NT tables
    nt_clade <- read_csv(paste0("data/", subtype, "/NT/nt_usage_by_HA_clade.csv"), show_col_types = FALSE, na = character()) %>%
      dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
      dplyr::rename_with(~ gsub("^HA_clade$", "Clade", .x), any_of("HA_clade")) %>%
      dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide")) %>%
      mutate(Group = subtype)
    
    nt_year <- read_csv(paste0("data/", subtype, "/NT/nt_usage_by_Year.csv"), show_col_types = FALSE, na = character()) %>%
      dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
      dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide")) %>%
      mutate(Group = subtype)
      
    nt_ym <- read_csv(paste0("data/", subtype, "/NT/nt_usage_by_Year_Month.csv"), show_col_types = FALSE, na = character()) %>%
      dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
      dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide")) %>%
      mutate(Group = subtype)
      
    nt_clade_list[[subtype]] <- nt_clade
    nt_year_list[[subtype]]  <- nt_year
    nt_ym_list[[subtype]]    <- nt_ym
  }
  
  aa_usage_by_clade      <- bind_rows(aa_clade_list)
  aa_usage_by_year       <- bind_rows(aa_year_list)
  aa_usage_by_year_month <- bind_rows(aa_ym_list)
  
  nt_usage_by_clade      <- bind_rows(nt_clade_list)
  nt_usage_by_year       <- bind_rows(nt_year_list)
  nt_usage_by_year_month <- bind_rows(nt_ym_list)

  # --- Important positions (Placeholder for Flu) ---
  important_pos_df <- data.frame(
    Gene = character(),
    Subtype = character(),
    Position = numeric(),
    Mutation = character(),
    Epitope = character(),
    Clinical_impact = character(),
    Source = character(),
    label = character(),
    stringsAsFactors = FALSE
  )
  
  # --- Save RDS cache ---
  message("Writing RDS cache to: ", RDS_CACHE)
  saveRDS(
    list(
      metadata_global        = metadata_global,
      total_raw              = total_raw,
      total_parsed           = total_parsed,
      aa_usage_by_clade      = aa_usage_by_clade,
      aa_usage_by_year       = aa_usage_by_year,
      aa_usage_by_year_month = aa_usage_by_year_month,
      nt_usage_by_clade      = nt_usage_by_clade,
      nt_usage_by_year       = nt_usage_by_year,
      nt_usage_by_year_month = nt_usage_by_year_month,
      important_pos_df       = important_pos_df
    ),
    file = RDS_CACHE
  )
}

# ---- Post-load steps ----
ALL_AAS <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*","X", "-")

aa_colors <- c(
  "A"="#E41A1C", "C"="#377EB8", "D"="#4DAF4A", "E"="#984EA3",
  "F"="#FF7F00", "G"="#FFFF33", "H"="#A65628", "I"="#F781BF",
  "K"="#999999", "L"="#66C2A5", "M"="#FC8D62", "N"="#8DA0CB",
  "P"="#E78AC3", "Q"="#A6D854", "R"="#FFD92F", "S"="#E5C494",
  "T"="#B3B3B3", "V"="#1B9E77", "W"="#D95F02", "Y"="#7570B3",
  "*"="#000000", "X"="#D3D3D3", "-"="#808080"
)

nt_colors <- c(
  "a"="#E41A1C", "c"="#377EB8", "g"="#4DAF4A", "t"="#984EA3",
  "A"="#E41A1C", "C"="#377EB8", "G"="#4DAF4A", "T"="#984EA3",
  "N"="#000000", "n"="#000000", "-"="#808080"
)

ggmsa_custom_colors <- data.frame(
  names = names(aa_colors),
  color = unname(aa_colors),
  stringsAsFactors = FALSE
)

# Generate rainbow palettes for clades
all_clades <- sort(unique(metadata_global$clade))
clade_colors_vec <- grDevices::rainbow(length(all_clades))
names(clade_colors_vec) <- all_clades

all_g_clades <- sort(unique(metadata_global$G_clade))
g_clade_colors_vec <- grDevices::rainbow(length(all_g_clades))
names(g_clade_colors_vec) <- all_g_clades
