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
# library(leaflet)             # Fixes: could not find function "leafletOutput"
# library(leaflet.minicharts)  # Required for the pie charts on the map
library(shinyWidgets)
library(shinyjs)

# Disable scientific notation for the session
options(scipen = 999)

# ==========================================
# 1. GLOBAL DATA LOADING & SETUP
# ==========================================
# Version 5: More clade columns in summary
RDS_CACHE <- "data/app_cache_flu.rds"

# Subtypes to load
SUBTYPES <- c("H1N1", "H3N2")

if (file.exists(RDS_CACHE)) {
  # ---- FAST PATH: load everything from the pre-built cache ----
  message("Loading data from RDS cache: ", RDS_CACHE)
  cache <- readRDS(RDS_CACHE)
  
  # Check if all required objects are present
  required_objects <- c("metadata_global", "aa_usage_by_clade", "nt_usage_by_clade", "metadata_summary_stats", "aa_usage_by_group")
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
    important_pos_df       = cache$important_pos_df

    # Pre-calculated groups
    aa_usage_by_group      <- cache$aa_usage_by_group
    nt_usage_by_group      <- cache$nt_usage_by_group

    # Pre-calculated stats

    metadata_summary_stats <- cache$metadata_summary_stats
    total_countries_val    <- cache$total_countries_val
    time_range_val         <- cache$time_range_val
    metadata_groups        <- cache$metadata_groups
    metadata_years         <- cache$metadata_years
    
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
                                      "HA_clade", "NA_clade", "HA_subclade", 
                                      "HA_proposedSubclade", "HA_short_clade", "HA_legacy_clade"))
      
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
  
  # Clean up empty clade names
  clade_cols_to_clean <- c("clade", "G_clade", "HA_subclade", "HA_proposedSubclade", "HA_short_clade", "HA_legacy_clade")
  for(col in clade_cols_to_clean) {
    if(col %in% colnames(metadata_global)) {
      metadata_global[[col]] <- ifelse(is.na(metadata_global[[col]]) | metadata_global[[col]] == "" | metadata_global[[col]] == "trace 0", "Unknown", metadata_global[[col]])
    }
  }

  # Robust Date Handling
  metadata_global$Year <- stringr::str_extract(metadata_global$date, "^\\d{4}")
  metadata_global$YM <- stringr::str_extract(metadata_global$date, "^\\d{4}-\\d{2}")
  
  total_raw    <- scales::comma(nrow(metadata_global))
  metadata_global <- metadata_global %>% filter(!is.na(Year))
  total_parsed <- scales::comma(nrow(metadata_global))
  
  # --- Usage tables ---
  aa_temp <- list()
  nt_temp <- list()
  
  for (subtype in SUBTYPES) {
    message("Loading usage tables for ", subtype)
    
    # AA tables
    aa_dir <- paste0("data/", subtype, "/AA/")
    if (dir.exists(aa_dir)) {
      aa_files <- list.files(aa_dir, pattern = "aa_usage_by_.*\\.csv", full.names = TRUE)
      for (f in aa_files) {
        # Extract group name: aa_usage_by_XYZ.csv -> XYZ
        group_name <- sub(".*_by_(.*)\\.csv$", "\\1", basename(f))
        
        df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
          dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
          mutate(Group = subtype)
          
        aa_temp[[group_name]][[subtype]] <- df
      }
    }
    
    # NT tables
    nt_dir <- paste0("data/", subtype, "/NT/")
    if (dir.exists(nt_dir)) {
      nt_files <- list.files(nt_dir, pattern = "nt_usage_by_.*\\.csv", full.names = TRUE)
      for (f in nt_files) {
        group_name <- sub(".*_by_(.*)\\.csv$", "\\1", basename(f))
        
        df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
          dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
          dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide")) %>%
          mutate(Group = subtype)
          
        nt_temp[[group_name]][[subtype]] <- df
      }
    }
  }
  
  # Bind rows for each group
  aa_usage_by_group <- list()
  for (gn in names(aa_temp)) {
    aa_usage_by_group[[gn]] <- bind_rows(aa_temp[[gn]])
  }
  
  nt_usage_by_group <- list()
  for (gn in names(nt_temp)) {
    nt_usage_by_group[[gn]] <- bind_rows(nt_temp[[gn]])
  }
  
  # Map to existing top-level variables for backward compatibility with server.R
  # Note: We rename HA_clade to Clade ONLY for these convenience variables
  aa_usage_by_clade      <- aa_usage_by_group[["HA_clade"]] %>% dplyr::rename_with(~ gsub("^HA_clade$", "Clade", .x), any_of("HA_clade"))
  aa_usage_by_year       <- aa_usage_by_group[["Year"]]
  aa_usage_by_year_month <- aa_usage_by_group[["Year_Month"]]
  
  nt_usage_by_clade      <- nt_usage_by_group[["HA_clade"]] %>% dplyr::rename_with(~ gsub("^HA_clade$", "Clade", .x), any_of("HA_clade"))
  nt_usage_by_year       <- nt_usage_by_group[["Year"]]
  nt_usage_by_year_month <- nt_usage_by_group[["Year_Month"]]

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
  
  # Pre-calculate these before saving if they weren't loaded from cache
  message("Pre-aggregating metadata statistics...")
  metadata_summary_stats <- metadata_global %>%
    group_by(Year, Group, region, country, clade, G_clade, 
             HA_subclade, HA_proposedSubclade, HA_short_clade, HA_legacy_clade) %>%
    summarise(n = n(), .groups = "drop")
  
  total_countries_val <- length(unique(metadata_global$country))
  time_range_val <- paste(min(metadata_global$Year, na.rm=T), "-", max(metadata_global$Year, na.rm=T))
  metadata_groups <- sort(na.omit(unique(metadata_global$Group)))
  metadata_years <- sort(na.omit(unique(metadata_global$Year)), decreasing = TRUE)

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
      important_pos_df       = important_pos_df,
      aa_usage_by_group      = aa_usage_by_group,
      nt_usage_by_group      = nt_usage_by_group,
      metadata_summary_stats = metadata_summary_stats,
      total_countries_val    = total_countries_val,
      time_range_val         = time_range_val,
      metadata_groups        = metadata_groups,
      metadata_years         = metadata_years
    ),
    file = RDS_CACHE
  )
}

# ==========================================
# 2. COORDINATE LOOKUP DATA (Pre-calculated for Performance)
# ==========================================
# region_coords <- data.frame(
#   region = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania"),
#   lat = c(1.0, 34.0, 48.0, 45.0, -15.0, -25.0),
#   lng = c(17.0, 100.0, 10.0, -100.0, -60.0, 135.0),
#   stringsAsFactors = FALSE
# )
# 
# # Move world_coords out of server.R to global.R to avoid recalculation on every session
# # ggplot2::map_data is slow, so we do it once here.
# message("Pre-calculating world coordinates...")
# world_coords <- ggplot2::map_data("world") %>%
#   dplyr::group_by(region) %>%
#   dplyr::summarise(lat = mean(lat), lng = mean(long), .groups = "drop") %>%
#   dplyr::rename(country = region)

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

# Generate rainbow palettes for all possible clades found in any clade column
all_possible_clades <- sort(unique(c(
  metadata_global$clade, 
  metadata_global$G_clade, 
  metadata_global$HA_subclade, 
  metadata_global$HA_proposedSubclade, 
  metadata_global$HA_short_clade, 
  metadata_global$HA_legacy_clade
)))

# Exclude 'Unknown' from the main rainbow generation to give it a neutral color
actual_clades <- setdiff(all_possible_clades, "Unknown")
master_clade_colors <- grDevices::rainbow(length(actual_clades))
names(master_clade_colors) <- actual_clades

# Add neutral color for Unknown
master_clade_colors["Unknown"] <- "#d3d3d3"

# For backward compatibility with server logic that expects specific names
clade_colors_vec   <- master_clade_colors
g_clade_colors_vec <- master_clade_colors
