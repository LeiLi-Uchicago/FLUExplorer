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
# Version 10: Optimize Metadata Size & Types
RDS_CACHE <- "data/app_cache_flu_v9.rds"

# Subtypes to load: Dynamically detect valid subtype folders in data/
possible_dirs <- list.dirs("data", full.names = FALSE, recursive = FALSE)
SUBTYPES <- possible_dirs[sapply(possible_dirs, function(d) file.exists(file.path("data", d, "metadata_merged_annotated.csv")))]
if (length(SUBTYPES) > 0) {
  SUBTYPES <- c(sort(SUBTYPES[grepl("^H", SUBTYPES)]), sort(SUBTYPES[!grepl("^H", SUBTYPES)]))
} else {
  # Fallback just in case
  SUBTYPES <- c("H1N1", "H3N2", "B_VIC", "B_YAM")
}

if (file.exists(RDS_CACHE)) {
  # ---- FAST PATH: load everything from the pre-built cache ----
  message("Loading data from RDS cache: ", RDS_CACHE)
  cache <- readRDS(RDS_CACHE)
  
  # Check if the lightweight cache is present
  required_objects <- c("metadata_summary_stats", "important_pos_df", "metadata_grouping_cols")
  if (all(required_objects %in% names(cache))) {
    total_raw          <- cache$total_raw
    total_parsed       <- cache$total_parsed
    important_pos_df       <- cache$important_pos_df
    metadata_summary_stats <- cache$metadata_summary_stats
    total_countries_val    <- cache$total_countries_val
    time_range_val         <- cache$time_range_val
    metadata_groups        <- cache$metadata_groups
    metadata_years         <- cache$metadata_years
    metadata_grouping_cols <- cache$metadata_grouping_cols
    
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
  # ---- SLOW PATH: process raw CSVs into lazy-load RDS format and build cache ----
  message("RDS cache not found. Processing raw CSV files into lazy-load RDS format...")
  
  # --- Metadata ---
  all_metadata <- list()
  for (subtype in SUBTYPES) {
    meta_path <- paste0("data/", subtype, "/metadata_merged_annotated.csv")
    if (file.exists(meta_path)) {
      message("Loading metadata for ", subtype)
      # CRITICAL: na = character() ensures 'NA' (Neuraminidase) is NOT treated as a missing value
      meta <- read_csv(meta_path, col_types = cols(.default = "c"), show_col_types = FALSE, na = character())
      
      # Rename columns to match app expectations
      # Group = Subtype (A / H1N1)
      meta <- meta %>% 
        dplyr::rename(any_of(c(Group = "Subtype", 
                                 date = "Collection_Date",
                                 clade = "HA_clade", 
                                 G_clade = "NA_clade")))
      
      # Force Group name to match the folder exactly (bypassing messy fasta labels like "B / Victoria")
      meta$Group <- subtype
      
      # Parse Location into region and country
      loc_split <- stringr::str_split(meta$Location, " / ", simplify = TRUE)
      meta$region <- loc_split[, 1]
      meta$country <- loc_split[, 2]
      
      all_metadata[[subtype]] <- meta
    }
  }
  metadata_global <- bind_rows(all_metadata)
  
  # Identify dynamic grouping columns by matching with actual usage file groupings
  usage_groups <- c()
  for (subtype in SUBTYPES) {
    aa_dir <- paste0("data/", subtype, "/AA/")
    if (dir.exists(aa_dir)) {
      files <- list.files(aa_dir, pattern = "^aa_usage_by_.*\\.csv")
      usage_groups <- unique(c(usage_groups, sub("^aa_usage_by_(.*)\\.csv$", "\\1", files)))
    }
  }
  mapped_usage_groups <- usage_groups
  mapped_usage_groups[mapped_usage_groups == "HA_clade"] <- "clade"
  mapped_usage_groups[mapped_usage_groups == "NA_clade"] <- "G_clade"
  
  metadata_grouping_cols <- intersect(colnames(metadata_global), setdiff(mapped_usage_groups, c("Year", "Year_Month")))
  
  # Clean up empty clade names
  for(col in metadata_grouping_cols) {
    if(col %in% colnames(metadata_global)) {
      metadata_global[[col]] <- ifelse(is.na(metadata_global[[col]]) | metadata_global[[col]] == "" | metadata_global[[col]] == "trace 0", "Unknown", metadata_global[[col]])
    }
  }

  # Robust Date Handling
  metadata_global$Year <- as.numeric(stringr::str_extract(metadata_global$date, "^\\d{4}"))
  metadata_global$YM <- stringr::str_extract(metadata_global$date, "^\\d{4}-\\d{2}")
  
  total_raw    <- scales::comma(nrow(metadata_global))
  metadata_global <- metadata_global %>% filter(!is.na(Year))
  total_parsed <- scales::comma(nrow(metadata_global))
  
  # --- Process Usage Tables into individual RDS files ---
  for (subtype in SUBTYPES) {
    message("Processing usage tables into RDS for ", subtype)
    
    # AA tables
    aa_dir <- paste0("data/", subtype, "/AA/")
    if (dir.exists(aa_dir)) {
      aa_files <- list.files(aa_dir, pattern = "aa_usage_by_.*\\.csv", full.names = TRUE)
      for (f in aa_files) {
        df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
          dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
          mutate(Group = subtype)
          
        rds_path <- sub("\\.csv$", ".rds", f)
        saveRDS(df, rds_path)
      }
    }
    
    # NT tables
    nt_dir <- paste0("data/", subtype, "/NT/")
    if (dir.exists(nt_dir)) {
      nt_files <- list.files(nt_dir, pattern = "nt_usage_by_.*\\.csv", full.names = TRUE)
      for (f in nt_files) {
        df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
          dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
          dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide")) %>%
          mutate(Group = subtype)
          
        rds_path <- sub("\\.csv$", ".rds", f)
        saveRDS(df, rds_path)
      }
    }
  }

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
    group_by(across(all_of(c("Year", "Group", "region", "country", metadata_grouping_cols)))) %>%
    summarise(n = n(), .groups = "drop")
  
  total_countries_val <- length(unique(metadata_global$country))
  time_range_val <- paste(min(metadata_global$Year, na.rm=T), "-", max(metadata_global$Year, na.rm=T))
  
  raw_groups <- na.omit(unique(metadata_global$Group))
  metadata_groups <- c(sort(raw_groups[grepl("^H", raw_groups)]), sort(raw_groups[!grepl("^H", raw_groups)]))
  
  metadata_years <- sort(na.omit(unique(metadata_global$Year)), decreasing = TRUE)

  saveRDS(
    list(
      total_raw              = total_raw,
      total_parsed           = total_parsed,
      important_pos_df       = important_pos_df,
      metadata_summary_stats = metadata_summary_stats,
      total_countries_val    = total_countries_val,
      time_range_val         = time_range_val,
      metadata_groups        = metadata_groups,
      metadata_years         = metadata_years,
      metadata_grouping_cols = metadata_grouping_cols
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
if (!exists("metadata_grouping_cols")) {
  metadata_grouping_cols <- c("clade", "G_clade") # fallback
}

all_possible_clades <- unique(unlist(lapply(metadata_grouping_cols, function(col) metadata_summary_stats[[col]])))
all_possible_clades <- sort(na.omit(all_possible_clades))

# Exclude 'Unknown' from the main rainbow generation to give it a neutral color
actual_clades <- setdiff(all_possible_clades, "Unknown")
master_clade_colors <- viridis::viridis(length(actual_clades))
names(master_clade_colors) <- actual_clades

# Add neutral color for Unknown
master_clade_colors["Unknown"] <- "#d3d3d3"

# For backward compatibility with server logic that expects specific names
clade_colors_vec   <- master_clade_colors
g_clade_colors_vec <- master_clade_colors

# ==========================================
# 3. LAZY-LOADING CACHE MECHANISM (LRU)
# ==========================================
lazy_cache <- new.env(parent = emptyenv())
lazy_cache$keys <- character(0)
lazy_cache$data <- list()

get_lazy_table <- function(rds_path, max_tables = 5) {
  if (!file.exists(rds_path)) return(NULL)
  
  if (rds_path %in% lazy_cache$keys) {
    # Move to the end of the line (most recently used)
    lazy_cache$keys <- c(setdiff(lazy_cache$keys, rds_path), rds_path)
    return(lazy_cache$data[[rds_path]])
  }
  
  # Read from disk and store in cache
  df <- readRDS(rds_path)
  lazy_cache$keys <- c(lazy_cache$keys, rds_path)
  lazy_cache$data[[rds_path]] <- df
  
  # Evict oldest if limit exceeded
  if (length(lazy_cache$keys) > max_tables) {
    evict <- lazy_cache$keys[1]
    lazy_cache$keys <- lazy_cache$keys[-1]
    lazy_cache$data[[evict]] <- NULL
  }
  
  return(df)
}
