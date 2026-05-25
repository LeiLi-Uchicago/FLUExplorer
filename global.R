library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(readr)
library(tidyr)
library(openxlsx)
library(plotly)
# library(msaR) 
# library(Biostrings)
# library(msa)
library(waiter)
library(lubridate)
library(tidyverse)
# library(leaflet)             # Fixes: could not find function "leafletOutput"
# library(leaflet.minicharts)  # Required for the pie charts on the map
library(shinyWidgets)
library(shinyjs)

USE_DUCKDB <- requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)
if (!USE_DUCKDB) {
  message("Package 'duckdb' is not installed. Falling back to legacy RDS lazy loading.")
}

# Disable scientific notation for the session
options(scipen = 999)

# ==========================================
# 1. GLOBAL DATA LOADING & SETUP
# ==========================================
# Version 11: Read important_positions.csv
RDS_CACHE <- "data/app_cache_flu.rds"
DUCKDB_CACHE <- "data/flu_explorer.duckdb"
DUCKDB_META_CACHE <- "data/flu_explorer_duckdb_meta.rds"
CACHE_SCHEMA_VERSION <- 4L
RAW_DATA_DIR <- "data/raw"
COUNT_RDS_CACHE_DIR <- "data/count_cache"
VALIDATION_ONLY_COUNT_COLS <- c("CodonStatus", "CodonSource")

metadata_file_path <- function(subtype) {
  file.path(RAW_DATA_DIR, subtype, "metadata_merged_annotated.csv")
}

count_root_path <- function(subtype, var_type = "AA") {
  raw_count_root <- file.path(RAW_DATA_DIR, subtype, "count")
  legacy_root <- file.path("data", subtype, var_type)
  if (dir.exists(raw_count_root)) raw_count_root else legacy_root
}

count_cache_root_path <- function(subtype, var_type = "AA") {
  file.path(COUNT_RDS_CACHE_DIR, subtype, var_type)
}

count_cache_gene_path <- function(subtype, var_type, gene) {
  file.path(count_cache_root_path(subtype, var_type), gene)
}

count_cache_file_path <- function(subtype, var_type, gene, group_by) {
  file.path(count_cache_gene_path(subtype, var_type, gene), paste0(tolower(var_type), "_usage_by_", group_by, ".rds"))
}

raw_count_file_path <- function(subtype, var_type, gene, group_by) {
  file.path(count_root_path(subtype, var_type), gene, paste0(tolower(var_type), "_usage_by_", group_by, ".csv"))
}

count_gene_dirs <- function(subtype, var_type = "AA", prefer_cache = FALSE) {
  prefix <- paste0(tolower(var_type), "_usage_by_")
  roots <- if (isTRUE(prefer_cache)) {
    c(count_cache_root_path(subtype, var_type), count_root_path(subtype, var_type))
  } else {
    c(count_root_path(subtype, var_type), count_cache_root_path(subtype, var_type))
  }
  roots <- roots[dir.exists(roots)]
  if (length(roots) == 0) return(character(0))
  dirs <- unlist(lapply(roots, list.dirs, full.names = TRUE, recursive = FALSE), use.names = FALSE)
  dirs <- dirs[basename(dirs) != ".DS_Store"]
  dirs <- dirs[vapply(dirs, function(path) {
    length(list.files(path, pattern = paste0("^", prefix, ".*\\.(csv|rds)$"))) > 0
  }, logical(1))]
  dirs[!duplicated(basename(dirs))]
}

available_count_genes <- function(subtype, var_type = "AA", prefer_cache = FALSE) {
  sort(unique(basename(count_gene_dirs(subtype, var_type, prefer_cache = prefer_cache))))
}

latest_raw_data_mtime <- function() {
  files <- list.files(RAW_DATA_DIR, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(as.POSIXct(NA))
  max(file.info(files)$mtime, na.rm = TRUE)
}

cache_older_than_raw_data <- function(cache_path) {
  if (!file.exists(cache_path)) return(TRUE)
  raw_mtime <- latest_raw_data_mtime()
  if (is.na(raw_mtime)) return(FALSE)
  cache_mtime <- file.info(cache_path)$mtime
  is.na(cache_mtime) || cache_mtime < raw_mtime
}

strip_validation_count_cols <- function(df) {
  df[, setdiff(names(df), VALIDATION_ONLY_COUNT_COLS), drop = FALSE]
}

parse_metadata_year_month <- function(metadata) {
  year_values <- if ("Year" %in% names(metadata)) as.character(metadata$Year) else rep(NA_character_, nrow(metadata))
  month_values <- if ("Month" %in% names(metadata)) as.character(metadata$Month) else rep(NA_character_, nrow(metadata))

  parsed_year <- suppressWarnings(as.numeric(trimws(year_values)))
  month_clean <- trimws(month_values)
  month_clean[is.na(month_clean) | month_clean == ""] <- NA_character_
  month_clean <- ifelse(
    grepl("^[0-9]+$", month_clean),
    stringr::str_pad(month_clean, width = 2, side = "left", pad = "0"),
    month_clean
  )

  metadata$Year <- parsed_year
  metadata$Month <- month_clean
  metadata$YM <- ifelse(!is.na(parsed_year) & !is.na(month_clean), paste0(parsed_year, "-", month_clean), NA_character_)
  metadata
}

empty_metadata_clade_explorer <- function() {
  list(
    month_totals = tibble::tibble(),
    summaries = tibble::tibble(),
    monthly = tibble::tibble(),
    breakdowns = tibble::tibble()
  )
}

is_unknown_metadata_value <- function(value) {
  text <- stringr::str_squish(as.character(value))
  is.na(text) | text == "" | stringr::str_to_lower(text) %in% c(
    "unknown", "na", "n/a", "none", "null", "?", "unassigned", "not assigned",
    "trace 0", "undetermined", "not determined"
  )
}

build_metadata_clade_explorer_summary <- function(metadata, annotation_cols) {
  annotation_cols <- annotation_cols[annotation_cols %in% names(metadata)]

  for (column in c(annotation_cols, "Group", "YM", "region", "country", "Host")) {
    if (!column %in% names(metadata)) metadata[[column]] <- NA_character_
  }

  month_totals <- metadata %>%
    filter(!is.na(.data$YM), .data$YM != "") %>%
    count(Group, YearMonth = .data$YM, name = "Total") %>%
    arrange(Group, YearMonth)

  if (length(annotation_cols) == 0 || nrow(metadata) == 0) {
    out <- empty_metadata_clade_explorer()
    out$month_totals <- month_totals
    return(out)
  }

  subtype_totals <- metadata %>%
    count(Group, name = "TotalSequences")

  base <- purrr::map_dfr(annotation_cols, function(annotation) {
    metadata %>%
      transmute(
        Group = as.character(.data$Group),
        Annotation = annotation,
        Clade = as.character(.data[[annotation]]),
        YearMonth = as.character(.data$YM),
        region = as.character(.data$region),
        country = as.character(.data$country),
        host = as.character(.data$Host)
      ) %>%
      filter(!is_unknown_metadata_value(.data$Clade))
  })

  if (nrow(base) == 0) {
    out <- empty_metadata_clade_explorer()
    out$month_totals <- month_totals
    return(out)
  }

  monthly <- base %>%
    filter(!is.na(.data$YearMonth), .data$YearMonth != "") %>%
    count(Group, Annotation, Clade, YearMonth, name = "Count") %>%
    left_join(month_totals, by = c("Group", "YearMonth")) %>%
    mutate(
      Total = dplyr::coalesce(.data$Total, 0L),
      Percent = dplyr::if_else(.data$Total > 0, (.data$Count / .data$Total) * 100, 0)
    ) %>%
    arrange(Group, Annotation, Clade, YearMonth)

  annotation_totals <- base %>%
    count(Group, Annotation, name = "AnnotatedTotal") %>%
    left_join(subtype_totals, by = "Group") %>%
    mutate(
      MissingAnnotationCount = pmax(.data$TotalSequences - .data$AnnotatedTotal, 0),
      AnnotationCoverage = dplyr::if_else(.data$TotalSequences > 0, (.data$AnnotatedTotal / .data$TotalSequences) * 100, 0)
    )

  totals <- base %>%
    count(Group, Annotation, Clade, name = "StrainCount") %>%
    left_join(annotation_totals, by = c("Group", "Annotation")) %>%
    arrange(Group, Annotation, desc(StrainCount), Clade) %>%
    group_by(Group, Annotation) %>%
    mutate(
      Rank = row_number(),
      DatasetShare = dplyr::if_else(.data$AnnotatedTotal > 0, (.data$StrainCount / .data$AnnotatedTotal) * 100, 0),
      TotalDatasetShare = dplyr::if_else(.data$TotalSequences > 0, (.data$StrainCount / .data$TotalSequences) * 100, 0)
    ) %>%
    ungroup()

  periods <- base %>%
    filter(!is.na(.data$YearMonth), .data$YearMonth != "") %>%
    group_by(Group, Annotation, Clade) %>%
    summarise(
      FirstMonth = min(.data$YearMonth),
      LastMonth = max(.data$YearMonth),
      .groups = "drop"
    )

  peaks <- monthly %>%
    group_by(Group, Annotation, Clade) %>%
    arrange(desc(.data$Percent), desc(.data$Count), .data$YearMonth, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      Group,
      Annotation,
      Clade,
      PeakMonth = YearMonth,
      PeakCount = Count,
      PeakPercent = Percent
    )

  summaries <- totals %>%
    left_join(periods, by = c("Group", "Annotation", "Clade")) %>%
    left_join(peaks, by = c("Group", "Annotation", "Clade")) %>%
    arrange(Group, Annotation, Rank)

  breakdowns <- purrr::map_dfr(c("country", "region", "host"), function(category) {
    base %>%
      filter(!is_unknown_metadata_value(.data[[category]])) %>%
      count(Group, Annotation, Clade, Category = category, Value = .data[[category]], name = "Count") %>%
      group_by(Group, Annotation, Clade, Category) %>%
      arrange(desc(.data$Count), .data$Value, .by_group = TRUE) %>%
      mutate(Rank = row_number()) %>%
      filter(.data$Rank <= 10) %>%
      ungroup()
  })

  list(
    month_totals = month_totals,
    summaries = summaries,
    monthly = monthly,
    breakdowns = breakdowns
  )
}

raw_metadata_available <- function(subtypes = SUBTYPES) {
  length(subtypes) > 0 && any(file.exists(vapply(subtypes, metadata_file_path, character(1))))
}

normalize_year_month_filter <- function(year, month) {
  year_chr <- trimws(as.character(year))
  month_chr <- trimws(as.character(month))
  month_chr <- ifelse(grepl("^[0-9]+$", month_chr), stringr::str_pad(month_chr, width = 2, side = "left", pad = "0"), month_chr)

  dplyr::case_when(
    year_chr %in% c("", "NA") | month_chr %in% c("", "NA") ~ NA_character_,
    year_chr %in% c("Unknown", "unassigned", "Unassigned") ~ year_chr,
    month_chr %in% c("Unknown", "unassigned", "Unassigned") ~ month_chr,
    TRUE ~ paste0(year_chr, "-", month_chr)
  )
}

empty_usage_duckdb_df <- function() {
  data.frame(
    Group = character(),
    Variation_Type = character(),
    Gene = character(),
    Grouping_Type = character(),
    Clade = character(),
    Position = numeric(),
    AminoAcid = character(),
    Count = numeric(),
    Year = character(),
    Month = character(),
    Year_Month = character(),
    Year_Month_Filter = character(),
    Codon_Usage = character(),
    stringsAsFactors = FALSE
  )
}

normalize_usage_table <- function(df, subtype, var_type, gene_name, group_name) {
  df <- strip_validation_count_cols(df)
  df <- df %>%
    dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
    mutate(Group = subtype)

  if (var_type == "NT") {
    df <- df %>%
      dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide"))
  }

  if (!"Gene" %in% names(df)) df$Gene <- gene_name
  if (!"AminoAcid" %in% names(df)) df$AminoAcid <- NA_character_
  if (!"Count" %in% names(df)) df$Count <- NA_real_
  if (!"Position" %in% names(df)) df$Position <- NA_real_
  if (!"Year" %in% names(df)) df$Year <- NA_character_
  if (!"Month" %in% names(df)) df$Month <- NA_character_
  if (!"Year_Month" %in% names(df)) {
    if (group_name == "Year_Month") {
      df$Year_Month <- normalize_year_month_filter(df$Year, df$Month)
    } else {
      df$Year_Month <- NA_character_
    }
  }
  if (!group_name %in% names(df)) df[[group_name]] <- "Unknown"
  if (group_name == "Year_Month") df[[group_name]] <- normalize_year_month_filter(df$Year, df$Month)
  if (!"Codon_Usage" %in% names(df) && "Codon" %in% names(df)) df$Codon_Usage <- df$Codon
  if (!"Codon_Usage" %in% names(df)) df$Codon_Usage <- NA_character_

  data.frame(
    Group = as.character(df$Group),
    Variation_Type = as.character(var_type),
    Gene = as.character(df$Gene),
    Grouping_Type = as.character(group_name),
    Clade = as.character(df[[group_name]]),
    Position = suppressWarnings(as.numeric(df$Position)),
    AminoAcid = as.character(df$AminoAcid),
    Count = suppressWarnings(as.numeric(df$Count)),
    Year = as.character(df$Year),
    Month = as.character(df$Month),
    Year_Month = as.character(df$Year_Month),
    Year_Month_Filter = normalize_year_month_filter(df$Year, df$Month),
    Codon_Usage = as.character(df$Codon_Usage),
    stringsAsFactors = FALSE
  )
}

build_usage_duckdb_cache <- function(subtypes = SUBTYPES, db_path = DUCKDB_CACHE) {
  if (!USE_DUCKDB) return(FALSE)

  tmp_path <- paste0(db_path, ".tmp")
  if (file.exists(tmp_path)) unlink(tmp_path)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp_path, read_only = FALSE)
  on.exit({
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  DBI::dbExecute(con, "PRAGMA memory_limit='700MB'")
  DBI::dbWriteTable(con, "usage", empty_usage_duckdb_df(), overwrite = TRUE)

  usage_groups <- character()

  for (subtype in subtypes) {
    message("Processing usage tables into DuckDB for ", subtype)

    for (var_type in c("AA", "NT")) {
      var_root <- count_root_path(subtype, var_type)
      if (!dir.exists(var_root)) next

      gene_dirs <- count_gene_dirs(subtype, var_type)
      for (g_dir in gene_dirs) {
        gene_name <- basename(g_dir)
        message("  Processing ", var_type, " / ", gene_name)

        prefix <- if (var_type == "AA") "aa" else "nt"
        pattern <- paste0("^", prefix, "_usage_by_.*\\.csv$")
        files <- list.files(g_dir, pattern = pattern, full.names = TRUE)

        for (f in files) {
          group_name <- sub(paste0("^", prefix, "_usage_by_(.*)\\.csv$"), "\\1", basename(f))
          usage_groups <- unique(c(usage_groups, group_name))

          df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
            normalize_usage_table(subtype, var_type, gene_name, group_name)

          DBI::dbWriteTable(con, "usage", df, append = TRUE)

          if (group_name == "Year_Month") {
            year_df <- df %>%
              group_by(Group, Variation_Type, Gene, Position, AminoAcid, Year) %>%
              summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
              mutate(
                Grouping_Type = "Year",
                Clade = as.character(Year),
                Month = NA_character_,
                Year_Month = NA_character_,
                Year_Month_Filter = NA_character_,
                Codon_Usage = NA_character_
              ) %>%
              dplyr::select(names(empty_usage_duckdb_df()))

            DBI::dbWriteTable(con, "usage", year_df, append = TRUE)
            usage_groups <- unique(c(usage_groups, "Year"))
          }

          rm(df)
          gc(FALSE)
        }
      }
    }
  }

  if (identical(tolower(Sys.getenv("FLUEXPLORER_DUCKDB_CREATE_INDEXES", "false")), "true")) {
    DBI::dbExecute(con, "CREATE INDEX idx_usage_main ON usage (\"Group\", Variation_Type, Gene, Grouping_Type, Position)")
    DBI::dbExecute(con, "CREATE INDEX idx_usage_clade ON usage (\"Group\", Variation_Type, Gene, Grouping_Type, Clade)")
    DBI::dbExecute(con, "CREATE INDEX idx_usage_time ON usage (\"Group\", Variation_Type, Gene, Grouping_Type, Position, Year_Month_Filter)")
  } else {
    message("Skipping DuckDB index creation. Set FLUEXPLORER_DUCKDB_CREATE_INDEXES=true to enable it during cache builds.")
  }

  saveRDS(list(usage_groups = usage_groups, built_at = Sys.time()), DUCKDB_META_CACHE)
  DBI::dbDisconnect(con, shutdown = TRUE)

  if (file.exists(db_path)) unlink(db_path)
  file.rename(tmp_path, db_path)
}

ensure_usage_duckdb_cache <- function() {
  if (!USE_DUCKDB) return(FALSE)
  if (file.exists(DUCKDB_CACHE)) return(TRUE)

  message("DuckDB cache missing or stale. Building: ", DUCKDB_CACHE)
  ok <- tryCatch(build_usage_duckdb_cache(), error = function(e) {
    message("DuckDB cache build failed: ", conditionMessage(e))
    FALSE
  })
  isTRUE(ok) && file.exists(DUCKDB_CACHE)
}

# Subtypes to load: Dynamically detect valid subtype folders in data/raw/
possible_dirs <- if (dir.exists(RAW_DATA_DIR)) {
  list.dirs(RAW_DATA_DIR, full.names = FALSE, recursive = FALSE)
} else {
  character(0)
}
has_metadata_file <- vapply(possible_dirs, function(d) file.exists(metadata_file_path(d)), logical(1))
SUBTYPES <- possible_dirs[has_metadata_file]
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
    if (
      (is.null(cache$cache_schema_version) || !identical(as.integer(cache$cache_schema_version), CACHE_SCHEMA_VERSION)) &&
        raw_metadata_available()
    ) {
      message("RDS cache schema is outdated. Rebuilding from raw metadata...")
      cache_loaded <- FALSE
    } else if (!"metadata_clade_explorer" %in% names(cache) && raw_metadata_available()) {
      message("RDS cache is missing Genetic Clade summaries. Rebuilding from raw metadata...")
      cache_loaded <- FALSE
    } else {
      total_raw          <- cache$total_raw
      total_parsed       <- cache$total_parsed
      important_pos_df       <- cache$important_pos_df
      metadata_summary_stats <- cache$metadata_summary_stats
      total_countries_val    <- cache$total_countries_val
      time_range_val         <- cache$time_range_val
      metadata_groups        <- cache$metadata_groups
      metadata_years         <- cache$metadata_years
      metadata_grouping_cols <- cache$metadata_grouping_cols
      metadata_clade_explorer <- cache$metadata_clade_explorer
      if (is.null(metadata_clade_explorer)) {
        metadata_clade_explorer <- empty_metadata_clade_explorer()
      }

      rm(cache)   # free the wrapper list from memory
      message("RDS cache loaded successfully.")
      cache_loaded <- TRUE
    }
  } else {
    message("RDS cache is outdated. Rebuilding...")
    cache_loaded <- FALSE
  }
} else {
  cache_loaded <- FALSE
}

if (!cache_loaded) {
  # ---- SLOW PATH: process raw CSVs and build cache ----
  message("RDS cache not found. Processing raw CSV files...")
  
  # --- Metadata ---
  all_metadata <- list()
  for (subtype in SUBTYPES) {
    meta_path <- metadata_file_path(subtype)
    if (file.exists(meta_path)) {
      message("Loading metadata for ", subtype)
      # CRITICAL: na = character() ensures 'NA' (Neuraminidase) is NOT treated as a missing value
      meta <- read_csv(meta_path, col_types = cols(.default = "c"), show_col_types = FALSE, na = character())
      
      # Rename columns to match app expectations
      # Group = Subtype (A / H1N1)
      meta <- meta %>% 
        dplyr::rename(any_of(c(Group = "Subtype",
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
    aa_root <- count_root_path(subtype, "AA")
    if (dir.exists(aa_root)) {
      # Look into the first gene subdirectory found
      gene_dirs <- count_gene_dirs(subtype, "AA")
      if (length(gene_dirs) > 0) {
        files <- list.files(gene_dirs[1], pattern = "^aa_usage_by_.*\\.csv")
        usage_groups <- unique(c(usage_groups, sub("^aa_usage_by_(.*)\\.csv$", "\\1", files)))
      }
    }
  }
  # Add Year explicitly since we will generate it from Year_Month
  if ("Year_Month" %in% usage_groups) {
    usage_groups <- unique(c(usage_groups, "Year"))
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

  # Metadata time handling: use the explicit parsed Year/Month columns from
  # metadata_merged_annotated.csv. Do not parse Collection_Date here.
  metadata_global <- parse_metadata_year_month(metadata_global)
  
  total_raw    <- scales::comma(nrow(metadata_global))
  metadata_global <- metadata_global %>% filter(!is.na(Year))
  total_parsed <- scales::comma(nrow(metadata_global))
  
  # --- Process Usage Tables ---
  duckdb_ready <- ensure_usage_duckdb_cache()
  if (!duckdb_ready) {
    message("DuckDB is unavailable. Processing usage tables into legacy RDS files.")

    for (subtype in SUBTYPES) {
      message("Processing usage tables into RDS for ", subtype)

      for (var_type in c("AA", "NT")) {
        var_root <- count_root_path(subtype, var_type)
        if (!dir.exists(var_root)) next

        gene_dirs <- count_gene_dirs(subtype, var_type)
        for (g_dir in gene_dirs) {
          gene_name <- basename(g_dir)
          message("  Processing ", var_type, " / ", gene_name)

          prefix <- if(var_type == "AA") "aa" else "nt"
          pattern <- paste0("^", prefix, "_usage_by_.*\\.csv$")
          files <- list.files(g_dir, pattern = pattern, full.names = TRUE)

          for (f in files) {
            is_ym <- grepl(paste0(prefix, "_usage_by_Year_Month\\.csv$"), f)
            group_name <- sub(paste0("^", prefix, "_usage_by_(.*)\\.csv$"), "\\1", basename(f))

            df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
              strip_validation_count_cols() %>%
              dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
              mutate(Group = subtype)

            if (var_type == "NT") {
              df <- df %>%
                dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide"))
            }

            if (is_ym) {
              df <- df %>%
                mutate(Year_Month = normalize_year_month_filter(Year, Month))
            }

            if (!"Codon_Usage" %in% names(df) && "Codon" %in% names(df)) df$Codon_Usage <- df$Codon

            out_dir <- count_cache_gene_path(subtype, var_type, gene_name)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            saveRDS(df, file.path(out_dir, paste0(prefix, "_usage_by_", group_name, ".rds")))

            if (is_ym) {
              year_df <- df %>%
                group_by(Group, Gene, Position, Year, AminoAcid) %>%
                summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

              saveRDS(year_df, file.path(out_dir, paste0(prefix, "_usage_by_Year.rds")))
            }
          }
        }
      }
    }
  }

  # --- Important positions ---
  important_pos_file <- "data/important_positions.csv"
  if (file.exists(important_pos_file)) {
    message("Loading important positions from CSV...")
    important_pos_df <- read_csv(important_pos_file, show_col_types = FALSE)
    if (!"label" %in% names(important_pos_df)) {
      important_pos_df$label <- paste(important_pos_df$Subtype, important_pos_df$Gene, important_pos_df$Position, sep = " - ")
    }
  } else {
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
  }
  
  # --- Save RDS cache ---
  message("Writing RDS cache to: ", RDS_CACHE)
  
  # Pre-calculate these before saving if they weren't loaded from cache
  message("Pre-aggregating metadata statistics...")
  metadata_summary_stats <- metadata_global %>%
    mutate(YearMonth = .data$YM) %>%
    group_by(across(all_of(c("Year", "YearMonth", "Group", "region", "country", metadata_grouping_cols)))) %>%
    summarise(n = n(), .groups = "drop")

  metadata_clade_explorer <- build_metadata_clade_explorer_summary(metadata_global, metadata_grouping_cols)
  
  total_countries_val <- length(unique(metadata_global$country))
  time_range_val <- paste(min(metadata_global$Year, na.rm=T), "-", max(metadata_global$Year, na.rm=T))
  
  raw_groups <- na.omit(unique(metadata_global$Group))
  metadata_groups <- c(sort(raw_groups[grepl("^H", raw_groups)]), sort(raw_groups[!grepl("^H", raw_groups)]))
  
  metadata_years <- sort(na.omit(unique(metadata_global$Year)), decreasing = TRUE)

  saveRDS(
    list(
      cache_schema_version   = CACHE_SCHEMA_VERSION,
      total_raw              = total_raw,
      total_parsed           = total_parsed,
      important_pos_df       = important_pos_df,
      metadata_summary_stats = metadata_summary_stats,
      total_countries_val    = total_countries_val,
      time_range_val         = time_range_val,
      metadata_groups        = metadata_groups,
      metadata_years         = metadata_years,
      metadata_grouping_cols = metadata_grouping_cols,
      metadata_clade_explorer = metadata_clade_explorer
    ),
    file = RDS_CACHE
  )
  
  # --- STARTUP MEMORY FLUSH ---
  suppressWarnings(rm(all_metadata, metadata_global, meta))
  gc(verbose = FALSE)
}

if (!exists("metadata_clade_explorer") || is.null(metadata_clade_explorer)) {
  metadata_clade_explorer <- empty_metadata_clade_explorer()
}

duckdb_cache_ready <- ensure_usage_duckdb_cache()

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

LAZY_CACHE_MAX_TABLES <- 2
LAZY_CACHE_MAX_MEM_MB <- 450

get_lazy_table <- function(rds_path, max_tables = LAZY_CACHE_MAX_TABLES, max_mem_mb = LAZY_CACHE_MAX_MEM_MB) {
  if (!file.exists(rds_path)) return(NULL)
  
  if (rds_path %in% lazy_cache$keys) {
    # Move to the end of the line (most recently used)
    lazy_cache$keys <- c(setdiff(lazy_cache$keys, rds_path), rds_path)
    return(lazy_cache$data[[rds_path]])
  }
  
  # Check memory usage against threshold before attempting to load new data
  if (sum(gc(verbose = FALSE)[, 2]) > max_mem_mb) {
    message("Memory usage exceeds ", max_mem_mb, " MB threshold. Clearing cache...")
    lazy_cache$keys <- character(0)
    lazy_cache$data <- list()
    gc(verbose = FALSE)
  }

  # Read from disk and store in cache
  df <- readRDS(rds_path)
  lazy_cache$keys <- c(lazy_cache$keys, rds_path)
  lazy_cache$data[[rds_path]] <- df
  
  # Evict oldest if limit exceeded
  while (length(lazy_cache$keys) > max_tables) {
    evict <- lazy_cache$keys[1]
    lazy_cache$keys <- lazy_cache$keys[-1]
    lazy_cache$data[[evict]] <- NULL
    
    # Force memory release back to the OS
    gc(verbose = FALSE)
  }
  
  return(df)
}

# ==========================================
# 4. DUCKDB QUERY HELPERS
# ==========================================
usage_db_env <- new.env(parent = emptyenv())
usage_db_env$con <- NULL

usage_duckdb_available <- function() {
  isTRUE(USE_DUCKDB) && isTRUE(duckdb_cache_ready) && file.exists(DUCKDB_CACHE)
}

usage_db_conn <- function() {
  if (!usage_duckdb_available()) return(NULL)

  if (is.null(usage_db_env$con) || !DBI::dbIsValid(usage_db_env$con)) {
    usage_db_env$con <- DBI::dbConnect(duckdb::duckdb(), dbdir = DUCKDB_CACHE, read_only = TRUE)
    DBI::dbExecute(usage_db_env$con, "PRAGMA memory_limit='700MB'")
  }

  usage_db_env$con
}

usage_query <- function(sql, params = NULL) {
  con <- usage_db_conn()
  if (is.null(con)) return(NULL)
  DBI::dbGetQuery(con, sql, params = params)
}

usage_sql_in_values <- function(values) {
  con <- usage_db_conn()
  if (is.null(con) || length(values) == 0) return(NULL)
  paste(DBI::dbQuoteString(con, values), collapse = ", ")
}

usage_file_groups <- function(subtype, var_type, gene) {
  dirs <- c(
    count_cache_gene_path(subtype, var_type, gene),
    file.path(count_root_path(subtype, var_type), gene)
  )
  files <- unlist(lapply(dirs[dir.exists(dirs)], function(dir_path) {
    list.files(dir_path, pattern = paste0("^", tolower(var_type), "_usage_by_.*\\.(rds|csv)$"))
  }), use.names = FALSE)
  sort(unique(sub(paste0("^", tolower(var_type), "_usage_by_(.*)\\.(rds|csv)$"), "\\1", files)))
}

usage_available_genes <- function(subtype, var_type) {
  if (!usage_duckdb_available()) return(available_count_genes(subtype, var_type, prefer_cache = TRUE))

  res <- usage_query(
    "SELECT DISTINCT Gene FROM usage WHERE \"Group\" = ? AND Variation_Type = ?",
    list(subtype, var_type)
  )
  if (is.null(res) || nrow(res) == 0) return(available_count_genes(subtype, var_type, prefer_cache = TRUE))
  sort(stats::na.omit(as.character(res$Gene)))
}

usage_available_groups <- function(subtype, var_type, gene) {
  if (!usage_duckdb_available()) return(usage_file_groups(subtype, var_type, gene))

  res <- usage_query(
    "SELECT DISTINCT Grouping_Type FROM usage WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ?",
    list(subtype, var_type, gene)
  )
  groups <- sort(as.character(res$Grouping_Type))

  if (length(groups) == 0) usage_file_groups(subtype, var_type, gene) else groups
}

usage_distinct_group_values <- function(subtype, var_type, gene, group_by) {
  res <- usage_query(
    "SELECT DISTINCT
       CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade
     FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?",
    list(subtype, var_type, gene, group_by)
  )
  if (is.null(res) || nrow(res) == 0) return(character(0))

  values <- sort(stats::na.omit(as.character(res$Clade)))
  special_values <- c("Unknown", "unassigned", "Unassigned")
  present_specials <- intersect(special_values, values)
  if (length(present_specials) > 0) values <- c(setdiff(values, present_specials), present_specials)
  values
}

usage_max_position <- function(subtype, var_type, gene) {
  res <- usage_query(
    "SELECT MAX(Position) AS max_position FROM usage WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ?",
    list(subtype, var_type, gene)
  )
  if (is.null(res) || nrow(res) == 0 || is.na(res$max_position[1])) return(NA_real_)
  as.numeric(res$max_position[1])
}

usage_year_month_choices <- function(subtype, var_type, gene, group_by, position) {
  res <- usage_query(
    "SELECT DISTINCT Year_Month_Filter FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
       AND Position = ? AND Year_Month_Filter IS NOT NULL",
    list(subtype, var_type, gene, group_by, as.numeric(position))
  )
  if (is.null(res) || nrow(res) == 0) return(character(0))

  ym_values <- stats::na.omit(as.character(res$Year_Month_Filter))
  special_values <- c("Unknown", "unassigned", "Unassigned")
  present_specials <- intersect(special_values, ym_values)
  chronological_yms <- sort(setdiff(ym_values, special_values))
  c(present_specials, chronological_yms)
}

usage_single_position <- function(subtype, var_type, gene, group_by, position, allowed_yms = NULL, min_seqs = 1, hide_empty_years = FALSE) {
  ym_filter <- ""
  if (!is.null(allowed_yms) && length(allowed_yms) > 0) {
    in_values <- usage_sql_in_values(allowed_yms)
    if (!is.null(in_values)) {
      ym_filter <- paste0(" AND Year_Month_Filter IN (", in_values, ")")
    }
  }

  sql <- paste0(
    "SELECT \"Group\", Gene, Position,
            CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
            AminoAcid, SUM(Count) AS Count,
            ANY_VALUE(Codon_Usage) AS Codon_Usage
     FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
       AND Position = ? AND AminoAcid NOT IN ('X', '-')",
    ym_filter,
    " GROUP BY \"Group\", Gene, Position,
       CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
       AminoAcid"
  )
  res <- usage_query(sql, list(subtype, var_type, gene, group_by, as.numeric(position)))
  if (is.null(res)) return(NULL)
  if (nrow(res) == 0) {
    out <- data.frame(Group=character(), Gene=character(), Position=numeric(), AminoAcid=character(), Count=numeric(), Valid_Total=numeric(), `Frequency(%)`=numeric(), check.names = FALSE)
    out[[group_by]] <- character()
    return(out)
  }

  res <- res %>%
    dplyr::rename(!!group_by := Clade) %>%
    group_by(.data[[group_by]]) %>%
    mutate(
      Valid_Total = sum(Count, na.rm = TRUE),
      `Frequency(%)` = (Count / Valid_Total) * 100
    ) %>%
    ungroup() %>%
    filter(Valid_Total >= min_seqs)

  if (group_by == "Year" && isTRUE(hide_empty_years)) {
    res <- res %>% filter(Valid_Total > 0)
  }

  if (all(is.na(res$Codon_Usage))) {
    res$Codon_Usage <- NULL
  }

  res
}

usage_pairwise_gene_data <- function(subtype, var_type, gene, group_by, clades = NULL) {
  clade_filter <- ""
  if (!is.null(clades) && length(clades) > 0) {
    in_values <- usage_sql_in_values(clades)
    if (!is.null(in_values)) {
      clade_filter <- paste0(
        " AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END IN (",
        in_values,
        ")"
      )
    }
  }

  sql <- paste0(
    "SELECT \"Group\", Gene,
            CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
            Position, AminoAcid, SUM(Count) AS Count,
            ANY_VALUE(Codon_Usage) AS Codon_Usage
     FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?",
    clade_filter,
    " GROUP BY \"Group\", Gene,
       CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
       Position, AminoAcid"
  )
  res <- usage_query(sql, list(subtype, var_type, gene, group_by))
  if (is.null(res)) return(NULL)
  if (all(is.na(res$Codon_Usage))) res$Codon_Usage <- NULL
  res
}

usage_pairwise_differences_for_gene <- function(subtype, var_type, gene, group_by, clade1, clade2, min_freq) {
  res <- usage_query(
    "WITH agg AS (
       SELECT Gene, Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
         AminoAcid, SUM(Count) AS Variant_Count
       FROM usage
       WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
         AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END IN (?, ?)
         AND AminoAcid NOT IN ('X', '-')
       GROUP BY Gene, Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
         AminoAcid
     ),
     freq AS (
       SELECT *,
         SUM(Variant_Count) OVER (PARTITION BY Gene, Position, Clade) AS Total_Seqs,
         100.0 * Variant_Count / SUM(Variant_Count) OVER (PARTITION BY Gene, Position, Clade) AS Freq
       FROM agg
     ),
     ranked AS (
       SELECT *,
         ROW_NUMBER() OVER (PARTITION BY Gene, Position, Clade ORDER BY Freq DESC, AminoAcid) AS rn
       FROM freq
     )
     SELECT Gene, Position, Clade, AminoAcid, Freq
     FROM ranked
     WHERE rn = 1 AND Freq >= ?",
    list(subtype, var_type, gene, group_by, clade1, clade2, as.numeric(min_freq))
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)

  c1_dom <- res %>%
    filter(Clade == clade1) %>%
    dplyr::select(Gene, Position, Clade1_AA = AminoAcid, Clade1_Freq = Freq)
  c2_dom <- res %>%
    filter(Clade == clade2) %>%
    dplyr::select(Gene, Position, Clade2_AA = AminoAcid, Clade2_Freq = Freq)

  inner_join(c1_dom, c2_dom, by = c("Gene", "Position")) %>%
    filter(Clade1_AA != Clade2_AA)
}

usage_position_distribution <- function(subtype, var_type, gene, group_by, position, hide_empty_years = FALSE) {
  res <- usage_pairwise_gene_data(subtype, var_type, gene, group_by)
  if (is.null(res)) return(NULL)

  res <- res %>%
    filter(Position == position, !(AminoAcid %in% c("X", "-")))

  if (is.null(res) || nrow(res) == 0) return(NULL)

  has_codon <- "Codon_Usage" %in% colnames(res)
  out <- res %>%
    group_by(Clade, AminoAcid) %>%
    summarise(
      Count = sum(Count, na.rm = TRUE),
      Codon_Usage = if (has_codon) dplyr::first(Codon_Usage) else NA_character_,
      .groups = "drop_last"
    ) %>%
    mutate(Total_in_Clade = sum(Count)) %>%
    mutate(`Frequency(%)` = (Count / Total_in_Clade) * 100) %>%
    ungroup()

  if (group_by == "Year" && isTRUE(hide_empty_years)) {
    out <- out %>% filter(Total_in_Clade > 0)
  }
  if (all(is.na(out$Codon_Usage))) out$Codon_Usage <- NULL
  out
}

usage_entropy_data <- function(subtype, var_type, gene, group_by, clade = "All") {
  clade_filter <- ""
  params <- list(subtype, var_type, gene, group_by)
  if (!identical(clade, "All")) {
    clade_filter <- " AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END = ?"
    params <- c(params, list(clade))
  }

  res <- usage_query(
    paste0(
      "SELECT Position, AminoAcid, SUM(Count) AS AA_Sum
       FROM usage
       WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?",
      clade_filter,
      " AND AminoAcid NOT IN ('X', '-')
       GROUP BY Position, AminoAcid"
    ),
    params
  )
  if (is.null(res)) return(NULL)

  res %>%
    group_by(Position) %>%
    mutate(Pos_Total = sum(AA_Sum), p = AA_Sum / Pos_Total) %>%
    filter(p > 0) %>%
    summarise(Entropy = -sum(p * log2(p)), .groups = "drop")
}

usage_lollipop_consensus <- function(subtype, var_type, gene, group_by, ref_group, tar_group, min_freq) {
  res <- usage_query(
    "WITH agg AS (
       SELECT Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
         AminoAcid, SUM(Count) AS Count
       FROM usage
       WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
         AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END IN (?, ?)
         AND AminoAcid NOT IN ('X', '-')
       GROUP BY Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
         AminoAcid
     ),
     freq AS (
       SELECT *,
         SUM(Count) OVER (PARTITION BY Position, Clade) AS Valid_Total,
         100.0 * Count / SUM(Count) OVER (PARTITION BY Position, Clade) AS New_Frequency
       FROM agg
     ),
     ranked AS (
       SELECT *,
         ROW_NUMBER() OVER (PARTITION BY Position, Clade ORDER BY New_Frequency DESC, AminoAcid) AS rn
       FROM freq
     )
     SELECT Position, Clade, AminoAcid, New_Frequency
     FROM ranked
     WHERE rn = 1 AND New_Frequency >= ?",
    list(subtype, var_type, gene, group_by, ref_group, tar_group, as.numeric(min_freq))
  )
  if (is.null(res)) return(NULL)
  res
}
