# ==========================================
# 3. SERVER LOGIC
# ==========================================
server <- function(input, output, session) {

  scalar_input <- function(x) {
    if (is.null(x) || length(x) == 0) return(NULL)
    if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
    x <- as.character(x)
    if (length(x) == 0 || is.na(x[[1]]) || identical(x[[1]], "")) return(NULL)
    x[[1]]
  }

  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x[[1]]) || identical(x[[1]], "")) y else x
  }

  selectize_state <- new.env(parent = emptyenv())

  update_selectize_if_changed <- function(input_id, choices, selected = NULL, server = TRUE) {
    choice_values <- as.character(choices)
    choice_names <- names(choices)
    choice_key <- if (!is.null(choice_names)) {
      paste(choice_names, choice_values, sep = "\t")
    } else {
      choice_values
    }
    selected_key <- scalar_input(selected) %||% "__NULL__"
    key <- paste(c(choice_key, selected_key), collapse = "\r")
    if (identical(selectize_state[[input_id]], key)) return(invisible(FALSE))
    selectize_state[[input_id]] <- key
    freezeReactiveValue(input, input_id)
    updateSelectizeInput(session, input_id, choices = choices, selected = selected, server = server)
    invisible(TRUE)
  }
  
  # --- REACTIVE DATA SWITCH ---
  
  current_colors <- reactive({
    if(input$variation_type == "AA") aa_colors else nt_colors
  })
  variant_label <- reactive({
    if(input$variation_type == "AA") "AA" else "NT"
  })

  # --- MEMORY MONITOR & CONTROL ---
  mem_timer <- reactiveTimer(5000) # Update every 5 seconds
  output$mem_usage <- renderText({
    mem_timer()
    # gc()[, 2] returns the used memory in MB for Ncells and Vcells
    mem_mb <- sum(gc()[, 2])
    paste0(round(mem_mb, 1), " MB")
  })
  
  observeEvent(input$free_mem, {
    # Clear global LRU cache
    lazy_cache$keys <- character(0)
    lazy_cache$data <- list()
    
    # Force aggressive garbage collection
    gc(verbose = FALSE)
    
    showNotification("Global data cache cleared and memory freed.", type = "message", duration = 3)
  })
  
  # Clean up session-specific memory when user disconnects
  session$onSessionEnded(function() {
    sp_data_val(NULL)
    pw_diff_val(NULL)
    gc(verbose = FALSE)
  })

  # --- GLOBAL LOADING INDICATOR FOR CONTEXT SWITCH ---
  observeEvent(list(input$global_subtype, input$variation_type), {
    waiter_show(
      html = tagList(
        spin_fading_circles(),
        tags$h3("Updating Dataset...", style = "color:white; margin-top: 20px;"),
        tags$p("Please wait while we update the application context.", style = "color:white;")
      ),
      color = "rgba(44, 62, 80, 0.9)"
    )
    
    session$onFlushed(function() {
      waiter_hide()
    }, once = TRUE)
  }, ignoreInit = TRUE, priority = 100)

  clade_explorer_available <- function() {
    !is.null(metadata_clade_explorer$summaries) && nrow(metadata_clade_explorer$summaries) > 0
  }

  clade_annotation_label <- function(annotation) {
    label <- gsub("_", " ", annotation)
    label <- gsub("^clade$", "HA clade", label, ignore.case = TRUE)
    label <- gsub("^G clade$", "NA clade", label, ignore.case = TRUE)
    label
  }

  available_genes <- reactive({
    req(input$global_subtype, input$variation_type)
    genes <- available_count_genes(input$global_subtype, input$variation_type, prefer_cache = TRUE)
    sort(genes)
  })
  
  current_usage_by_clade <- reactive({
    req(input$global_subtype, input$variation_type, input$sp_gene)
    var_lower <- tolower(input$variation_type)

    if (usage_duckdb_available()) {
      groups <- usage_available_groups(input$global_subtype, input$variation_type, input$sp_gene)
      group_by <- if ("HA_clade" %in% groups) "HA_clade" else setdiff(groups, c("Year", "Year_Month"))[1]
      if (is.na(group_by) || length(group_by) == 0) group_by <- groups[1]

      res <- usage_query(
        "SELECT DISTINCT \"Group\", Gene, Clade, Position
         FROM usage
         WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?",
        list(input$global_subtype, input$variation_type, input$sp_gene, group_by)
      )
      if (!is.null(res) && nrow(res) > 0) return(res)
    }
    
    dir_path <- count_cache_gene_path(input$global_subtype, input$variation_type, input$sp_gene)
    rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, input$sp_gene, "HA_clade")
    
    # Fallback to the first available group file if HA_clade is not available
    if (!file.exists(rds_file)) {
      files <- list.files(dir_path, pattern = paste0("^", var_lower, "_usage_by_.*\\.rds$"))
      # Exclude Year/Month files for the 'by_clade' fallback
      files <- setdiff(files, c(paste0(var_lower, "_usage_by_Year.rds"), paste0(var_lower, "_usage_by_Year_Month.rds")))
      if (length(files) > 0) {
        rds_file <- file.path(dir_path, files[1])
      }
    }
    
    df <- get_lazy_table(rds_file)
    if (!is.null(df)) {
      # Try to find a suitable column to rename to Clade
      possible_cols <- c("HA_clade", "NA_clade", "clade", "G_clade")
      found_col <- intersect(possible_cols, colnames(df))
      
      if (length(found_col) > 0) {
        df <- df %>% dplyr::rename(Clade = !!sym(found_col[1]))
      } else {
        # If no clade column found, and it's not a Year/Month table, 
        # try to rename the first column that matches the filename grouping
        group_col <- sub(paste0("^", var_lower, "_usage_by_(.*)\\.rds$"), "\\1", basename(rds_file))
        if (group_col %in% colnames(df)) {
          df <- df %>% dplyr::rename(Clade = !!sym(group_col))
        } else if (!"Clade" %in% colnames(df)) {
          # Last resort: if Clade still missing, create it as 'Unknown'
          df$Clade <- "Unknown"
        }
      }
      return(df)
    }
    return(data.frame(Group=character(), Gene=character(), Clade=character(), Position=numeric(), AminoAcid=character(), Count=numeric()))
  })

  # --- HELPER: Update Grouping Choices based on Loaded Data ---
  observe({
    req(input$global_subtype, input$variation_type, input$sp_gene)
    available_groups <- usage_available_groups(input$global_subtype, input$variation_type, input$sp_gene)
    
    if (length(available_groups) == 0) return()

    group_map <- setNames(available_groups, available_groups)
    
    if ("Year_Month" %in% names(group_map)) names(group_map)[group_map == "Year_Month"] <- "Year-Month"
    other_indices <- which(!(group_map %in% c("Year", "Year_Month")))
    for (i in other_indices) {
      names(group_map)[i] <- gsub("_", " ", group_map[i])
    }
    
    priority_keys <- intersect(c("Year", "Year_Month"), available_groups)
    other_keys <- sort(setdiff(available_groups, priority_keys))
    final_ordered_keys <- c(priority_keys, other_keys)
    
    final_choices <- setNames(final_ordered_keys, names(group_map)[match(final_ordered_keys, group_map)])
    
    # Robustly determine the next selection to avoid reactive loops
    current_sel <- input$sp_group_by
    if (is.null(current_sel) || !(current_sel %in% available_groups)) {
      if ("Year" %in% available_groups) {
        current_sel <- "Year"
      } else if (length(available_groups) > 0) {
        current_sel <- available_groups[1]
      } else {
        current_sel <- NULL # No valid selection
      }
    }
    updateSelectInput(session, "sp_group_by", choices = final_choices, selected = current_sel)
  })

  add_year_month_filter_column <- function(data) {
    if (is.null(data) || !is.data.frame(data)) return(data)

    year_chr <- if ("Year" %in% colnames(data)) trimws(as.character(data$Year)) else rep(NA_character_, nrow(data))
    month_chr <- if ("Month" %in% colnames(data)) trimws(as.character(data$Month)) else rep(NA_character_, nrow(data))
    month_chr <- ifelse(grepl("^[0-9]+$", month_chr), stringr::str_pad(month_chr, width = 2, side = "left", pad = "0"), month_chr)

    data$Year_Month_Filter <- dplyr::case_when(
      year_chr %in% c("", "NA") | month_chr %in% c("", "NA") ~ NA_character_,
      year_chr %in% c("Unknown", "unassigned", "Unassigned") ~ year_chr,
      month_chr %in% c("Unknown", "unassigned", "Unassigned") ~ month_chr,
      TRUE ~ paste0(year_chr, "-", month_chr)
    )

    data
  }

  sp_year_month_choices <- reactive({
    req(input$global_subtype, input$variation_type, input$sp_gene, input$sp_group_by, sp_position_debounced())
    pos <- sp_position_debounced()

    if (usage_duckdb_available()) {
      return(usage_year_month_choices(input$global_subtype, input$variation_type, input$sp_gene, input$sp_group_by, pos))
    }

    rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, input$sp_gene, input$sp_group_by)
    data <- get_lazy_table(rds_file)

    if (is.null(data) || !all(c("Year", "Month") %in% colnames(data))) return(character(0))

    data <- add_year_month_filter_column(data)

    ym_values <- data %>%
      filter(Group == input$global_subtype, Gene == input$sp_gene, Position == pos) %>%
      pull(Year_Month_Filter) %>%
      unique() %>%
      stats::na.omit() %>%
      as.character()

    special_values <- c("Unknown", "unassigned", "Unassigned")
    present_specials <- intersect(special_values, ym_values)
    chronological_yms <- sort(setdiff(ym_values, special_values))
    c(present_specials, chronological_yms)
  })

  get_selected_year_months <- function(ym_values, start_value, end_value) {
    if (length(ym_values) == 0 || is.null(start_value) || is.null(end_value)) return(NULL)
    if (!(start_value %in% ym_values) || !(end_value %in% ym_values)) return(NULL)

    selected_idx <- match(c(start_value, end_value), ym_values)
    ym_values[seq(min(selected_idx), max(selected_idx))]
  }

  every_nth_value <- function(values, n = 6) {
    values <- as.character(values)
    if (length(values) == 0) return(character(0))
    values[seq.int(1, length(values), by = n)]
  }

  output$sp_year_month_range_ui <- renderUI({
    ym_values <- sp_year_month_choices()

    if (length(ym_values) == 0) {
      return(
        helpText("Time range slider is unavailable for this grouping.")
      )
    }

    start_sel <- if (!is.null(input$sp_year_month_start) && input$sp_year_month_start %in% ym_values) {
      input$sp_year_month_start
    } else {
      ym_values[1]
    }
    end_sel <- if (!is.null(input$sp_year_month_end) && input$sp_year_month_end %in% ym_values) {
      input$sp_year_month_end
    } else {
      ym_values[length(ym_values)]
    }

    tagList(
      tags$label("Filter Year-Month:", style = "font-weight: bold; color: #2c3e50;"),
      fluidRow(
        column(6, selectInput("sp_year_month_start", "Start:", choices = ym_values, selected = start_sel, width = "100%")),
        column(6, selectInput("sp_year_month_end", "End:", choices = ym_values, selected = end_sel, width = "100%"))
      )
    )
  })
  
  # --- HELPER: Update Gene Dropdowns based on Subtype ---
  observeEvent(available_genes(), {
    genes <- available_genes()
    
    # Smart mapping for gene selection when switching between AA and NT
    # Function to pick the best matching gene from the new list
    get_best_gene <- function(current_gene, available_list) {
      if (is.null(current_gene) || current_gene == "") return(if ("HA" %in% available_list) "HA" else available_list[1])
      if (current_gene %in% available_list) return(current_gene)
      
      # Mapping AA -> NT
      if (current_gene %in% c("M1", "M2") && "M" %in% available_list) return("M")
      if (current_gene %in% c("NS1", "NEP") && "NS" %in% available_list) return("NS")
      
      # Mapping NT -> AA
      if (current_gene == "M") {
        if ("M1" %in% available_list) return("M1")
        if ("M2" %in% available_list) return("M2")
      }
      if (current_gene == "NS") {
        if ("NS1" %in% available_list) return("NS1")
        if ("NEP" %in% available_list) return("NEP")
      }
      
      # Default fallback
      if ("HA" %in% available_list) return("HA")
      return(available_list[1])
    }
    
    sel_sp   <- get_best_gene(input$sp_gene, genes)
    sel_ent  <- get_best_gene(input$ent_gene, genes)
    sel_lol  <- get_best_gene(input$lol_gene, genes)
    # sel_heat <- get_best_gene(input$heat_gene, genes)
    
    updateSelectInput(session, "sp_gene", choices = genes, selected = sel_sp)
    updateSelectInput(session, "ent_gene", choices = genes, selected = sel_ent)
    updateSelectInput(session, "lol_gene", choices = genes, selected = sel_lol)
    # updateSelectInput(session, "heat_gene", choices = genes, selected = sel_heat)
  })
  
  # --- HELPER: Update Pairwise & Landscape Grouping Choices ---
  observeEvent(list(input$global_subtype, input$variation_type, input$sp_gene), {
    req(input$sp_gene)
    available_groups <- usage_available_groups(input$global_subtype, input$variation_type, input$sp_gene)
    
    if (length(available_groups) == 0) return()
    
    group_map <- setNames(available_groups, available_groups)
    if ("Year_Month" %in% names(group_map)) names(group_map)[group_map == "Year_Month"] <- "Year-Month"
    for (i in which(!(group_map %in% c("Year", "Year_Month")))) {
      names(group_map)[i] <- gsub("_", " ", group_map[i])
    }
    
    final_choices <- setNames(available_groups, names(group_map))
    
    current_sel <- if (!is.null(input$sp_group_by) && input$sp_group_by %in% available_groups) {
      input$sp_group_by
    } else if (!is.null(input$pw_group_by) && input$pw_group_by %in% available_groups) {
      input$pw_group_by
    } else {
      if ("HA_clade" %in% available_groups) "HA_clade" else available_groups[1]
    }
    updateSelectInput(session, "pw_group_by", choices = final_choices, selected = current_sel)
    updateSelectInput(session, "ent_group_by", choices = final_choices, selected = current_sel)
    updateSelectInput(session, "lol_group_by", choices = final_choices, selected = current_sel)
    # updateSelectInput(session, "heat_group_by", choices = final_choices, selected = current_sel)
  })

  empty_pairwise_usage_df <- function() {
    data.frame(
      Group = character(),
      Gene = character(),
      Clade = character(),
      Position = numeric(),
      AminoAcid = character(),
      Count = numeric(),
      stringsAsFactors = FALSE
    )
  }

  empty_pairwise_position_df <- function(include_codon = TRUE) {
    out <- data.frame(
      Clade = character(),
      AminoAcid = character(),
      Count = numeric(),
      Total_in_Clade = numeric(),
      `Frequency(%)` = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (include_codon) out$Codon_Usage <- character()
    out
  }

  get_pairwise_rds_file <- function(gene, group_by = input$pw_group_by) {
    count_cache_file_path(input$global_subtype, input$variation_type, gene, group_by)
  }

  load_pairwise_gene_data <- function(gene, group_by = input$pw_group_by) {
    req(input$global_subtype, input$variation_type, group_by)

    if (usage_duckdb_available()) {
      df <- usage_pairwise_gene_data(input$global_subtype, input$variation_type, gene, group_by)
      if (is.null(df)) return(empty_pairwise_usage_df())
      return(df)
    }

    df <- get_lazy_table(get_pairwise_rds_file(gene, group_by), max_tables = 2, max_mem_mb = 450)
    if (is.null(df)) return(empty_pairwise_usage_df())

    if (group_by %in% colnames(df)) {
      df <- df %>% dplyr::rename(Clade = !!sym(group_by))
    } else if (!"Clade" %in% colnames(df)) {
      df$Clade <- "Unknown"
    }

    df %>%
      filter(Group == input$global_subtype, Gene == gene)
  }

  get_dominant_variants_for_clade <- function(gene_data, clade_name, min_freq) {
    gene_data %>%
      filter(Clade == clade_name) %>%
      filter(!(AminoAcid %in% c("X", "-"))) %>%
      group_by(Gene, Position, AminoAcid) %>%
      summarise(Variant_Count = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
      mutate(
        Total_Seqs = sum(Variant_Count),
        Freq = (Variant_Count / Total_Seqs) * 100
      ) %>%
      filter(Freq == max(Freq)) %>%
      filter(row_number() == 1, Freq >= min_freq) %>%
      ungroup()
  }

  get_pairwise_position_distribution <- function(gene, position, group_by = input$pw_group_by) {
    req(position)

    if (usage_duckdb_available()) {
      res <- usage_position_distribution(input$global_subtype, input$variation_type, gene, group_by, position, input$pw_hide_empty_years)
      if (is.null(res)) return(empty_pairwise_position_df(include_codon = FALSE))
      return(res)
    }

    res <- load_pairwise_gene_data(gene, group_by) %>%
      filter(Position == position) %>%
      filter(!(AminoAcid %in% c("X", "-")))

    if (nrow(res) == 0) return(empty_pairwise_position_df(include_codon = "Codon_Usage" %in% colnames(res)))

    has_codon <- "Codon_Usage" %in% colnames(res)

    res <- res %>%
      group_by(Clade, AminoAcid) %>%
      summarise(
        Count = sum(Count, na.rm = TRUE),
        Codon_Usage = if (has_codon) dplyr::first(Codon_Usage) else NA_character_,
        .groups = "drop_last"
      ) %>%
      mutate(Total_in_Clade = sum(Count)) %>%
      mutate(`Frequency(%)` = (Count / Total_in_Clade) * 100) %>%
      ungroup()

    if (group_by == "Year" && isTRUE(input$pw_hide_empty_years)) {
      res <- res %>% filter(Total_in_Clade > 0)
    }

    res
  }

  compute_pairwise_differences <- function(clade1, clade2, min_freq) {
    genes <- available_genes()
    diff_list <- vector("list", length(genes))
    diff_idx <- 0

    message("Pairwise comparison started. Memory before loop: ", round(sum(gc(verbose = FALSE)[, 2]), 1), " MB")

    for (gene in genes) {
      if (usage_duckdb_available()) {
        gene_diffs <- usage_pairwise_differences_for_gene(
          input$global_subtype,
          input$variation_type,
          gene,
          input$pw_group_by,
          clade1,
          clade2,
          min_freq
        )

        if (!is.null(gene_diffs) && nrow(gene_diffs) > 0) {
          diff_idx <- diff_idx + 1
          diff_list[[diff_idx]] <- gene_diffs
        }

        rm(gene_diffs)
        gc(FALSE)
        next
      }

      gene_data <- load_pairwise_gene_data(gene) %>%
        filter(Clade %in% c(clade1, clade2))

      if (nrow(gene_data) == 0) {
        rm(gene_data)
        gc(FALSE)
        next
      }

      c1_dom <- get_dominant_variants_for_clade(gene_data, clade1, min_freq) %>%
        dplyr::select(Gene, Position, Clade1_AA = AminoAcid, Clade1_Freq = Freq)

      c2_dom <- get_dominant_variants_for_clade(gene_data, clade2, min_freq) %>%
        dplyr::select(Gene, Position, Clade2_AA = AminoAcid, Clade2_Freq = Freq)

      gene_diffs <- inner_join(c1_dom, c2_dom, by = c("Gene", "Position")) %>%
        filter(Clade1_AA != Clade2_AA)

      if (nrow(gene_diffs) > 0) {
        diff_idx <- diff_idx + 1
        diff_list[[diff_idx]] <- gene_diffs
      }

      rm(gene_data, c1_dom, c2_dom, gene_diffs)
      gc(FALSE)
    }

    out <- if (diff_idx == 0) {
      data.frame(
        Gene = character(),
        Position = numeric(),
        Clade1_AA = character(),
        Clade1_Freq = numeric(),
        Clade2_AA = character(),
        Clade2_Freq = numeric(),
        stringsAsFactors = FALSE
      )
    } else {
      bind_rows(diff_list[seq_len(diff_idx)]) %>%
        arrange(Gene, Position)
    }

    message("Pairwise comparison finished. Memory after loop: ", round(sum(gc(verbose = FALSE)[, 2]), 1), " MB")
    out
  }

  ent_usage_data <- reactive({
    req(input$global_subtype, input$variation_type, input$ent_group_by, input$ent_gene)
    if (usage_duckdb_available()) {
      df <- usage_pairwise_gene_data(input$global_subtype, input$variation_type, input$ent_gene, input$ent_group_by)
      if (is.null(df)) {
        return(data.frame(Group=character(), Gene=character(), Clade=character(), Position=numeric(), AminoAcid=character(), Count=numeric()))
      }
      return(df)
    }

    var_lower <- tolower(input$variation_type)
    rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, input$ent_gene, input$ent_group_by)
    
    df <- get_lazy_table(rds_file)
    if (!is.null(df)) {
      if (input$ent_group_by %in% colnames(df)) {
        df <- df %>% dplyr::rename(Clade = !!sym(input$ent_group_by))
      } else if (!"Clade" %in% colnames(df)) {
        df$Clade <- "Unknown"
      }
      return(df)
    }
    return(data.frame(Group=character(), Gene=character(), Clade=character(), Position=numeric(), AminoAcid=character(), Count=numeric()))
  })

  lol_usage_data <- reactive({
    req(input$global_subtype, input$variation_type, input$lol_group_by, input$lol_gene)
    if (usage_duckdb_available()) {
      df <- usage_pairwise_gene_data(input$global_subtype, input$variation_type, input$lol_gene, input$lol_group_by)
      if (is.null(df)) {
        return(data.frame(Group=character(), Gene=character(), Clade=character(), Position=numeric(), AminoAcid=character(), Count=numeric()))
      }
      return(df)
    }

    var_lower <- tolower(input$variation_type)
    rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, input$lol_gene, input$lol_group_by)
    
    df <- get_lazy_table(rds_file)
    if (!is.null(df)) {
      if (input$lol_group_by %in% colnames(df)) {
        df <- df %>% dplyr::rename(Clade = !!sym(input$lol_group_by))
      } else if (!"Clade" %in% colnames(df)) {
        df$Clade <- "Unknown"
      }
      return(df)
    }
    return(data.frame(Group=character(), Gene=character(), Clade=character(), Position=numeric(), AminoAcid=character(), Count=numeric()))
  })

  # heat_usage_data <- reactive({
  #   req(input$global_subtype, input$variation_type, input$heat_group_by, input$heat_gene)
  #   var_lower <- tolower(input$variation_type)
  #   rds_file <- paste0("data/", input$global_subtype, "/", input$variation_type, "/", input$heat_gene, "/", var_lower, "_usage_by_", input$heat_group_by, ".rds")
  #   
  #   df <- get_lazy_table(rds_file)
  #   if (!is.null(df)) {
  #     if (input$heat_group_by %in% colnames(df)) {
  #       df <- df %>% dplyr::rename(Clade = !!sym(input$heat_group_by))
  #     } else if (!"Clade" %in% colnames(df)) {
  #       df$Clade <- "Unknown"
  #     }
  #     return(df)
  #   }
  #   return(data.frame(Group=character(), Gene=character(), Clade=character(), Position=numeric(), AminoAcid=character(), Count=numeric()))
  # })

  # --- HELPER: Update Clade Dropdowns based on Subtype ---
  observeEvent(list(input$global_subtype, input$variation_type, input$pw_group_by), {
    req(input$global_subtype, input$variation_type, input$pw_group_by)
    
    # Fast path: Check one gene instead of evaluating pairwise_usage_data() which eagerly loads all genes
    genes <- available_count_genes(input$global_subtype, input$variation_type, prefer_cache = TRUE)
    if(length(genes) == 0) return()
    gene <- if("HA" %in% genes) "HA" else genes[1]

    if (usage_duckdb_available()) {
      clades <- usage_distinct_group_values(input$global_subtype, input$variation_type, gene, input$pw_group_by)
      if (length(clades) == 0) return()
      updateSelectInput(session, "pw_clade1", choices = c("Select Group..." = "", clades), selected = "")
      updateSelectInput(session, "pw_clade2", choices = c("Select Group..." = "", clades), selected = "")
      return()
    }
    
    rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, gene, input$pw_group_by)
    df <- get_lazy_table(rds_file)
    if(is.null(df)) return()
    
    clades <- if(input$pw_group_by %in% colnames(df)) unique(as.character(df[[input$pw_group_by]])) else if("Clade" %in% colnames(df)) unique(as.character(df$Clade)) else "Unknown"
    clades <- sort(clades)
    
    special_values <- c("Unknown", "unassigned", "Unassigned")
    present_specials <- intersect(special_values, clades)
    if (length(present_specials) > 0) {
      clades <- c(setdiff(clades, present_specials), present_specials)
    }
    
    clades_choices <- c("Select Group..." = "", clades)
    
    updateSelectInput(session, "pw_clade1", choices = clades_choices, selected = "")
    updateSelectInput(session, "pw_clade2", choices = clades_choices, selected = "")
  }, ignoreInit = TRUE)
  
  observeEvent(list(input$global_subtype, input$variation_type, input$ent_group_by, input$ent_gene), {
    req(input$global_subtype, input$variation_type, input$ent_group_by, input$ent_gene)

    if (usage_duckdb_available()) {
      clade_choices <- usage_distinct_group_values(input$global_subtype, input$variation_type, input$ent_gene, input$ent_group_by)
      if (length(clade_choices) == 0) return()
      updateSelectInput(session, "ent_group", choices = c("All", clade_choices), selected = "All")
      return()
    }
    
    rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, input$ent_gene, input$ent_group_by)
    df <- get_lazy_table(rds_file)
    if(is.null(df)) return()
    
    clades <- if(input$ent_group_by %in% colnames(df)) unique(as.character(df[[input$ent_group_by]])) else if("Clade" %in% colnames(df)) unique(as.character(df$Clade)) else "Unknown"
    clade_choices <- sort(clades)
    
    updateSelectInput(session, "ent_group", choices = c("All", clade_choices), selected = "All")
  }, ignoreInit = TRUE)
  
  observeEvent(list(input$global_subtype, input$variation_type, input$lol_group_by, input$lol_gene), {
    req(input$global_subtype, input$variation_type, input$lol_group_by, input$lol_gene)

    if (usage_duckdb_available()) {
      clades <- usage_distinct_group_values(input$global_subtype, input$variation_type, input$lol_gene, input$lol_group_by)
      if (length(clades) == 0) return()
      updateSelectInput(session, "lol_ref_group", choices = clades, selected = clades[1])
      updateSelectInput(session, "lol_tar_group", choices = clades, selected = if(length(clades)>1) clades[2] else clades[1])
      return()
    }
    
    rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, input$lol_gene, input$lol_group_by)
    df <- get_lazy_table(rds_file)
    if(is.null(df)) return()
    
    clades <- if(input$lol_group_by %in% colnames(df)) unique(as.character(df[[input$lol_group_by]])) else if("Clade" %in% colnames(df)) unique(as.character(df$Clade)) else "Unknown"
    clades <- sort(clades)
    
    special_values <- c("Unknown", "unassigned", "Unassigned")
    present_specials <- intersect(special_values, clades)
    if (length(present_specials) > 0) {
      clades <- c(setdiff(clades, present_specials), present_specials)
    }
    updateSelectInput(session, "lol_ref_group", choices = clades, selected = clades[1])
    updateSelectInput(session, "lol_tar_group", choices = clades, selected = if(length(clades)>1) clades[2] else if(length(clades)>0) clades[1] else NULL)
  }, ignoreInit = TRUE)
  
  # --- HELPER: Dynamically Update Position Limit based on Gene Length (Tab 1) ---
  observeEvent(list(input$sp_gene, input$variation_type), {
    req(input$global_subtype, input$sp_gene)
    if (usage_duckdb_available()) {
      max_val <- usage_max_position(input$global_subtype, input$variation_type, input$sp_gene)
      if (!is.na(max_val)) updateNumericInput(session, "sp_position", max = max_val)
      return()
    }

    g_data <- current_usage_by_clade() %>% filter(Group == input$global_subtype, Gene == input$sp_gene)
    if(nrow(g_data) > 0) updateNumericInput(session, "sp_position", max = max(g_data$Position, na.rm=TRUE))
  })
  
  # --- HELPER: Disable/Hide Quick Access in NT mode ---
  observeEvent(input$variation_type, {
    if (input$variation_type == "NT") {
      shinyjs::hide("sp_quick_access_section")
      updateSelectInput(session, "sp_quick_visit", selected = "None")
    } else {
      shinyjs::show("sp_quick_access_section")
    }
  })

  # ==========================================
  # SERVER: GENETIC CLADE
  # ==========================================

  output$gc_status_notice <- renderUI({
    if (clade_explorer_available()) return(NULL)
    div(
      class = "alert alert-warning",
      tags$strong("Genetic Clade summary is not available. "),
      "This view requires a refreshed metadata cache built from the subtype metadata CSV files. ",
      "The rest of FLUExplorer can still run from the existing cache data."
    )
  })

  gc_annotation_choices <- reactive({
    if (!clade_explorer_available()) return(character(0))
    req(input$global_subtype)
    annotations <- metadata_clade_explorer$summaries %>%
      filter(.data$Group == input$global_subtype) %>%
      pull(.data$Annotation) %>%
      unique() %>%
      sort()
    stats::setNames(annotations, vapply(annotations, clade_annotation_label, character(1)))
  })

  observeEvent(list(input$global_subtype, gc_annotation_choices()), {
    choices <- gc_annotation_choices()
    if (length(choices) == 0) {
      updateSelectInput(session, "gc_annotation", choices = character(0), selected = NULL)
      update_selectize_if_changed("gc_clade", choices = character(0), selected = NULL, server = FALSE)
      return(invisible(NULL))
    }
    current <- scalar_input(input$gc_annotation)
    selected <- if (!is.null(current) && current %in% unname(choices)) current else unname(choices)[[1]]
    updateSelectInput(session, "gc_annotation", choices = choices, selected = selected)
  }, ignoreInit = FALSE, priority = 100)

  observeEvent(list(input$global_subtype, input$gc_annotation), {
    if (!clade_explorer_available()) {
      update_selectize_if_changed("gc_clade", choices = character(0), selected = NULL, server = FALSE)
      return(invisible(NULL))
    }

    req(input$global_subtype)
    annotation <- scalar_input(input$gc_annotation)
    if (is.null(annotation)) {
      choices <- gc_annotation_choices()
      annotation <- if (length(choices) > 0) unname(choices)[[1]] else NULL
    }
    req(annotation)

    choices_df <- metadata_clade_explorer$summaries %>%
      filter(.data$Group == input$global_subtype, .data$Annotation == annotation) %>%
      arrange(.data$Rank, .data$Clade)

    if (nrow(choices_df) == 0) {
      update_selectize_if_changed("gc_clade", choices = character(0), selected = NULL, server = FALSE)
      return(invisible(NULL))
    }

    labels <- paste0(
      choices_df$Clade,
      " | rank #", choices_df$Rank,
      " | n=", scales::comma(choices_df$StrainCount)
    )
    choices <- stats::setNames(choices_df$Clade, labels)
    current_clade <- scalar_input(input$gc_clade)
    selected <- if (!is.null(current_clade) && current_clade %in% choices_df$Clade) {
      current_clade
    } else {
      choices_df$Clade[[1]]
    }

    update_selectize_if_changed("gc_clade", choices = choices, selected = selected, server = FALSE)
  }, ignoreInit = FALSE, priority = 100)

  gc_selected_summary <- reactive({
    validate(need(clade_explorer_available(), "Refresh the metadata summary cache to use Genetic Clade."))
    req(input$global_subtype)
    annotation <- scalar_input(input$gc_annotation)
    clade <- scalar_input(input$gc_clade)
    req(annotation, clade)
    out <- metadata_clade_explorer$summaries %>%
      filter(.data$Group == input$global_subtype, .data$Annotation == annotation, .data$Clade == clade)
    validate(need(nrow(out) > 0, "Choose a clade to explore."))
    out %>% slice_head(n = 1)
  })

  gc_timeseries <- reactive({
    summary_row <- gc_selected_summary()
    selected_months <- metadata_clade_explorer$monthly %>%
      filter(
        .data$Group == summary_row$Group[[1]],
        .data$Annotation == summary_row$Annotation[[1]],
        .data$Clade == summary_row$Clade[[1]]
      ) %>%
      select(YearMonth, Count)

    metadata_clade_explorer$month_totals %>%
      filter(.data$Group == summary_row$Group[[1]]) %>%
      left_join(selected_months, by = "YearMonth") %>%
      mutate(
        Count = dplyr::coalesce(.data$Count, 0L),
        Percent = dplyr::if_else(.data$Total > 0, (.data$Count / .data$Total) * 100, 0)
      ) %>%
      arrange(.data$YearMonth)
  })

  output$gc_summary_cards <- renderUI({
    summary_row <- gc_selected_summary()
    period <- if (is.na(summary_row$FirstMonth[[1]]) || is.na(summary_row$LastMonth[[1]])) {
      "Not available"
    } else {
      paste(summary_row$FirstMonth[[1]], "to", summary_row$LastMonth[[1]])
    }
    peak <- if (is.na(summary_row$PeakMonth[[1]])) {
      "Not available"
    } else {
      paste0(summary_row$PeakMonth[[1]], " (", round(summary_row$PeakPercent[[1]], 2), "%)")
    }

    fluidRow(
      column(
        3,
        div(class = "summary-card gc-summary-card",
            div(class = "summary-card-title", "Total Sequences"),
            div(class = "summary-card-value", scales::comma(summary_row$StrainCount[[1]])),
            tags$p(class = "gc-summary-card-note", clade_annotation_label(summary_row$Annotation[[1]])))
      ),
      column(
        3,
        div(class = "summary-card gc-summary-card",
            div(class = "summary-card-title", "Rank / Share"),
            div(class = "summary-card-value", paste0("#", summary_row$Rank[[1]])),
            tags$p(class = "gc-summary-card-note", paste0(round(summary_row$DatasetShare[[1]], 2), "% of subtype")))
      ),
      column(
        3,
        div(class = "summary-card gc-summary-card",
            div(class = "summary-card-title", "Active Period"),
            div(class = "summary-card-value", period),
            tags$p(class = "gc-summary-card-note", "First to last month"))
      ),
      column(
        3,
        div(class = "summary-card gc-summary-card",
            div(class = "summary-card-title", "Peak Month"),
            div(class = "summary-card-value", peak),
            tags$p(class = "gc-summary-card-note", "Highest monthly prevalence"))
      )
    )
  })

  output$gc_prevalence_plot <- renderPlotly({
    summary_row <- gc_selected_summary()
    df <- gc_timeseries()
    validate(need(nrow(df) > 0, "No monthly prevalence data are available."))
    annotation_label <- clade_annotation_label(summary_row$Annotation[[1]])
    clade_label <- summary_row$Clade[[1]]

    hover_text <- paste0(
      "Month: ", df$YearMonth,
      "<br>", annotation_label, ": ", clade_label,
      "<br>Clade count: ", scales::comma(df$Count),
      "<br>Total monthly sequences: ", scales::comma(df$Total),
      "<br>Percent: ", round(df$Percent, 3), "%"
    )

    plot_ly() %>%
      add_bars(
        data = df,
        x = ~YearMonth,
        y = ~Count,
        name = "Clade count",
        yaxis = "y2",
        marker = list(color = "rgba(89, 163, 168, 0.28)", line = list(color = "rgba(89, 163, 168, 0.45)", width = 0.5)),
        text = hover_text,
        hoverinfo = "text"
      ) %>%
      add_trace(
        data = df,
        x = ~YearMonth,
        y = ~Percent,
        name = "Monthly prevalence",
        type = "scatter",
        mode = "lines+markers",
        fill = "tozeroy",
        fillcolor = "rgba(41, 128, 185, 0.12)",
        line = list(color = "#2980b9", width = 2.5),
        marker = list(color = "#2980b9", size = 6),
        text = hover_text,
        hoverinfo = "text"
      ) %>%
      layout(
        barmode = "overlay",
        xaxis = list(title = "Year-Month", rangeslider = list(visible = nrow(df) > 24)),
        yaxis = list(title = "Percent of monthly sequences", ticksuffix = "%", rangemode = "tozero"),
        yaxis2 = list(title = "Clade count", overlaying = "y", side = "right", showgrid = FALSE, rangemode = "tozero"),
        legend = list(orientation = "h", x = 0, y = 1.12),
        margin = list(l = 70, r = 75, b = 80, t = 40)
      ) %>%
      config(displayModeBar = FALSE)
  })

  gc_breakdown_data <- function(category) {
    summary_row <- gc_selected_summary()
    metadata_clade_explorer$breakdowns %>%
      filter(
        .data$Group == summary_row$Group[[1]],
        .data$Annotation == summary_row$Annotation[[1]],
        .data$Clade == summary_row$Clade[[1]],
        .data$Category == category
      ) %>%
      arrange(.data$Count) %>%
      mutate(Value = factor(.data$Value, levels = .data$Value))
  }

  render_gc_breakdown_plot <- function(category, label, color) {
    renderPlotly({
      df <- gc_breakdown_data(category)
      validate(need(nrow(df) > 0, paste("No", tolower(label), "summary is available.")))
      plot_ly(
        df,
        x = ~Count,
        y = ~Value,
        type = "bar",
        orientation = "h",
        marker = list(color = color),
        text = ~paste0(label, ": ", Value, "<br>Count: ", scales::comma(Count)),
        hoverinfo = "text"
      ) %>%
        layout(
          xaxis = list(title = "Sequences"),
          yaxis = list(title = ""),
          margin = list(l = 110, r = 15, b = 45, t = 10)
        ) %>%
        config(displayModeBar = FALSE)
    })
  }

  output$gc_country_plot <- render_gc_breakdown_plot("country", "Country", "#4E79A7")
  output$gc_region_plot <- render_gc_breakdown_plot("region", "Region", "#59A3A8")
  output$gc_host_plot <- render_gc_breakdown_plot("host", "Host", "#F28E2B")

  output$gc_monthly_table <- renderDT({
    gc_timeseries() %>%
      transmute(
        `Year-Month` = YearMonth,
        `Clade Count` = Count,
        `Total Monthly Sequences` = Total,
        `Percent (%)` = round(Percent, 3)
      ) %>%
      datatable(options = list(pageLength = 12, autoWidth = TRUE), rownames = FALSE)
  })
  
  # ==========================================
  # SERVER: TAB 1 - STATS (Map & Static Plots)
  # ==========================================
  
  output$total_seqs <- renderText({ paste0(total_raw) })
  output$total_countries <- renderText({ paste0(total_countries_val) })
  output$time_range <- renderText({ time_range_val })

  metadata_plot_label <- function(col, subtype = NULL) {
    if (identical(col, "clade") && identical(subtype, "B_YAM")) return("Yamagata Clade")
    if (identical(col, "clade") && identical(subtype, "B_VIC")) return("Victoria Clade")
    if (identical(col, "clade")) return("HA Clade")
    if (identical(col, "G_clade")) return("NA Clade")
    label <- gsub("_", " ", col)
    label <- gsub("HA ", "HA ", label, fixed = TRUE)
    label <- gsub("NA ", "NA ", label, fixed = TRUE)
    label <- gsub("clade", "Clade", label, ignore.case = TRUE)
    label
  }

  # --- DYNAMIC CLADE PLOT FILL DROPDOWN ---
  observeEvent(input$clade_plot_subtype, {
    req(input$clade_plot_subtype)

    subtype_summary <- metadata_summary_stats %>%
      filter(.data$Group == input$clade_plot_subtype)

    informative_meta_cols <- metadata_grouping_cols[
      vapply(metadata_grouping_cols, function(col) {
        if (!col %in% names(subtype_summary)) return(FALSE)
        values <- unique(as.character(subtype_summary[[col]]))
        any(!is_unknown_metadata_value(values))
      }, logical(1))
    ]

    choices <- c(
      stats::setNames(
        informative_meta_cols,
        vapply(informative_meta_cols, metadata_plot_label, character(1), subtype = input$clade_plot_subtype)
      ),
      c("Region" = "region", "Country" = "country")
    )

    current_sel <- if (!is.null(input$clade_plot_fill) && input$clade_plot_fill %in% unname(choices)) {
      input$clade_plot_fill
    } else if ("clade" %in% unname(choices)) {
      "clade"
    } else {
      unname(choices)[1]
    }
    updateSelectInput(session, "clade_plot_fill", choices = choices, selected = current_sel)
  })

  # --- REACTIVE MAP DATA ---
  # map_data_filtered <- reactive({
  #   req(input$global_subtype, input$map_geo_level, input$map_clade_type, input$map_year)
  #   
  #   # PERFORMANCE: use pre-aggregated summary instead of full metadata_global
  #   plot_df <- metadata_summary_stats
  #   
  #   if(input$global_subtype != "All") {
  #     plot_df <- plot_df %>% filter(Group == input$global_subtype)
  #   }
  #   if(input$map_year != "All") {
  #     plot_df <- plot_df %>% filter(Year == input$map_year)
  #   }
  #   
  #   geo_col <- if(input$map_geo_level == "Region") "region" else "country"
  #   clade_col <- if(input$map_clade_type == "clade") "clade" else "G_clade"
  #   
  #   summary_df <- plot_df %>%
  #     filter(!!sym(clade_col) != "" & !is.na(!!sym(clade_col))) %>%
  #     group_by(!!sym(geo_col), !!sym(clade_col)) %>%
  #     summarise(n = sum(n), .groups = "drop") %>% # use sum(n) because it's pre-aggregated
  #     tidyr::pivot_wider(names_from = !!sym(clade_col), values_from = n, values_fill = 0)
  #   
  #   if(input$map_geo_level == "Region") {
  #     res <- inner_join(summary_df, region_coords, by = "region")
  #   } else {
  #     res <- inner_join(summary_df, world_coords, by = c("country" = "country"))
  #   }
  #   
  #   return(as.data.frame(res))
  # })
  # 
  # # --- RENDER MAP ---
  # output$world_map <- renderLeaflet({
  #   data <- map_data_filtered()
  #   validate(need(nrow(data) > 0, "No data available for the selected filters."))
  # 
  #   geo_col_name <- if(input$map_geo_level == "Region") "region" else "country"
  # 
  #   chart_cols <- sort(setdiff(colnames(data), c(geo_col_name, "lat", "lng")))
  # 
  #   active_colors <- if(input$map_clade_type == "clade") {
  #     as.character(clade_colors_vec[chart_cols])
  #   } else {
  #     as.character(g_clade_colors_vec[chart_cols])
  #   }
  # 
  #   leaflet(data) %>%
  #     addProviderTiles(providers$CartoDB.Positron) %>%
  #     setView(lng = 10, lat = 15, zoom = 2) %>%
  #     addMinicharts(
  #       data$lng, data$lat,
  #       type = "pie",
  #       chartdata = data[, chart_cols],
  #       colorPalette = active_colors,
  #       width = 45,
  #       transitionTime = 0,
  #       showLabels = FALSE
  #     )
  # })

  stats_metadata_filtered <- reactive({
    req(input$stats_year_range)
    
    major_continents <- c("Africa", "Asia", "Europe", "North America", "South America", "Oceania")
    
    metadata_summary_stats %>%
      filter(Year >= input$stats_year_range[1], Year <= input$stats_year_range[2]) %>%
      filter(region %in% major_continents)
  })

  output$stats_time_plot <- renderPlotly({
    plot_data <- stats_metadata_filtered() %>%
      group_by(Year, Group) %>%
      summarise(Count = sum(n), .groups = "drop")
      
    n_groups <- length(unique(plot_data$Group))
    my_colors <- setNames(viridis::viridis(n_groups, option = "turbo", begin = 0.1, end = 0.9), sort(unique(plot_data$Group)))
    
    plot_ly(plot_data, x = ~Year, y = ~Count, color = ~Group, colors = my_colors,
            type = "bar", hoverinfo = "text",
            text = ~paste0("Year: ", Year, "<br>Subtype: ", Group, "<br>Sequences: ", scales::comma(Count)),
            marker = list(line = list(color = 'white', width = 0.5))) %>%
      layout(barmode = 'stack',
             xaxis = list(title = "Year", tickangle = -45, tickfont = list(family = "Arial", size = 12), tickformat = "d"),
             yaxis = list(title = "Sequence Count", tickformat = ","),
             legend = list(orientation = 'h', x = 0.5, xanchor = 'center', y = -0.2),
             margin = list(b = 50)) %>%
      config(displayModeBar = FALSE)
  })

  output$stats_geo_plot <- renderPlotly({
    plot_data <- stats_metadata_filtered() %>%
      group_by(region) %>%
      summarise(Count = sum(n), .groups = "drop")
      
    n_regions <- length(unique(plot_data$region))
    my_colors <- setNames(viridis::viridis(n_regions, option = "mako"), sort(unique(plot_data$region)))
    
    plot_ly(plot_data, x = ~reorder(region, Count), y = ~Count, color = ~region, colors = my_colors,
            type = "bar", hoverinfo = "text",
            text = ~paste0("Region: ", region, "<br>Count: ", scales::comma(Count))) %>%
      layout(showlegend = FALSE,
             xaxis = list(title = "Region", tickfont = list(family = "Arial", size = 12)),
             yaxis = list(title = "Count", tickformat = ",")) %>%
      config(displayModeBar = FALSE)
  })

  output$stats_clade_plot <- renderPlotly({
    req(input$clade_plot_fill, input$clade_plot_subtype, input$clade_plot_palette)

    plot_df <- stats_metadata_filtered()
    # No longer check for "All" because it's removed from UI
    plot_df <- plot_df %>% filter(Group == input$clade_plot_subtype)

    summary_df <- plot_df %>%
      group_by(Year, fill_val = !!sym(input$clade_plot_fill)) %>%
      summarise(Count = sum(n), .groups = "drop")
    
    validate(need(nrow(summary_df) > 0, "No data available for the current filters."))
    
    fill_items <- sort(unique(summary_df$fill_val))
    
    actual_items <- setdiff(fill_items, "Unknown")
    if (input$clade_plot_palette == "rainbow") {
      my_colors <- setNames(grDevices::rainbow(length(actual_items)), actual_items)
    } else {
      my_colors <- setNames(viridis::viridis(length(actual_items), option = input$clade_plot_palette), actual_items)
    }
    if ("Unknown" %in% fill_items) {
      my_colors["Unknown"] <- "#d3d3d3"
    }
    
    plot_ly(summary_df, x = ~Year, y = ~Count, color = ~fill_val, colors = my_colors,
            type = "bar", hoverinfo = "text",
            text = ~paste0("Year: ", Year, "<br>", input$clade_plot_fill, ": ", fill_val, "<br>Count: ", scales::comma(Count)),
            marker = list(line = list(color = 'white', width = 0.5))) %>%
      layout(barmode = 'stack',
             xaxis = list(title = "Year", tickangle = -45, tickfont = list(family = "Arial", size = 12), tickformat = "d"),
             yaxis = list(title = "Sequence Count", tickformat = ","),
             legend = list(title = list(text = ""))) %>%
      config(displayModeBar = FALSE)
  })
  
  # ==========================================
  # SERVER: TAB 2 - SINGLE POSITION EXPLORER
  # ==========================================
  
  observeEvent(input$sp_quick_visit, {
    req(input$sp_quick_visit != "None")
    idx <- as.numeric(input$sp_quick_visit)
    row_data <- important_pos_df[idx, ]
    
    freezeReactiveValue(input, "global_subtype")
    freezeReactiveValue(input, "sp_gene")
    freezeReactiveValue(input, "sp_position")
    
    updateSelectInput(session, "global_subtype", selected = as.character(row_data$Subtype))
    updateSelectInput(session, "sp_gene", selected = as.character(row_data$Gene))
    updateNumericInput(session, "sp_position", value = as.numeric(row_data$Position))
  })
  
  # This reactiveVal will hold the data for the Single Position Explorer.
  # It acts as a buffer, allowing us to show a waiter during calculation.
  sp_data_val <- reactiveVal()
  sp_position_debounced <- debounce(reactive(input$sp_position), 800)

  # This observer triggers when any relevant input changes. It performs the heavy
  # calculation and shows a full-screen waiter while doing so.
  observeEvent(list(input$global_subtype, input$sp_gene, sp_position_debounced(), input$sp_group_by, input$variation_type, input$sp_min_seqs, input$sp_hide_empty_years, input$sp_year_month_start, input$sp_year_month_end), {
    subtype   <- input$global_subtype
    gene      <- input$sp_gene
    pos       <- sp_position_debounced()
    group_col <- input$sp_group_by 
    var_type  <- input$variation_type
    
    req(subtype, gene, pos, group_col, var_type)

    # Show the full-screen, modal waiter
    waiter_show(
      html = tagList(
        spin_fading_circles(),
        tags$h3("Loading Data...", style = "color:white; margin-top: 20px;"),
        tags$p("Please wait while we fetch and process the records.", style = "color:white;")
      ),
      color = "rgba(44, 62, 80, 0.9)"
    )
    on.exit(waiter_hide(), add = TRUE)

    if (usage_duckdb_available()) {
      ym_values <- sp_year_month_choices()
      allowed_yms <- get_selected_year_months(ym_values, input$sp_year_month_start, input$sp_year_month_end)

      filtered <- usage_single_position(
        subtype,
        var_type,
        gene,
        group_col,
        pos,
        allowed_yms = allowed_yms,
        min_seqs = input$sp_min_seqs,
        hide_empty_years = input$sp_hide_empty_years
      )

      if (is.null(filtered)) {
        sp_data_val(paste("Table not found for group:", group_col))
      } else {
        sp_data_val(filtered)
      }
      return()
    }
    
    rds_file <- count_cache_file_path(subtype, var_type, gene, group_col)
    data <- get_lazy_table(rds_file)
    
    # Handle cases where data is missing or inconsistent
    if (is.null(data)) { sp_data_val(paste("Table not found for group:", group_col)); return() }
    if (!(group_col %in% colnames(data))) { sp_data_val("Updating data..."); return() }

    data <- add_year_month_filter_column(data)
    
    # We always keep Group, Gene, Position, and the primary grouping column
    group_cols <- c("Group", "Gene", "Position", group_col)
    
    # 1. Basic filtering by gene, position, subtype
    filtered <- data %>% 
      filter(Group == subtype, Gene == gene, Position == pos) %>%
      {
        ym_values <- sp_year_month_choices()
        allowed_yms <- get_selected_year_months(ym_values, input$sp_year_month_start, input$sp_year_month_end)

        if (!is.null(allowed_yms)) {
          dplyr::filter(., Year_Month_Filter %in% allowed_yms)
        } else {
          .
        }
      } %>%
      # 2. Filter out "X" and "-"
      filter(!(AminoAcid %in% c("X", "-"))) %>%
      # NEW Step: Aggregate Counts by the grouping column and AminoAcid
      group_by(across(all_of(c(group_cols, "AminoAcid")))) %>%
      summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
      # 3. Recalculate totals and frequencies based on remaining valid sequences
      mutate(
        Valid_Total = sum(Count), 
        `Frequency(%)` = (Count / Valid_Total) * 100
      ) %>%
      ungroup()
    
    # 4. Apply minimum sequences filter based on Valid_Total
    filtered <- filtered %>% filter(Valid_Total >= input$sp_min_seqs)
    
    # 5. Clean up temporal metadata if necessary (KEEP 'Unknown' if requested)
    # if (group_col == "Year_Month") {
    #   filtered <- filtered %>% filter(Year_Month != "Unknown")
    # } else if (group_col == "Year") {
    #   filtered <- filtered %>% filter(Year != "Unknown")
    # }
    
    # 6. Optionally hide years without records (Valid_Total > 0)
    if (group_col == "Year" && input$sp_hide_empty_years) {
      filtered <- filtered %>% filter(Valid_Total > 0)
    }
    
    # Once calculation is done, update the reactiveVal
    sp_data_val(filtered)

  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  sp_filtered_data <- reactive({
    # This reactive now simply returns the value calculated in the observer.
    # Downstream consumers will automatically update when sp_data_val() changes.
    sp_data_val()
  })
  
  sp_plot_ggplot <- reactive({
    data <- sp_filtered_data()
    
    # If data is not a dataframe, it's a message from the observer (e.g., "Updating...").
    validate(need(is.data.frame(data), data))
    validate(need(nrow(data) > 0, "No data available for the current selection."))
    
    req(input$sp_font_size, input$sp_group_by)
    group_col <- input$sp_group_by
    show_counts <- isTRUE(input$sp_show_counts)
    y_col <- if(show_counts) "Count" else "Frequency(%)"
    y_lab <- if(show_counts) "Sequence Count" else "Frequency (%)"
    is_aa <- (input$variation_type == "AA")
    has_codon <- "Codon_Usage" %in% colnames(data)
    
    y_scale <- if(show_counts) {
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 
    } else {
      scale_y_continuous(expand = c(0, 0), limits = c(0, 105))  
    }
    
    # Enforce correct data type for X-axis to properly show or hide gaps
    # Define "null-ish" values to move to the front
    special_values <- c("Unknown", "unassigned", "Unassigned")

    if (group_col == "Year") {
      present_specials <- intersect(special_values, as.character(data[[group_col]]))
      has_specials <- length(present_specials) > 0
      
      if (input$sp_hide_empty_years || has_specials) {
        # Treat as categorical to hide gaps or handle Unknown/unassigned
        all_years <- sort(unique(as.character(data[[group_col]])))
        if (has_specials) {
          all_years <- c(present_specials, setdiff(all_years, present_specials))
        }
        data[[group_col]] <- factor(data[[group_col]], levels = all_years)
        x_scale <- scale_x_discrete()
      } else {
        # Treat as continuous to show gaps naturally
        data[[group_col]] <- as.numeric(as.character(data[[group_col]]))
        x_scale <- scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE))))))
      }
    } else if (group_col == "Year_Month") {
      all_yms <- sort(unique(as.character(data[[group_col]])))
      present_specials <- intersect(special_values, all_yms)
      has_specials <- length(present_specials) > 0
      
      if (has_specials) {
        all_yms <- c(present_specials, setdiff(all_yms, present_specials))
      }
      data[[group_col]] <- factor(data[[group_col]], levels = all_yms)
      
      # Select breaks (every 6th to avoid overlap), ensuring specials are shown if present
      if (has_specials) {
        # Always include specials, then sample every 6th from the chronological months
        chronological_yms <- setdiff(all_yms, special_values)
        sampled_yms <- every_nth_value(chronological_yms)
        x_scale <- scale_x_discrete(breaks = c(present_specials, sampled_yms))
      } else {
        x_scale <- scale_x_discrete(breaks = every_nth_value(all_yms))
      }
    } else {
      # Generic handling for other groups (Clades, etc.): Ensure Unknown/unassigned are first
      all_vals <- sort(unique(as.character(data[[group_col]])))
      present_specials <- intersect(special_values, all_vals)
      if (length(present_specials) > 0) {
        all_vals <- c(present_specials, setdiff(all_vals, present_specials))
      }
      data[[group_col]] <- factor(data[[group_col]], levels = all_vals)
      x_scale <- scale_x_discrete()
    }
    
    # Pre-calculate tooltip text
    data <- data %>%
      mutate(
        numbering_text = case_when(
          is_aa & input$sp_gene == "HA" & input$global_subtype == "H3N2" & Position <= 16 ~ " (Signal Peptide)",
          is_aa & input$sp_gene == "HA" & input$global_subtype == "H3N2" & Position > 16 & Position <= 345 ~ paste0(" (H3 HA1: ", Position - 16, ")"),
          is_aa & input$sp_gene == "HA" & input$global_subtype == "H3N2" & Position > 345 ~ paste0(" (H3 HA2: ", Position - 345, ")"),
          is_aa & input$sp_gene == "HA" & input$global_subtype == "H1N1" & Position <= 17 ~ " (Signal Peptide)",
          is_aa & input$sp_gene == "HA" & input$global_subtype == "H1N1" & Position > 17 & Position <= 344 ~ paste0(" (H1 HA1: ", Position - 17, ")"),
          is_aa & input$sp_gene == "HA" & input$global_subtype == "H1N1" & Position > 344 ~ paste0(" (H1 HA2: ", Position - 344, ")"),
          TRUE ~ ""
        ),
        tooltip_text = paste0(
          group_col, ": ", !!sym(group_col), 
          "<br>Position: ", Position, numbering_text,
          "<br>", if(is_aa) "Amino Acid: " else "Nucleotide: ", AminoAcid, 
          "<br>Count: ", Count, " / ", Valid_Total,
          "<br>Frequency: ", round(`Frequency(%)`, 2), "%"
        )
      )
    
    if (has_codon) {
      data <- data %>%
        mutate(tooltip_text = paste0(tooltip_text, "<br>Codons: ", Codon_Usage))
    }
    
    # FIX: Replaced .data[[]] with !!sym() and added group = AminoAcid
    ggplot(data, aes(x = !!sym(group_col), y = !!sym(y_col), fill = AminoAcid, group = AminoAcid,
                     text = tooltip_text)) + 
      geom_col(color = "black", size = 0.2) + # FIX: Changed linewidth to size
      scale_fill_manual(values = current_colors(), drop = FALSE) + 
      x_scale +
      y_scale + 
      theme_minimal(base_size = input$sp_font_size) + 
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        panel.grid.major.x = element_blank()
      ) +
      labs(x = group_col, y = y_lab, fill = variant_label())
  })
  
  output$sp_aa_plot <- renderPlotly({
    p <- sp_plot_ggplot()
    req(p)
    ggplotly(p, tooltip = "text") %>%
      config(displayModeBar = FALSE)
  })
  
  output$downloadSpPlot <- downloadHandler(
    filename = function() { 
      paste0(input$global_subtype, "_", input$sp_gene, "_Pos_", input$sp_position, "_Plot.", tolower(input$sp_plot_format)) 
    },
    content = function(file) { 
      ggsave(file, plot = sp_plot_ggplot(), 
             device = tolower(input$sp_plot_format), 
             width = 10, height = 5, dpi = 300) 
    }
  )
  
  output$sp_aa_table <- renderDT({
    data <- sp_filtered_data()
    
    # Handle messages from the observer (e.g., "Updating...") or empty dataframes
    validate(need(is.data.frame(data), data))
    validate(need(nrow(data) > 0, "No data to display for this selection."))
    
    req(input$sp_group_by)
    
    # Identify which columns are actually available in the data
    actual_cols <- colnames(data)
    
    # Select columns to show, prioritizing the current grouping column
    cols_to_show <- c(input$sp_group_by, "AminoAcid", "Count", "Valid_Total", "Frequency(%)")
    if("Codon_Usage" %in% actual_cols) cols_to_show <- c(cols_to_show, "Codon_Usage")
    
    # Final safety check: only use columns that exist in the dataframe
    final_cols <- intersect(cols_to_show, actual_cols)
    
    # Ensure we have at least one column to display
    validate(need(length(final_cols) > 0, "Processing data table..."))
    
    # Sort the table by the grouping column and then by frequency
    table_data <- data %>% 
      dplyr::select(all_of(final_cols))
    
    if (input$sp_group_by %in% colnames(table_data)) {
      table_data <- table_data %>% arrange(!!sym(input$sp_group_by), desc(`Frequency(%)`))
    }
    
    datatable(
      table_data, 
      options = list(pageLength = 10, autoWidth = TRUE, order = list()), 
      rownames = FALSE
    ) %>% formatRound("Frequency(%)", digits = 2)
  })
  
  output$sp_position_details <- renderUI({
    req(input$global_subtype, input$sp_gene, sp_position_debounced())
    # Only show important sites information in Amino Acid mode
    if (input$variation_type == "NT") return(NULL)
    
    match <- important_pos_df %>% 
      filter(Subtype == as.character(input$global_subtype), 
             Gene == as.character(input$sp_gene), 
             Position == sp_position_debounced())
    
    if(nrow(match) > 0) {
      wellPanel(style = "background-color: #e3f2fd; border-left: 5px solid #2196f3;",
                fluidRow(
                  column(2, strong("Mutation: "), match$Mutation),
                  column(2, strong("Epitope: "), match$Epitope),
                  column(2, strong("Impact: "), match$Clinical_impact),
                  column(4, strong("Source: "), em(match$Source))
                )
      )
    }
  })
  
  output$sp_range_label <- renderUI({
    req(input$global_subtype, input$sp_gene)

    if (usage_duckdb_available()) {
      gene_max <- usage_max_position(input$global_subtype, input$variation_type, input$sp_gene)
      req(!is.na(gene_max))
      return(tags$label(paste0(variant_label(), " Position (1 - ", gene_max, "):"),
                 `for` = "sp_position",
                 style = "display: block; margin-bottom: 5px; font-weight: bold; color: #2c3e50;"))
    }
    
    gene_max <- current_usage_by_clade() %>% 
      filter(Group == as.character(input$global_subtype), 
             Gene == as.character(input$sp_gene)) %>% 
      pull(Position) %>% 
      max(na.rm = TRUE)
    
    tags$label(paste0(variant_label(), " Position (1 - ", gene_max, "):"), 
               `for` = "sp_position", 
               style = "display: block; margin-bottom: 5px; font-weight: bold; color: #2c3e50;")
  })
  
  output$sp_numbering_label <- renderUI({
    req(input$global_subtype, input$sp_gene, input$sp_position, input$variation_type)
    
    # Only calculate structural numbering for Amino Acids in the HA gene
    if (input$variation_type == "AA" && input$sp_gene == "HA") {
      pos <- input$sp_position
      if (input$global_subtype == "H3N2") {
        if (pos <= 16) {
          return(span("(Signal Peptide)", style = "margin-left: 10px; color: #7f8c8d; font-style: italic;"))
        } else if (pos <= 345) {
          return(span(paste0("(H3 HA1: ", pos - 16, ")"), style = "margin-left: 10px; color: #e74c3c; font-weight: bold;"))
        } else {
          return(span(paste0("(H3 HA2: ", pos - 345, ")"), style = "margin-left: 10px; color: #e74c3c; font-weight: bold;"))
        }
      } else if (input$global_subtype == "H1N1") {
        if (pos <= 17) {
          return(span("(Signal Peptide)", style = "margin-left: 10px; color: #7f8c8d; font-style: italic;"))
        } else if (pos <= 344) {
          return(span(paste0("(H1 HA1: ", pos - 17, ")"), style = "margin-left: 10px; color: #e74c3c; font-weight: bold;"))
        } else {
          return(span(paste0("(H1 HA2: ", pos - 344, ")"), style = "margin-left: 10px; color: #e74c3c; font-weight: bold;"))
        }
      }
    }
    return(NULL)
  })
  
  observeEvent(list(input$sp_gene, input$variation_type), {
    req(input$global_subtype, input$sp_gene)

    if (usage_duckdb_available()) {
      max_val <- usage_max_position(input$global_subtype, input$variation_type, input$sp_gene)
      if (!is.na(max_val)) {
        updateNumericInput(session, "sp_position", max = max_val)
        updateSliderInput(session, "sp_pos_slider", max = max_val)
        if(input$sp_position > max_val) {
          updateNumericInput(session, "sp_position", value = max_val)
        }
      }
      return()
    }
    
    gene_data <- current_usage_by_clade() %>% 
      filter(Group == as.character(input$global_subtype), 
             Gene == as.character(input$sp_gene))
    
    if(nrow(gene_data) > 0) {
      max_val <- max(gene_data$Position, na.rm = TRUE)
      updateNumericInput(session, "sp_position", max = max_val)
      updateSliderInput(session, "sp_pos_slider", max = max_val)
      if(input$sp_position > max_val) {
        updateNumericInput(session, "sp_position", value = max_val)
      }
    }
  })
  
  observeEvent(input$sp_pos_plus, {
    req(input$global_subtype, input$sp_gene)

    if (usage_duckdb_available()) {
      gene_max <- usage_max_position(input$global_subtype, input$variation_type, input$sp_gene)
      req(!is.na(gene_max))
      updateNumericInput(session, "sp_position", value = min(gene_max, input$sp_position + 1))
      return()
    }

    gene_max <- current_usage_by_clade() %>% 
      filter(Group == as.character(input$global_subtype), Gene == as.character(input$sp_gene)) %>% 
      pull(Position) %>% max(na.rm = TRUE)
    
    new_val <- min(gene_max, input$sp_position + 1)
    updateNumericInput(session, "sp_position", value = new_val)
  })
  
  observeEvent(input$sp_pos_minus, {
    new_val <- max(1, input$sp_position - 1)
    updateNumericInput(session, "sp_position", value = new_val)
  })
  
  observeEvent(list(input$sp_position, input$variation_type), {
    req(input$global_subtype, input$sp_gene)

    if (usage_duckdb_available()) {
      gene_max <- usage_max_position(input$global_subtype, input$variation_type, input$sp_gene)
      req(!is.na(gene_max))
      if (!is.na(input$sp_position) && input$sp_position > gene_max) {
        updateNumericInput(session, "sp_position", value = gene_max)
        showNotification(
          paste("Maximum position for", input$sp_gene, "is", gene_max),
          type = "warning",
          duration = 2
        )
      }
      return()
    }
    
    gene_max <- current_usage_by_clade() %>% 
      filter(Group == as.character(input$global_subtype), 
             Gene == as.character(input$sp_gene)) %>% 
      pull(Position) %>% 
      max(na.rm = TRUE)
    
    if (!is.na(input$sp_position) && input$sp_position > gene_max) {
      updateNumericInput(session, "sp_position", value = gene_max)
      showNotification(
        paste("Maximum position for", input$sp_gene, "is", gene_max),
        type = "warning",
        duration = 2
      )
    }
  })
  
  # ==========================================
  # SERVER: TAB 2 - PAIRWISE COMPARISON 
  # ==========================================
  
  clicked_data_val <- reactiveValues(gene = NULL, pos = NULL, ready = FALSE, plot_id = NULL, table_id = NULL)
  
  # This reactiveVal will hold the data for the Pairwise Comparison table.
  pw_diff_val <- reactiveVal(NULL)

  observeEvent(list(input$global_subtype, input$variation_type, input$pw_group_by), {
    pw_diff_val(NULL)
    clicked_data_val$gene <- NULL
    clicked_data_val$pos <- NULL
    clicked_data_val$ready <- FALSE
    clicked_data_val$plot_id <- NULL
    clicked_data_val$table_id <- NULL
  }, ignoreInit = TRUE)

  # This observer runs the heavy comparison only after both groups are chosen.
  # Changing "Group by" is handled separately above to just refresh the dropdown options.
  observeEvent(list(
    input$pw_clade1, input$pw_clade2, input$pw_min_freq,
    session$clientData$output_pw_diff_table_hidden
  ), {
    # Prevent execution if the Pairwise Comparison tab is not currently visible
    if (!isFALSE(session$clientData$output_pw_diff_table_hidden)) return()
    
    req(input$global_subtype, input$variation_type, input$pw_group_by)
    
    if (is.null(input$pw_clade1) || is.null(input$pw_clade2) || input$pw_clade1 == "" || input$pw_clade2 == "") {
      pw_diff_val(NULL)
      return()
    }
    
    # --- FAST SANITY CHECK ---
    # Prevent premature execution: If the group changed, wait for the clade dropdowns to update first.
    genes <- available_count_genes(input$global_subtype, input$variation_type, prefer_cache = TRUE)
    if(length(genes) > 0) {
      gene <- if("HA" %in% genes) "HA" else genes[1]
      if (usage_duckdb_available()) {
        valid_clades <- usage_distinct_group_values(input$global_subtype, input$variation_type, gene, input$pw_group_by)
        if(length(valid_clades) > 0 && (!(input$pw_clade1 %in% c("", valid_clades)) || !(input$pw_clade2 %in% c("", valid_clades)))) return()
      } else {
      rds_file <- count_cache_file_path(input$global_subtype, input$variation_type, gene, input$pw_group_by)
      df <- get_lazy_table(rds_file)
      if(!is.null(df)) {
        valid_clades <- if(input$pw_group_by %in% colnames(df)) unique(as.character(df[[input$pw_group_by]])) else if("Clade" %in% colnames(df)) unique(as.character(df$Clade)) else "Unknown"
        if(!(input$pw_clade1 %in% c("", valid_clades)) || !(input$pw_clade2 %in% c("", valid_clades))) return()
      }
      }
    }
    
    if(input$pw_clade1 == input$pw_clade2) {
      pw_diff_val(data.frame(Message = "Reference and Target groups are identical. Please select different groups."))
      return()
    }

    # Show the full-screen, modal waiter
    waiter_show(
      html = tagList(
        spin_fading_circles(),
        tags$h3("Comparing Groups...", style = "color:white; margin-top: 20px;"),
        tags$p("Please wait while we identify robust differences.", style = "color:white;")
      ),
      color = "rgba(44, 62, 80, 0.9)"
    )
    on.exit(waiter_hide(), add = TRUE)
    
    pw_diff_val(compute_pairwise_differences(input$pw_clade1, input$pw_clade2, input$pw_min_freq))
    gc(FALSE)

  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # The original reactive now just returns the value from the observer
  pairwise_differences <- reactive({
    pw_diff_val()
  })
  
  output$pw_diff_table <- renderDT({
    data <- pairwise_differences()

    if (is.null(data)) {
      return(datatable(data.frame(Message = "Please select both Group 1 and Group 2 to begin comparison."), rownames = FALSE, options = list(dom = 't')))
    }
    
    # If the dataframe has a 'Message' column, it's an error/info message from the observer
    if ("Message" %in% colnames(data)) {
      return(datatable(data, rownames = FALSE, options = list(dom = 't')))
    }
    
    if(nrow(data) == 0) return(datatable(data.frame(Message = "No robust differences found for the current selection."), rownames = FALSE, options = list(dom = 't')))
    
    display_data <- data %>% 
      mutate(Position = sprintf('<a href="#" onclick="Shiny.setInputValue(\'modal_clicked\', \'%s|%s\', {priority: \'event\'});"><strong>%s</strong></a>', Gene, Position, Position)) %>% 
      dplyr::rename(`Group 1 AA` = Clade1_AA, `Group 1 Freq (%)` = Clade1_Freq, `Group 2 AA` = Clade2_AA, `Group 2 Freq (%)` = Clade2_Freq)
    
    datatable(display_data, escape = FALSE, options = list(pageLength = 15, autoWidth = TRUE), rownames = FALSE) %>% 
      formatRound(c("Group 1 Freq (%)", "Group 2 Freq (%)"), digits = 2)
  })
  
  observeEvent(input$modal_clicked, {
    parts <- strsplit(input$modal_clicked, "\\|")[[1]]
    next_gene <- parts[1]
    next_pos <- as.numeric(parts[2])
    modal_token <- paste0(format(Sys.time(), "%Y%m%d%H%M%OS6"), "_", sample.int(1000000, 1))
    next_plot_id <- paste0("modal_plot_", modal_token)
    next_table_id <- paste0("modal_table_", modal_token)
    
    clicked_data_val$gene <- NULL
    clicked_data_val$pos <- NULL
    clicked_data_val$ready <- FALSE
    clicked_data_val$plot_id <- NULL
    clicked_data_val$table_id <- NULL
    
    showModal(modalDialog(
      title = paste(variant_label(), "Usage: Gene", next_gene, "- Position", next_pos),
      size = "l", easyClose = TRUE,
      fluidRow(
        column(4, sliderInput("modal_font_size", "Plot Font Size:", min = 10, max = 24, value = 14, step = 1)),
        column(4, radioButtons("modal_plot_format", "Format:", choices = c("PDF", "PNG"), inline = TRUE)),
        column(4, downloadButton("downloadModalPlot", "Download Plot", class = "btn-info", style="margin-top: 25px; width: 100%;"))
      ),
      plotlyOutput(next_plot_id, height = "500px"), 
      hr(), 
      DTOutput(next_table_id)
    ))
    
    session$onFlushed(function() {
      clicked_data_val$gene <- next_gene
      clicked_data_val$pos <- next_pos
      clicked_data_val$plot_id <- next_plot_id
      clicked_data_val$table_id <- next_table_id
      clicked_data_val$ready <- TRUE
    }, once = TRUE)
  })
  
  modal_data <- reactive({ 
    req(isTRUE(clicked_data_val$ready), clicked_data_val$gene, clicked_data_val$pos)
    
    # Capture variation type once
    is_aa <- (input$variation_type == "AA")
    
    res <- get_pairwise_position_distribution(clicked_data_val$gene, clicked_data_val$pos)

    if (nrow(res) == 0) return(res)

    group_col <- input$pw_group_by
    has_codon <- "Codon_Usage" %in% colnames(res)
    
    # Pre-calculate a clean tooltip text to avoid complex logic inside aes()
    res <- res %>%
      mutate(
        numbering_text = case_when(
          is_aa & clicked_data_val$gene == "HA" & input$global_subtype == "H3N2" & clicked_data_val$pos <= 16 ~ " (Signal Peptide)",
          is_aa & clicked_data_val$gene == "HA" & input$global_subtype == "H3N2" & clicked_data_val$pos > 16 & clicked_data_val$pos <= 345 ~ paste0(" (H3 HA1: ", clicked_data_val$pos - 16, ")"),
          is_aa & clicked_data_val$gene == "HA" & input$global_subtype == "H3N2" & clicked_data_val$pos > 345 ~ paste0(" (H3 HA2: ", clicked_data_val$pos - 345, ")"),
          is_aa & clicked_data_val$gene == "HA" & input$global_subtype == "H1N1" & clicked_data_val$pos <= 17 ~ " (Signal Peptide)",
          is_aa & clicked_data_val$gene == "HA" & input$global_subtype == "H1N1" & clicked_data_val$pos > 17 & clicked_data_val$pos <= 344 ~ paste0(" (H1 HA1: ", clicked_data_val$pos - 17, ")"),
          is_aa & clicked_data_val$gene == "HA" & input$global_subtype == "H1N1" & clicked_data_val$pos > 344 ~ paste0(" (H1 HA2: ", clicked_data_val$pos - 344, ")"),
          TRUE ~ ""
        ),
        tooltip_text = as.character(paste0(
          group_col, ": ", Clade, 
          "<br>Position: ", clicked_data_val$pos, numbering_text,
          "<br>", if(is_aa) "Amino Acid: " else "Nucleotide: ", AminoAcid, 
          "<br>Frequency: ", round(!!sym("Frequency(%)"), 2), "%",
          "<br>Count: ", Count, " / ", Total_in_Clade
        ))
      )
      
    if (has_codon) {
      res <- res %>%
        mutate(tooltip_text = as.character(paste0(tooltip_text, "<br>Codons: ", Codon_Usage)))
    }
      
    return(res)
  })
  
  library(ggtext) 
  
  modal_plot_ggplot <- reactive({
    req(isTRUE(clicked_data_val$ready), input$modal_font_size, input$pw_clade1, input$pw_clade2)
    data <- modal_data()
    validate(need(nrow(data) > 0, "No data available."))
    
    selected_clades <- c(input$pw_clade1, input$pw_clade2)
    plot_clades <- sort(unique(data$Clade))
    
    html_labels <- ifelse(
      plot_clades %in% selected_clades, 
      paste0("<b style='color:red;'>", plot_clades, "</b>"), 
      plot_clades
    )
    names(html_labels) <- plot_clades 
    
    group_col <- input$pw_group_by
    special_values <- c("Unknown", "unassigned", "Unassigned")
    
    # Enforce correct data type to hide/show gaps dynamically
    if (group_col == "Year") {
      present_specials <- intersect(special_values, as.character(data$Clade))
      has_specials <- length(present_specials) > 0
      
      if (isTRUE(input$pw_hide_empty_years) || has_specials) {
        all_years <- sort(unique(as.character(data$Clade)))
        if (has_specials) all_years <- c(present_specials, setdiff(all_years, present_specials))
        data$Clade <- factor(data$Clade, levels = all_years)
        x_scale <- scale_x_discrete(labels = html_labels)
      } else {
        data$Clade <- as.numeric(as.character(data$Clade))
        safe_sel_clades <- suppressWarnings(as.numeric(selected_clades))
        x_scale <- scale_x_continuous(
          breaks = function(x) {
            bks <- unique(floor(pretty(seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE)))))
            sort(unique(c(bks, safe_sel_clades[!is.na(safe_sel_clades)])))
          },
          labels = function(b) {
            ifelse(b %in% safe_sel_clades, paste0("<b style='color:red;'>", b, "</b>"), as.character(b))
          }
        )
      }
    } else if (group_col == "Year_Month") {
      all_yms <- sort(unique(as.character(data$Clade)))
      present_specials <- intersect(special_values, all_yms)
      has_specials <- length(present_specials) > 0
      if (has_specials) all_yms <- c(present_specials, setdiff(all_yms, present_specials))
      data$Clade <- factor(data$Clade, levels = all_yms)
      if (has_specials) {
        chronological_yms <- setdiff(all_yms, special_values)
        sampled_yms <- every_nth_value(chronological_yms)
        bks <- unique(c(present_specials, sampled_yms, selected_clades[selected_clades %in% all_yms]))
        # Retain original chronological ordering for the breaks
        bks <- all_yms[all_yms %in% bks] 
        x_scale <- scale_x_discrete(breaks = bks, labels = unname(html_labels[bks]))
      } else {
        bks <- every_nth_value(all_yms)
        bks <- unique(c(bks, selected_clades[selected_clades %in% all_yms]))
        bks <- all_yms[all_yms %in% bks] 
        x_scale <- scale_x_discrete(breaks = bks, labels = unname(html_labels[bks]))
      }
    } else {
      all_vals <- sort(unique(as.character(data$Clade)))
      present_specials <- intersect(special_values, all_vals)
      if (length(present_specials) > 0) all_vals <- c(present_specials, setdiff(all_vals, present_specials))
      data$Clade <- factor(data$Clade, levels = all_vals)
      x_scale <- scale_x_discrete(labels = html_labels)
    }
    
    ggplot(data, aes(x = Clade, y = !!sym("Frequency(%)"), fill = AminoAcid, group = AminoAcid, 
                     text = tooltip_text)) + 
      geom_col(color = "black", size = 0.2) + 
      scale_fill_manual(values = current_colors(), drop = FALSE) + 
      scale_y_continuous(expand = c(0,0), limits = c(0, 105)) + 
      x_scale + 
      labs(x = input$pw_group_by, y = "Frequency (%)", fill = variant_label()) + 
      theme_minimal(base_size = input$modal_font_size) + 
      theme(
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1), 
        axis.title = element_text(face = "bold"), 
        panel.grid.major.x = element_blank()
      )
  })
  
  observe({
    req(clicked_data_val$plot_id)
    local({
      plot_id <- clicked_data_val$plot_id
      output[[plot_id]] <- renderPlotly({
        if (!isTRUE(clicked_data_val$ready) || !identical(clicked_data_val$plot_id, plot_id)) {
          return(plotly_empty() %>% layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
        }
        
        ggplotly(modal_plot_ggplot(), tooltip = "text") %>%
          config(displayModeBar = FALSE)
      })
    })
  })
  
  output$downloadModalPlot <- downloadHandler(
    filename = function() { 
      paste0(input$global_subtype, "_", clicked_data_val$gene, "_Pos_", clicked_data_val$pos, "_Plot.", tolower(input$modal_plot_format)) 
    },
    content = function(file) { 
      ggsave(file, plot = modal_plot_ggplot(), 
             device = tolower(input$modal_plot_format), 
             width = 10, height = 5, dpi = 300) 
    }
  )
  
  observe({
    req(clicked_data_val$table_id)
    local({
      table_id <- clicked_data_val$table_id
      output[[table_id]] <- renderDT({
        if (!isTRUE(clicked_data_val$ready) || !identical(clicked_data_val$table_id, table_id)) {
          return(datatable(data.frame(Message = "Loading selected position..."), rownames = FALSE, options = list(dom = 't')))
        }
        
        data <- modal_data()
        cols_to_show <- c("Clade", "AminoAcid", "Count", "Total_in_Clade", "Frequency(%)")
        if("Codon_Usage" %in% colnames(data)) cols_to_show <- c(cols_to_show, "Codon_Usage")
        
        display_data <- data %>%
          dplyr::select(all_of(cols_to_show)) %>%
          arrange(Clade, desc(`Frequency(%)`)) %>%
          dplyr::rename(`Group` = Clade)
        
        datatable(display_data, options = list(pageLength = 5, autoWidth = TRUE), rownames = FALSE) %>% formatRound("Frequency(%)", digits = 2) 
      })
    })
  })
  
  output$downloadPairwiseCSV <- downloadHandler(
    filename = function() { paste0("Differences_", input$pw_clade1, "_vs_", input$pw_clade2, ".csv") },
    content = function(file) { 
      data <- pairwise_differences()
      if(nrow(data)>0 && !("Message" %in% colnames(data))) data <- data %>% mutate(Group = input$global_subtype) %>% dplyr::select(Group, everything())
      write_csv(data, file) 
    }
  )
  
  output$downloadPairwiseExcel <- downloadHandler(
    filename = function() { paste0("Matrices_", input$pw_clade1, "_vs_", input$pw_clade2, ".xlsx") },
    content = function(file) {
      diffs <- pairwise_differences(); wb <- createWorkbook()
      if(nrow(diffs) == 0 || "Message" %in% colnames(diffs)) { 
        addWorksheet(wb, "No Differences"); writeData(wb, "No Differences", "No differences found."); saveWorkbook(wb, file, overwrite = TRUE); return() 
      }
      base_df <- data.frame(AminoAcid = if(input$variation_type == "AA") ALL_AAS else c("a","c","g","t","A","C","G","T","N","n","-"))
      current_gene <- NULL
      current_gene_data <- NULL

      for(i in 1:nrow(diffs)) {
        r_gene <- diffs$Gene[i]; r_pos <- diffs$Position[i]

        if (usage_duckdb_available()) {
          pos_data <- get_pairwise_position_distribution(r_gene, r_pos)
        } else {
          if (!identical(current_gene, r_gene)) {
          current_gene <- r_gene
          current_gene_data <- load_pairwise_gene_data(r_gene)
          }

          pos_data <- current_gene_data %>%
            filter(Position == r_pos) %>%
            filter(!(AminoAcid %in% c("X", "-"))) %>%
            group_by(Clade, AminoAcid) %>%
            summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
            mutate(`Frequency(%)` = (Count / sum(Count)) * 100) %>%
            ungroup()
        }
        
        sorted_clades <- sort(unique(pos_data$Clade))
        pct_matrix <- left_join(base_df, pos_data %>% dplyr::select(Clade, AminoAcid, `Frequency(%)`) %>% pivot_wider(names_from = Clade, values_from = `Frequency(%)`, values_fill = 0), by = "AminoAcid")
        pct_matrix[is.na(pct_matrix)] <- 0; pct_matrix <- pct_matrix[, c("AminoAcid", sorted_clades)]
        cnt_matrix <- left_join(base_df, pos_data %>% dplyr::select(Clade, AminoAcid, Count) %>% pivot_wider(names_from = Clade, values_from = Count, values_fill = 0), by = "AminoAcid")
        cnt_matrix[is.na(cnt_matrix)] <- 0; cnt_matrix <- cnt_matrix[, c("AminoAcid", sorted_clades)]
        
        sheet_name <- substr(paste(r_gene, "Pos", r_pos), 1, 31); addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, "Percentage (%)", startRow=1, startCol=1); writeData(wb, sheet_name, pct_matrix, startRow=2, startCol=1)
        start_count_row <- 2 + nrow(pct_matrix) + 2; writeData(wb, sheet_name, "Count", startRow=start_count_row, startCol=1); writeData(wb, sheet_name, cnt_matrix, startRow=start_count_row+1, startCol=1)
        
        num_cols <- ncol(pct_matrix)
        addStyle(wb, sheet_name, style = createStyle(numFmt = "0.00"), rows = 3:(2 + nrow(pct_matrix)), cols = 2:num_cols, gridExpand = TRUE)
        addStyle(wb, sheet_name, style = createStyle(numFmt = "0"), rows = (start_count_row + 2):(start_count_row + 1 + nrow(cnt_matrix)), cols = 2:num_cols, gridExpand = TRUE)
        conditionalFormatting(wb, sheet_name, cols = 2:num_cols, rows = 3:(2 + nrow(pct_matrix)), style = c("#FFFFFF", "#238B45"), type = "colourScale")
        headerStyle <- createStyle(textDecoration = "bold"); addStyle(wb, sheet_name, style = headerStyle, rows = c(1, start_count_row), cols = 1); addStyle(wb, sheet_name, style = headerStyle, rows = c(2, start_count_row + 1), cols = 1:num_cols, gridExpand = TRUE)
        highlight_cols <- which(colnames(pct_matrix) %in% c(input$pw_clade1, input$pw_clade2))
        if(length(highlight_cols) > 0) addStyle(wb, sheet_name, style = createStyle(fontColour = "#FF0000", textDecoration = "bold"), rows = c(2, start_count_row + 1), cols = highlight_cols, gridExpand = TRUE)
      }

      rm(current_gene, current_gene_data)
      gc(FALSE)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ==========================================
  # SERVER: TAB 3 - ENTROPY LANDSCAPE
  # ==========================================
  output$ent_plot_title <- renderText({ 
    clade_text <- if(input$ent_group == "All") paste("All", input$ent_group_by) else paste(input$ent_group_by, input$ent_group)
    mode_text <- if(input$variation_type == "AA") "Amino Acid" else "Nucleotide"
    paste(mode_text, "Shannon Entropy Landscape - Subtype", input$global_subtype, "| Gene", input$ent_gene, "|", clade_text) 
  })
  
  output$ent_plot <- renderPlotly({
    req(input$global_subtype, input$ent_gene, input$ent_group)

    if (usage_duckdb_available()) {
      ent_data <- usage_entropy_data(input$global_subtype, input$variation_type, input$ent_gene, input$ent_group_by, input$ent_group)
    } else {
      tmp <- ent_usage_data() %>%
        filter(Group == input$global_subtype, Gene == input$ent_gene)

      if (input$ent_group != "All") {
        tmp <- tmp %>% filter(Clade == input$ent_group)
      }

      ent_data <- tmp %>%
        # NEW: Remove "X" and "-" to ensure entropy only measures valid biological variation
        filter(!(AminoAcid %in% c("X", "-"))) %>%
        group_by(Position, AminoAcid) %>%
        summarise(AA_Sum = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
        mutate(
          Pos_Total = sum(AA_Sum),
          p = AA_Sum / Pos_Total
        ) %>%
        filter(p > 0) %>%
        summarise(Entropy = -sum(p * log2(p)), .groups = "drop")
    }
    
    validate(need(nrow(ent_data) > 0, "No data available for these selections after filtering unknowns."))
    
    # NT max entropy is 2 bits (log2(4)), AA max is ~4.39 bits (log2(21))
    y_max <- if(input$variation_type == "AA") log2(21) else 2.0
    
    # Define thresholds based on mode
    mid_thresh <- if(input$variation_type == "AA") 0.2 else 0.1
    high_thresh <- if(input$variation_type == "AA") 1.0 else 0.5
    
    plot_ly(data = ent_data, x = ~Position, y = ~Entropy) %>%
      add_bars(  
        marker = list(color = '#2980b9'), 
        name = "Entropy",
        hoverinfo = "text",
        text = ~paste0("Position: ", Position, "<br>Entropy: ", round(Entropy, 4), " bits")
      ) %>%
      layout(
        xaxis = list(title = paste(variant_label(), "Position"), automargin = TRUE),
        yaxis = list(title = "Shannon Entropy (Bits)", range = c(0, y_max)), 
        hovermode = "closest",
        font = list(size = input$ent_font_size),
        margin = list(l = 60, r = 40, b = 60, t = 40),
        autosize = TRUE,
        
        # ADD HORIZONTAL LINES
        shapes = list(
          list(
            type = "line",
            x0 = 0, x1 = 1, xref = "paper", # Spans the whole width
            y0 = mid_thresh, y1 = mid_thresh, yref = "y",
            line = list(color = "orange", dash = "dash", width = 1.5)
          ),
          list(
            type = "line",
            x0 = 0, x1 = 1, xref = "paper", 
            y0 = high_thresh, y1 = high_thresh, yref = "y",
            line = list(color = "red", dash = "dash", width = 1.5)
          )
        ),
        
        # ADD LABELS FOR THE LINES
        annotations = list(
          list(
            x = 1, y = mid_thresh + 0.05, xref = "paper", yref = "y",
            text = "Mid Variant", showarrow = FALSE,
            xanchor = "right", yanchor = "bottom", 
            font = list(color = "orange", size = 10)
          ),
          list(
            x = 1, y = high_thresh + 0.05, xref = "paper", yref = "y",
            text = "High Variant", showarrow = FALSE,
            xanchor = "right", yanchor = "bottom", 
            font = list(color = "red", size = 10)
          )
        )
      ) %>%
      config(displayModeBar = FALSE)
  })
  
  # ==========================================
  # SERVER: TAB 4 - MUTATION LOLLIPOP
  # ==========================================
  output$lol_plot_title <- renderText({ 
    mode_text <- if(input$variation_type == "AA") "Amino Acid" else "Nucleotide"
    paste(mode_text, "Mutations in", input$lol_tar_group, "vs Reference", input$lol_ref_group, "(Gene", input$lol_gene, ")") 
  })
  
  lol_plot_object <- reactive({
    req(input$global_subtype, input$lol_gene, input$lol_ref_group, input$lol_tar_group)
    validate(need(input$lol_ref_group != input$lol_tar_group, paste("Reference and Target groups are identical. Please select different", input$lol_group_by, "groups.")))

    if (usage_duckdb_available()) {
      consensus <- usage_lollipop_consensus(
        input$global_subtype,
        input$variation_type,
        input$lol_gene,
        input$lol_group_by,
        input$lol_ref_group,
        input$lol_tar_group,
        input$lol_min_freq
      )
      validate(need(!is.null(consensus), "No fixed mutations found between these groups."))

      c1 <- consensus %>%
        filter(Clade == input$lol_ref_group) %>%
        dplyr::select(Position, Ref_AA = AminoAcid)
      c2 <- consensus %>%
        filter(Clade == input$lol_tar_group) %>%
        dplyr::select(Position, Tar_AA = AminoAcid)
    } else {
      # 1. Fetch Reference Clade Data (c1)
      c1 <- lol_usage_data() %>%
        filter(Group == input$global_subtype, Clade == input$lol_ref_group, Gene == input$lol_gene) %>%
        # Step A: Exclude "X" and "-"
        filter(!(AminoAcid %in% c("X", "-"))) %>%
        # NEW Step: Aggregate Counts across sub-groups (Year, Month, etc.) to get clade-wide totals per position
        group_by(Position, AminoAcid) %>%
        summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
        # Step B: Recalculate frequencies based on valid amino acids
        group_by(Position) %>%
        mutate(
          Valid_Total = sum(Count),
          New_Frequency = (Count / Valid_Total) * 100
        ) %>%
        # Step C: Find the consensus residue
        filter(New_Frequency == max(New_Frequency)) %>%
        filter(row_number() == 1, New_Frequency >= input$lol_min_freq) %>%
        ungroup() %>%
        dplyr::select(Position, Ref_AA = AminoAcid)

      # 2. Fetch Target Clade Data (c2)
      c2 <- lol_usage_data() %>%
        filter(Group == input$global_subtype, Clade == input$lol_tar_group, Gene == input$lol_gene) %>%
        # Step A: Exclude "X" and "-"
        filter(!(AminoAcid %in% c("X", "-"))) %>%
        # NEW Step: Aggregate Counts across sub-groups
        group_by(Position, AminoAcid) %>%
        summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
        # Step B: Recalculate frequencies based on valid amino acids
        group_by(Position) %>%
        mutate(
          Valid_Total = sum(Count),
          New_Frequency = (Count / Valid_Total) * 100
        ) %>%
        # Step C: Find the consensus residue
        filter(New_Frequency == max(New_Frequency)) %>%
        filter(row_number() == 1, New_Frequency >= input$lol_min_freq) %>%
        ungroup() %>%
        dplyr::select(Position, Tar_AA = AminoAcid)
    }
    
    muts <- inner_join(c1, c2, by = "Position") %>% 
      filter(Ref_AA != Tar_AA) %>% 
      arrange(Position) %>%
      mutate(
        Label = paste0(Ref_AA, Position, Tar_AA),
        Y_Level = rep(c(1.0, 1.4, 1.8, 2.2), length.out = n()),
        HoverText = paste("Position:", Position, "<br>Mutation:", Label)
      )
    
    validate(need(nrow(muts) > 0, "No fixed mutations found between these groups."))
    
    # FIX: Added 'text' to all geoms so Plotly tooltips don't crash when scanning layers
    ggplot(muts, aes(x = Position, y = Y_Level)) +
      geom_segment(aes(xend = Position, yend = 0, text = HoverText), color = "gray60", size = 1) +
      geom_point(aes(fill = Tar_AA, text = HoverText), size = 5, shape = 21, color = "black") +
      geom_text(aes(y = Y_Level + 0.2, label = Label, text = HoverText), size = input$lol_font_size / 3) +
      scale_fill_manual(values = current_colors(), drop = FALSE) +
      scale_y_continuous(limits = c(0, 3.0), breaks = NULL) + 
      labs(x = paste(variant_label(), "Position"), y = "", fill = paste("New", variant_label())) +
      theme_minimal(base_size = input$lol_font_size) +
      theme(axis.title.x = element_text(face = "bold"), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
  })
  
  output$lol_plot <- renderPlotly({
    suppressWarnings(ggplotly(lol_plot_object(), tooltip = "text"))
  })
  
  output$downloadLolPlot <- downloadHandler(
    filename = function() { paste0("Lollipop_", input$lol_ref_group, "_vs_", input$lol_tar_group, "_", input$lol_gene, ".", tolower(input$lol_plot_format)) },
    content = function(file) { ggsave(file, plot = lol_plot_object(), device = tolower(input$lol_plot_format), width = 12, height = 6, dpi = 300) }
  )
  
  # ==========================================
  # SERVER: TAB 5 - CONSENSUS HEATMAP (msaR)
  # ==========================================
  # output$heat_plot_title <- renderText({ 
  #   mode_text <- if(input$variation_type == "AA") "Amino Acid" else "Nucleotide"
  #   paste("Interactive Consensus", mode_text, "MSA - Subtype", input$global_subtype, "| Gene", input$heat_gene) 
  # })
  # 
  # output$msa_dynamic_container <- renderUI({
  #   req(input$global_subtype)
  #   clade_count <- heat_usage_data() %>% filter(Group == input$global_subtype) %>% pull(Clade) %>% unique() %>% length()
  #   
  #   if (!is.null(input$show_mut_only) && input$show_mut_only) {
  #     clade_count <- clade_count + 1
  #   }
  #   
  #   outer_height <- (clade_count * 20) + 150
  #   msaROutput("heat_plot", width = "100%", height = paste0(outer_height, "px"))
  # })
  # 
  # output$heat_plot <- renderMsaR({
  #   req(input$global_subtype, input$heat_gene)
  #   
  #   cache_dir <- "data"
  #   if (!dir.exists(cache_dir)) dir.create(cache_dir, showWarnings = FALSE)
  #   
  #   safe_subtype <- gsub("[^A-Za-z0-9_]", "_", input$global_subtype)
  #   safe_gene <- gsub("[^A-Za-z0-9_]", "_", input$heat_gene)
  #   safe_group <- gsub("[^A-Za-z0-9_]", "_", input$heat_group_by)
  #   
  #   prefix <- if(input$variation_type == "AA") "MSA_AA_" else "MSA_NT_"
  #   # Filename no longer includes minFreq as it is removed
  #   aln_filename <- file.path(cache_dir, paste0(prefix, safe_subtype, "_", safe_gene, "_", safe_group, ".fasta"))
  #   
  #   if (file.exists(aln_filename)) {
  #     if(input$variation_type == "AA") {
  #       aligned_strings <- Biostrings::readAAStringSet(aln_filename)
  #     } else {
  #       aligned_strings <- Biostrings::readDNAStringSet(aln_filename)
  #     }
  #     original_clade_order <- names(aligned_strings) 
  #     
  #   } else {
  #     # 1. Prepare Base Grid of All Positions and Clades to ensure NO GAPS
  #     all_pos_in_gene <- heat_usage_data() %>% 
  #       filter(Group == input$global_subtype, Gene == input$heat_gene) %>% 
  #       pull(Position) %>% unique() %>% sort()
  #     
  #     all_clades_in_group <- heat_usage_data() %>% 
  #       filter(Group == input$global_subtype, Gene == input$heat_gene) %>% 
  #       pull(Clade) %>% unique() %>% sort()
  #     
  #     grid <- expand.grid(Clade = all_clades_in_group, Position = all_pos_in_gene, stringsAsFactors = FALSE)
  #     
  #     # 2. Filter and Identify Majority Character per Position/Clade
  #     exclude_chars <- if(input$variation_type == "AA") c("X", "-") else c("N", "n", "-")
  #     
  #     raw_data <- heat_usage_data() %>% 
  #       filter(Group == input$global_subtype, Gene == input$heat_gene) %>%
  #       # Filter out ambiguous/gaps for calculation
  #       filter(!(AminoAcid %in% exclude_chars)) %>% 
  #       group_by(Clade, Position, AminoAcid) %>%
  #       summarise(Total_Count = sum(Count), .groups = "drop") %>%
  #       group_by(Clade, Position) %>%
  #       filter(Total_Count == max(Total_Count)) %>%
  #       filter(row_number() == 1) %>%
  #       ungroup()
  #     
  #     # 3. Merge with Grid and Fill Gaps with "-"
  #     complete_data <- grid %>%
  #       left_join(raw_data, by = c("Clade", "Position")) %>%
  #       mutate(AminoAcid = ifelse(is.na(AminoAcid), "-", AminoAcid)) %>%
  #       arrange(Clade, Position)
  #     
  #     # 4. Reconstruct the consensus sequences
  #     raw_seqs <- complete_data %>%
  #       group_by(Clade) %>%
  #       summarise(seq = toupper(paste(AminoAcid, collapse = "")), .groups = "drop") %>%
  #       arrange(Clade)
  #     
  #     original_clade_order <- raw_seqs$Clade
  #     
  #     if(input$variation_type == "AA") {
  #       unaligned_strings <- Biostrings::AAStringSet(setNames(raw_seqs$seq, original_clade_order))
  #     } else {
  #       unaligned_strings <- Biostrings::DNAStringSet(setNames(raw_seqs$seq, original_clade_order))
  #     }
  #     
  #     waiter_show(
  #       html = tagList(
  #         spin_fading_circles(), 
  #         h3(paste("Running ClustalW Alignment for", input$heat_gene, "sequences..."), style = "color: white; margin-top: 20px;")
  #       ),
  #       color = "rgba(44, 62, 80, 0.9)"
  #     )
  #     
  #     on.exit(waiter_hide(), add = TRUE) 
  #     
  #     aligned_msa <- suppressMessages(msa::msa(unaligned_strings, method = "ClustalW", order = "input"))
  #     
  #     if(input$variation_type == "AA") {
  #       aligned_strings <- as(aligned_msa, "AAStringSet")
  #     } else {
  #       aligned_strings <- as(aligned_msa, "DNAStringSet")
  #     }
  #     aligned_strings <- aligned_strings[original_clade_order]
  #     
  #     Biostrings::writeXStringSet(aligned_strings, filepath = aln_filename)
  #     
  #   } 
  #   
  #   if (input$show_mut_only) {
  #     seq_char_matrix <- as.matrix(aligned_strings)
  #     exclude_chars <- if(input$variation_type == "AA") c("X", "-") else c("N", "n", "-")
  #     
  #     consensus_seq <- apply(seq_char_matrix, 2, function(col) {
  #       # Filter out X/- or N/- before finding the majority
  #       valid_col <- col[!(toupper(col) %in% toupper(exclude_chars))]
  #       if(length(valid_col) == 0) return("-") # Fallback if all are ambiguous
  #       freqs <- table(valid_col)
  #       names(freqs)[which.max(freqs)] 
  #     })
  #     
  #     for (i in 1:nrow(seq_char_matrix)) {
  #       match_idx <- seq_char_matrix[i, ] == consensus_seq
  #       seq_char_matrix[i, match_idx] <- "."
  #     }
  #     
  #     consensus_string <- paste(consensus_seq, collapse = "")
  #     new_seqs <- apply(seq_char_matrix, 1, paste, collapse = "")
  #     
  #     if(input$variation_type == "AA") {
  #       final_strings <- Biostrings::AAStringSet(c(Consensus = consensus_string, new_seqs[original_clade_order]))
  #     } else {
  #       final_strings <- Biostrings::DNAStringSet(c(Consensus = consensus_string, new_seqs[original_clade_order]))
  #     }
  #     
  #   } else {
  #     final_strings <- aligned_strings
  #   }
  #   
  #   inner_align_height <- (length(final_strings) * 20) + 20
  #   
  #   msaR(
  #     final_strings, 
  #     menu = TRUE, 
  #     overviewbox = FALSE, 
  #     seqlogo = !isTRUE(input$show_mut_only), 
  #     colorscheme = if(input$variation_type == "AA") "clustal" else "nucleotide",
  #     alignmentHeight = inner_align_height
  #     )
  #     })

      # --- HIDE LOADING CURTAIN WHEN READY ---
      session$onFlushed(function() {
        waiter_hide()
      }, once = TRUE)

      # --- PRE-RENDER HIDDEN TABS FOR INSTANT UX ---
      outputOptions(output, "stats_time_plot", suspendWhenHidden = FALSE)
      outputOptions(output, "stats_geo_plot", suspendWhenHidden = FALSE)
      outputOptions(output, "stats_clade_plot", suspendWhenHidden = FALSE)
      }
