# ==========================================
# 3. SERVER LOGIC
# ==========================================
server <- function(input, output, session) {
  
  # --- REACTIVE DATA SWITCH ---
  
  current_colors <- reactive({
    if(input$variation_type == "AA") aa_colors else nt_colors
  })
  variant_label <- reactive({
    if(input$variation_type == "AA") "AA" else "NT"
  })
  
  current_usage_by_clade <- reactive({
    req(input$global_subtype, input$variation_type)
    var_lower <- tolower(input$variation_type)
    
    dir_path <- paste0("data/", input$global_subtype, "/", input$variation_type, "/")
    rds_file <- paste0(dir_path, var_lower, "_usage_by_HA_clade.rds")
    
    # Fallback to the first available group file if HA_clade is not available for this subtype
    if (!file.exists(rds_file)) {
      files <- list.files(dir_path, pattern = paste0("^", var_lower, "_usage_by_.*\\.rds$"))
      files <- setdiff(files, c(paste0(var_lower, "_usage_by_Year.rds"), paste0(var_lower, "_usage_by_Year_Month.rds")))
      if (length(files) > 0) {
        rds_file <- paste0(dir_path, files[1])
      }
    }
    
    df <- get_lazy_table(rds_file)
    if (!is.null(df)) {
      if ("HA_clade" %in% colnames(df)) {
        df <- df %>% dplyr::rename(Clade = HA_clade)
      } else {
        group_col <- sub(paste0("^", var_lower, "_usage_by_(.*)\\.rds$"), "\\1", basename(rds_file))
        if (group_col %in% colnames(df)) {
          df <- df %>% dplyr::rename(Clade = !!sym(group_col))
        }
      }
      return(df)
    }
    return(data.frame(Group=character(), Gene=character(), Clade=character(), Position=numeric(), AminoAcid=character(), Count=numeric()))
  })

  # --- HELPER: Update Grouping Choices based on Loaded Data ---
  observe({
    req(input$global_subtype, input$variation_type)
    dir_path <- paste0("data/", input$global_subtype, "/", input$variation_type, "/")
    files <- list.files(dir_path, pattern = paste0("^", tolower(input$variation_type), "_usage_by_.*\\.rds$"))
    available_groups <- sub(paste0("^", tolower(input$variation_type), "_usage_by_(.*)\\.rds$"), "\\1", files)
    
    if (length(available_groups) == 0) return()

    # Create a mapping for display names
    # Key = Internal Key (Year), Value = Display Label (Year)
    group_map <- setNames(available_groups, available_groups)
    
    # Custom labels for better readability
    if ("Year_Month" %in% names(group_map)) names(group_map)[group_map == "Year_Month"] <- "Year-Month"
    # Convert underscores to spaces and capitalize for others (e.g., HA_clade -> HA Clade)
    other_indices <- which(!(group_map %in% c("Year", "Year_Month")))
    for (i in other_indices) {
      names(group_map)[i] <- gsub("_", " ", group_map[i])
    }
    
    # Desired priority order: 1. Year, 2. Year_Month, then others
    priority_keys <- intersect(c("Year", "Year_Month"), available_groups)
    other_keys <- sort(setdiff(available_groups, priority_keys))
    
    final_ordered_keys <- c(priority_keys, other_keys)
    
    # Create the named list for choices: list("Display" = "internal_key")
    final_choices <- setNames(final_ordered_keys, names(group_map)[match(final_ordered_keys, group_map)])
    
    # Determine selection: keep current if still valid, else default to Year
    current_sel <- if (input$sp_group_by %in% available_groups) input$sp_group_by else "Year"
    
    updateSelectInput(session, "sp_group_by", choices = final_choices, selected = current_sel)
  })
  
  # --- HELPER: Update Gene Dropdowns based on Subtype ---
  observeEvent(list(input$global_subtype, input$variation_type), { updateSelectInput(session, "sp_gene", choices = current_usage_by_clade() %>% filter(Group == input$global_subtype) %>% pull(Gene) %>% unique() %>% sort()) })
  observeEvent(list(input$global_subtype, input$variation_type), { updateSelectInput(session, "ent_gene", choices = current_usage_by_clade() %>% filter(Group == input$global_subtype) %>% pull(Gene) %>% unique() %>% sort()) })
  observeEvent(list(input$global_subtype, input$variation_type), { updateSelectInput(session, "lol_gene", choices = current_usage_by_clade() %>% filter(Group == input$global_subtype) %>% pull(Gene) %>% unique() %>% sort()) })
  observeEvent(list(input$global_subtype, input$variation_type), { updateSelectInput(session, "heat_gene", choices = current_usage_by_clade() %>% filter(Group == input$global_subtype) %>% pull(Gene) %>% unique() %>% sort()) })
  
  # --- HELPER: Update Clade Dropdowns based on Subtype ---
  observeEvent(list(input$global_subtype, input$variation_type), {
    clades <- current_usage_by_clade() %>% filter(Group == input$global_subtype) %>% pull(Clade) %>% unique() %>% sort()
    updateSelectInput(session, "pw_clade1", choices = clades, selected = clades[1])
    updateSelectInput(session, "pw_clade2", choices = clades, selected = if(length(clades)>1) clades[2] else clades[1])
  })
  observeEvent(list(input$global_subtype, input$ent_gene, input$variation_type), {
    req(input$global_subtype, input$ent_gene)
    clade_choices <- current_usage_by_clade() %>% 
      filter(Group == input$global_subtype, Gene == input$ent_gene) %>% 
      pull(Clade) %>% unique() %>% sort()
    
    updateSelectInput(session, "ent_clade", choices = c("All", clade_choices), selected = "All")
  })
  observeEvent(list(input$global_subtype, input$variation_type), {
    clades <- current_usage_by_clade() %>% filter(Group == input$global_subtype) %>% pull(Clade) %>% unique() %>% sort()
    updateSelectInput(session, "lol_ref_clade", choices = clades, selected = clades[1])
    updateSelectInput(session, "lol_tar_clade", choices = clades, selected = if(length(clades)>1) clades[2] else clades[1])
  })
  
  # --- HELPER: Dynamically Update Position Limit based on Gene Length (Tab 1) ---
  observeEvent(list(input$sp_gene, input$variation_type), {
    req(input$global_subtype, input$sp_gene)
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
  # SERVER: TAB 1 - STATS (Map & Static Plots)
  # ==========================================
  
  output$total_seqs <- renderText({ paste0(total_raw) })
  output$total_countries <- renderText({ paste0(total_countries_val) })
  output$time_range <- renderText({ time_range_val })

  # --- DYNAMIC CLADE PLOT FILL DROPDOWN ---
  observeEvent(input$clade_plot_subtype, {
    req(input$clade_plot_subtype)
    
    # Find valid groups for this subtype by scanning the directory
    dir_path <- paste0("data/", input$clade_plot_subtype, "/AA/")
    files <- list.files(dir_path, pattern = "^aa_usage_by_.*\\.rds$")
    all_usage_groups <- sub("^aa_usage_by_(.*)\\.rds$", "\\1", files)
    
    valid_groups <- setdiff(all_usage_groups, c("Year", "Year_Month"))
    
    # Map usage group names to their respective metadata columns
    meta_cols <- valid_groups
    meta_cols[meta_cols == "HA_clade"] <- "clade"
    meta_cols[meta_cols == "NA_clade"] <- "G_clade"
    
    display_names <- gsub("_", " ", valid_groups)
    display_names <- gsub("clade", "Clade", display_names, ignore.case = TRUE)
    
    choices <- c(setNames(meta_cols, display_names), c("Region" = "region", "Country" = "country"))
    
    current_sel <- if (!is.null(input$clade_plot_fill) && input$clade_plot_fill %in% choices) input$clade_plot_fill else choices[1]
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
  
  sp_filtered_data <- reactive({
    subtype   <- input$global_subtype
    gene      <- input$sp_gene
    pos       <- input$sp_position
    group_col <- input$sp_group_by 
    var_type  <- input$variation_type
    
    req(subtype, gene, pos, group_col, var_type)
    
    var_lower <- tolower(var_type)
    rds_file <- paste0("data/", subtype, "/", var_type, "/", var_lower, "_usage_by_", group_col, ".rds")
    data <- get_lazy_table(rds_file)
    
    validate(need(!is.null(data), paste("Table not found for group:", group_col)))
    validate(need(group_col %in% colnames(data), "Updating data..."))
    
    # 1. Basic filtering by gene, position, subtype
    filtered <- data %>% 
      filter(Group == subtype, Gene == gene, Position == pos) %>%
      # 2. Filter out "X" and "-"
      filter(!(AminoAcid %in% c("X", "-")))
    
    # 3. Recalculate totals and frequencies based on remaining valid sequences
    # This ensures your plot bars correctly represent 100% of the available data
    filtered <- filtered %>%
      group_by(!!sym(group_col)) %>%
      mutate(
        Valid_Total = sum(Count), # This is the total count of valid AAs/NTs for the group
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
    
    return(filtered)
  })
  
  sp_plot_ggplot <- reactive({
    data <- sp_filtered_data()
    
    if(is.character(data) && grepl("Data Error", data)) return(NULL)
    req(input$sp_font_size, input$sp_group_by)
    validate(need(!is.null(data) && nrow(data) > 0, "No data available."))
    
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
        sampled_yms <- chronological_yms[seq(1, length(chronological_yms), by = 6)]
        x_scale <- scale_x_discrete(breaks = c(present_specials, sampled_yms))
      } else {
        x_scale <- scale_x_discrete(breaks = all_yms[seq(1, length(all_yms), by = 6)])
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
    req(input$sp_group_by)
    
    # Select columns to show, including Codon_Usage if present
    cols_to_show <- c(input$sp_group_by, "AminoAcid", "Count", "Valid_Total", "Frequency(%)")
    if("Codon_Usage" %in% colnames(data)) cols_to_show <- c(cols_to_show, "Codon_Usage")
    
    # Safely intersect to completely prevent "Can't subset elements that don't exist" errors
    cols_to_show <- intersect(cols_to_show, colnames(data))
    req(input$sp_group_by %in% colnames(data))
    
    datatable(
      data %>% dplyr::select(all_of(cols_to_show)) %>% arrange(!!sym(input$sp_group_by), desc(`Frequency(%)`)), 
      options = list(pageLength = 10, autoWidth = TRUE), 
      rownames = FALSE
    ) %>% formatRound("Frequency(%)", digits = 2)
  })
  
  output$sp_position_details <- renderUI({
    req(input$global_subtype, input$sp_gene, input$sp_position)
    # Only show important sites information in Amino Acid mode
    if (input$variation_type == "NT") return(NULL)
    
    match <- important_pos_df %>% 
      filter(Subtype == as.character(input$global_subtype), 
             Gene == as.character(input$sp_gene), 
             Position == input$sp_position)
    
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
  
  clicked_data_val <- reactiveValues(gene = NULL, pos = NULL)
  
  get_aggregated_clade <- function(subtype, clade_name, min_freq) {
    current_usage_by_clade() %>% 
      filter(Group == subtype, Clade == clade_name) %>% 
      # NEW: Exclude unknowns and gaps before any aggregation
      filter(!(AminoAcid %in% c("X", "-"))) %>% 
      group_by(Gene, Position, AminoAcid) %>%
      summarise(Variant_Count = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
      # Recalculate totals and frequencies based ONLY on valid sequences
      mutate(
        Total_Seqs = sum(Variant_Count),
        Freq = (Variant_Count / Total_Seqs) * 100
      ) %>%
      # Select the dominant valid amino acid
      filter(Freq == max(Freq)) %>%
      filter(row_number() == 1, Freq >= min_freq) %>%
      ungroup()
  }
  
  pairwise_differences <- reactive({
    req(input$pw_clade1, input$pw_clade2, input$global_subtype)
    if(input$pw_clade1 == input$pw_clade2) return(data.frame())
    
    # Fetch dominant AA for Clade 1 (Cleaned of X/-)
    c1_dom <- get_aggregated_clade(input$global_subtype, input$pw_clade1, input$pw_min_freq) %>%
      dplyr::select(Gene, Position, Clade1_AA = AminoAcid, Clade1_Freq = Freq)
    
    # Fetch dominant AA for Clade 2 (Cleaned of X/-)
    c2_dom <- get_aggregated_clade(input$global_subtype, input$pw_clade2, input$pw_min_freq) %>%
      dplyr::select(Gene, Position, Clade2_AA = AminoAcid, Clade2_Freq = Freq)
    
    # Join and filter for actual substitutions (e.g., A -> V)
    inner_join(c1_dom, c2_dom, by = c("Gene", "Position")) %>% 
      filter(Clade1_AA != Clade2_AA) %>% 
      arrange(Gene, Position)
  })
  
  output$pw_diff_table <- renderDT({
    data <- pairwise_differences()
    if(nrow(data) == 0) return(datatable(data.frame(Message = "No robust differences found."), rownames = FALSE))
    
    display_data <- data %>% 
      mutate(Position = sprintf('<a href="#" onclick="Shiny.setInputValue(\'modal_clicked\', \'%s|%s\', {priority: \'event\'});"><strong>%s</strong></a>', Gene, Position, Position)) %>% 
      dplyr::rename(`Clade 1 AA` = Clade1_AA, `Clade 1 Freq (%)` = Clade1_Freq, `Clade 2 AA` = Clade2_AA, `Clade 2 Freq (%)` = Clade2_Freq)
    
    datatable(display_data, escape = FALSE, options = list(pageLength = 15, autoWidth = TRUE), rownames = FALSE) %>% 
      formatRound(c("Clade 1 Freq (%)", "Clade 2 Freq (%)"), digits = 2)
  })
  
  observeEvent(input$modal_clicked, {
    parts <- strsplit(input$modal_clicked, "\\|")[[1]]
    clicked_data_val$gene <- parts[1]
    clicked_data_val$pos <- as.numeric(parts[2])
    
    showModal(modalDialog(
      title = paste(variant_label(), "Usage: Gene", clicked_data_val$gene, "- Position", clicked_data_val$pos),
      size = "l", easyClose = TRUE,
      fluidRow(
        column(4, sliderInput("modal_font_size", "Plot Font Size:", min = 10, max = 24, value = 14, step = 1)),
        column(4, radioButtons("modal_plot_format", "Format:", choices = c("PDF", "PNG"), inline = TRUE)),
        column(4, downloadButton("downloadModalPlot", "Download Plot", class = "btn-info", style="margin-top: 25px; width: 100%;"))
      ),
      plotlyOutput("modal_plot", height = "500px"), 
      hr(), 
      DTOutput("modal_table")
    ))
  })
  
  modal_data <- reactive({ 
    req(clicked_data_val$gene, clicked_data_val$pos)
    
    # Capture variation type once
    is_aa <- (input$variation_type == "AA")
    
    res <- current_usage_by_clade() %>% 
      filter(Group == input$global_subtype, Gene == clicked_data_val$gene, Position == clicked_data_val$pos) %>%
      # NEW: Remove "X" and "-" before any calculations
      filter(!(AminoAcid %in% c("X", "-")))
      
    if (nrow(res) == 0) return(res)
    
    # Check if Codon_Usage column exists before pipe to avoid dot-reference issues
    has_codon <- "Codon_Usage" %in% colnames(res)
    
    res <- res %>%
      group_by(Clade, AminoAcid) %>%
      summarise(
        Count = sum(Count, na.rm = TRUE), 
        Codon_Usage = if(has_codon) dplyr::first(Codon_Usage) else NA_character_,
        .groups = "drop_last"
      ) %>%
      # Recalculate the denominator based ONLY on valid amino acids
      mutate(Total_in_Clade = sum(Count)) %>%
      mutate(`Frequency(%)` = (Count / Total_in_Clade) * 100) %>%
      ungroup()
      
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
          "Clade: ", Clade, 
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
    req(input$modal_font_size, input$pw_clade1, input$pw_clade2)
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
    
    ggplot(data, aes(x = Clade, y = !!sym("Frequency(%)"), fill = AminoAcid, group = AminoAcid, 
                     text = tooltip_text)) + 
      geom_col(color = "black", size = 0.2) + 
      scale_fill_manual(values = current_colors(), drop = FALSE) + 
      scale_y_continuous(expand = c(0,0), limits = c(0, 105)) + 
      scale_x_discrete(labels = html_labels) + 
      labs(x = "Clade", y = "Frequency (%)", fill = variant_label()) + 
      theme_minimal(base_size = input$modal_font_size) + 
      theme(
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1), 
        axis.title = element_text(face = "bold"), 
        panel.grid.major.x = element_blank()
      )
  })
  
  output$modal_plot <- renderPlotly({
    ggplotly(modal_plot_ggplot(), tooltip = "text") %>%
      config(displayModeBar = FALSE)
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
  
  output$modal_table <- renderDT({ 
    data <- modal_data()
    cols_to_show <- c("Clade", "AminoAcid", "Count", "Total_in_Clade", "Frequency(%)")
    if("Codon_Usage" %in% colnames(data)) cols_to_show <- c(cols_to_show, "Codon_Usage")
    
    datatable(data %>% dplyr::select(all_of(cols_to_show)) %>% arrange(Clade, desc(`Frequency(%)`)), 
              options = list(pageLength = 5, autoWidth = TRUE), rownames = FALSE) %>% formatRound("Frequency(%)", digits = 2) 
  })
  
  output$downloadPairwiseCSV <- downloadHandler(
    filename = function() { paste0("Differences_", input$pw_clade1, "_vs_", input$pw_clade2, ".csv") },
    content = function(file) { 
      data <- pairwise_differences()
      if(nrow(data)>0) data <- data %>% mutate(Group = input$global_subtype) %>% dplyr::select(Group, everything())
      write_csv(data, file) 
    }
  )
  
  output$downloadPairwiseExcel <- downloadHandler(
    filename = function() { paste0("Matrices_", input$pw_clade1, "_vs_", input$pw_clade2, ".xlsx") },
    content = function(file) {
      diffs <- pairwise_differences(); wb <- createWorkbook()
      if(nrow(diffs) == 0) { 
        addWorksheet(wb, "No Differences"); writeData(wb, "No Differences", "No differences found."); saveWorkbook(wb, file, overwrite = TRUE); return() 
      }
      base_df <- data.frame(AminoAcid = if(input$variation_type == "AA") ALL_AAS else c("a","c","g","t","A","C","G","T","N","n","-"))
      for(i in 1:nrow(diffs)) {
        r_gene <- diffs$Gene[i]; r_pos <- diffs$Position[i]
        
        pos_data <- current_usage_by_clade() %>% 
          filter(Group == input$global_subtype, Gene == r_gene, Position == r_pos) %>%
          # Step 1: Remove "X" and "-" before any calculations
          filter(!(AminoAcid %in% c("X", "-"))) %>%
          group_by(Clade, AminoAcid) %>% 
          summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop_last") %>%
          # Step 2: Recalculate Frequency based only on the sum of valid amino acids per Clade
          mutate(`Frequency(%)` = (Count / sum(Count)) * 100) %>% 
          ungroup()
        
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
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ==========================================
  # SERVER: TAB 3 - ENTROPY LANDSCAPE
  # ==========================================
  output$ent_plot_title <- renderText({ 
    clade_text <- if(input$ent_clade == "All") "All Clades" else paste("Clade", input$ent_clade)
    mode_text <- if(input$variation_type == "AA") "Amino Acid" else "Nucleotide"
    paste(mode_text, "Shannon Entropy Landscape - Subtype", input$global_subtype, "| Gene", input$ent_gene, "|", clade_text) 
  })
  
  output$ent_plot <- renderPlotly({
    req(input$global_subtype, input$ent_gene, input$ent_clade)
    
    tmp <- current_usage_by_clade() %>% 
      filter(Group == input$global_subtype, Gene == input$ent_gene)
    
    if (input$ent_clade != "All") {
      tmp <- tmp %>% filter(Clade == input$ent_clade)
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
    paste(mode_text, "Mutations in", input$lol_tar_clade, "vs Reference", input$lol_ref_clade, "(Gene", input$lol_gene, ")") 
  })
  
  lol_plot_object <- reactive({
    req(input$global_subtype, input$lol_gene, input$lol_ref_clade, input$lol_tar_clade)
    validate(need(input$lol_ref_clade != input$lol_tar_clade, "Reference and Target clades are identical. Please select different clades."))
    
    # 1. Fetch Reference Clade Data (c1)
    c1 <- current_usage_by_clade() %>% 
      filter(Group == input$global_subtype, Clade == input$lol_ref_clade, Gene == input$lol_gene) %>% 
      # Step A: Exclude "X" and "-"
      filter(!(AminoAcid %in% c("X", "-"))) %>%
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
    c2 <- current_usage_by_clade() %>% 
      filter(Group == input$global_subtype, Clade == input$lol_tar_clade, Gene == input$lol_gene) %>% 
      # Step A: Exclude "X" and "-"
      filter(!(AminoAcid %in% c("X", "-"))) %>%
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
    
    muts <- inner_join(c1, c2, by = "Position") %>% 
      filter(Ref_AA != Tar_AA) %>% 
      arrange(Position) %>%
      mutate(
        Label = paste0(Ref_AA, Position, Tar_AA),
        Y_Level = rep(c(1.0, 1.4, 1.8, 2.2), length.out = n()),
        HoverText = paste("Position:", Position, "<br>Mutation:", Label)
      )
    
    validate(need(nrow(muts) > 0, "No fixed mutations found between these clades."))
    
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
    filename = function() { paste0("Lollipop_", input$lol_ref_clade, "_vs_", input$lol_tar_clade, "_", input$lol_gene, ".", tolower(input$lol_plot_format)) },
    content = function(file) { ggsave(file, plot = lol_plot_object(), device = tolower(input$lol_plot_format), width = 12, height = 6, dpi = 300) }
  )
  
  # ==========================================
  # SERVER: TAB 5 - CONSENSUS HEATMAP (msaR)
  # ==========================================
  output$heat_plot_title <- renderText({ 
    mode_text <- if(input$variation_type == "AA") "Amino Acid" else "Nucleotide"
    paste("Interactive Consensus", mode_text, "MSA - Subtype", input$global_subtype, "| Gene", input$heat_gene) 
  })
  
  output$msa_dynamic_container <- renderUI({
    req(input$global_subtype)
    clade_count <- current_usage_by_clade() %>% filter(Group == input$global_subtype) %>% pull(Clade) %>% unique() %>% length()
    
    if (!is.null(input$show_mut_only) && input$show_mut_only) {
      clade_count <- clade_count + 1
    }
    
    outer_height <- (clade_count * 20) + 150
    msaROutput("heat_plot", width = "100%", height = paste0(outer_height, "px"))
  })
  
  output$heat_plot <- renderMsaR({
    req(input$global_subtype, input$heat_gene)
    
    cache_dir <- "data"
    if (!dir.exists(cache_dir)) dir.create(cache_dir, showWarnings = FALSE)
    
    safe_subtype <- gsub("[^A-Za-z0-9_]", "_", input$global_subtype)
    safe_gene <- gsub("[^A-Za-z0-9_]", "_", input$heat_gene)
    
    prefix <- if(input$variation_type == "AA") "MSA_AA_" else "MSA_NT_"
    # Filename no longer includes minFreq as it is removed
    aln_filename <- file.path(cache_dir, paste0(prefix, safe_subtype, "_", safe_gene, ".fasta"))
    
    if (file.exists(aln_filename)) {
      if(input$variation_type == "AA") {
        aligned_strings <- Biostrings::readAAStringSet(aln_filename)
      } else {
        aligned_strings <- Biostrings::readDNAStringSet(aln_filename)
      }
      original_clade_order <- names(aligned_strings) 
      
    } else {
      # 1. Prepare Base Grid of All Positions and Clades to ensure NO GAPS
      all_pos_in_gene <- current_usage_by_clade() %>% 
        filter(Group == input$global_subtype, Gene == input$heat_gene) %>% 
        pull(Position) %>% unique() %>% sort()
      
      all_clades_in_group <- current_usage_by_clade() %>% 
        filter(Group == input$global_subtype, Gene == input$heat_gene) %>% 
        pull(Clade) %>% unique() %>% sort()
      
      grid <- expand.grid(Clade = all_clades_in_group, Position = all_pos_in_gene, stringsAsFactors = FALSE)
      
      # 2. Filter and Identify Majority Character per Position/Clade
      exclude_chars <- if(input$variation_type == "AA") c("X", "-") else c("N", "n", "-")
      
      raw_data <- current_usage_by_clade() %>% 
        filter(Group == input$global_subtype, Gene == input$heat_gene) %>%
        # Filter out ambiguous/gaps for calculation
        filter(!(AminoAcid %in% exclude_chars)) %>% 
        group_by(Clade, Position, AminoAcid) %>%
        summarise(Total_Count = sum(Count), .groups = "drop") %>%
        group_by(Clade, Position) %>%
        filter(Total_Count == max(Total_Count)) %>%
        filter(row_number() == 1) %>%
        ungroup()
      
      # 3. Merge with Grid and Fill Gaps with "-"
      complete_data <- grid %>%
        left_join(raw_data, by = c("Clade", "Position")) %>%
        mutate(AminoAcid = ifelse(is.na(AminoAcid), "-", AminoAcid)) %>%
        arrange(Clade, Position)
      
      # 4. Reconstruct the consensus sequences
      raw_seqs <- complete_data %>%
        group_by(Clade) %>%
        summarise(seq = toupper(paste(AminoAcid, collapse = "")), .groups = "drop") %>%
        arrange(Clade)
      
      original_clade_order <- raw_seqs$Clade
      
      if(input$variation_type == "AA") {
        unaligned_strings <- Biostrings::AAStringSet(setNames(raw_seqs$seq, original_clade_order))
      } else {
        unaligned_strings <- Biostrings::DNAStringSet(setNames(raw_seqs$seq, original_clade_order))
      }
      
      waiter_show(
        html = tagList(
          spin_fading_circles(), 
          h3(paste("Running ClustalW Alignment for", input$heat_gene, "sequences..."), style = "color: white; margin-top: 20px;")
        ),
        color = "rgba(44, 62, 80, 0.9)"
      )
      
      on.exit(waiter_hide(), add = TRUE) 
      
      aligned_msa <- suppressMessages(msa::msa(unaligned_strings, method = "ClustalW", order = "input"))
      
      if(input$variation_type == "AA") {
        aligned_strings <- as(aligned_msa, "AAStringSet")
      } else {
        aligned_strings <- as(aligned_msa, "DNAStringSet")
      }
      aligned_strings <- aligned_strings[original_clade_order]
      
      Biostrings::writeXStringSet(aligned_strings, filepath = aln_filename)
      
    } 
    
    if (input$show_mut_only) {
      seq_char_matrix <- as.matrix(aligned_strings)
      exclude_chars <- if(input$variation_type == "AA") c("X", "-") else c("N", "n", "-")
      
      consensus_seq <- apply(seq_char_matrix, 2, function(col) {
        # Filter out X/- or N/- before finding the majority
        valid_col <- col[!(toupper(col) %in% toupper(exclude_chars))]
        if(length(valid_col) == 0) return("-") # Fallback if all are ambiguous
        freqs <- table(valid_col)
        names(freqs)[which.max(freqs)] 
      })
      
      for (i in 1:nrow(seq_char_matrix)) {
        match_idx <- seq_char_matrix[i, ] == consensus_seq
        seq_char_matrix[i, match_idx] <- "."
      }
      
      consensus_string <- paste(consensus_seq, collapse = "")
      new_seqs <- apply(seq_char_matrix, 1, paste, collapse = "")
      
      if(input$variation_type == "AA") {
        final_strings <- Biostrings::AAStringSet(c(Consensus = consensus_string, new_seqs[original_clade_order]))
      } else {
        final_strings <- Biostrings::DNAStringSet(c(Consensus = consensus_string, new_seqs[original_clade_order]))
      }
      
    } else {
      final_strings <- aligned_strings
    }
    
    inner_align_height <- (length(final_strings) * 20) + 20
    
    msaR(
      final_strings, 
      menu = TRUE, 
      overviewbox = FALSE, 
      seqlogo = !isTRUE(input$show_mut_only), 
      colorscheme = if(input$variation_type == "AA") "clustal" else "nucleotide",
      alignmentHeight = inner_align_height
      )
      })

      # --- HIDE LOADING CURTAIN WHEN READY ---
      session$onFlushed(function() {
        waiter_hide()
      }, once = TRUE)
      }