# ==========================================
# 2. UI (User Interface)
# ==========================================

library(shinyjs)

ui <- navbarPage(
  
  title = div(
    tags$img(src = "app_icon_round.png", height = "30px", style = "margin-right: 10px; vertical-align: middle;"),
    "FLU Divergence Explorer"
  ),
  
  header = tags$head(
    tags$link(rel = "shortcut icon", href = "app_icon_round.png"), 
    use_waiter(), 
    shinyjs::useShinyjs(),
    tags$style(HTML("
      @media (min-width: 768px) { 
        .modal-dialog { width: 80vw !important; max-width: 80vw !important; } 
      }
      .jumbotron { background-color: #f8f9fa; padding: 2rem 2rem; border-radius: 0.5rem; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-top: 20px; }
      .navbar-brand { padding-top: 10px; padding-bottom: 10px; height: auto; }
      .tab-content { min-height: calc(100vh - 160px); }
      
      /* Fixed Navbar Layout Fix */
      .navbar-header { float: left; }
      .navbar-nav { float: left; }
      
      /* Significant Switch Styling - FIXED position */
      .variation-switch-container {
        position: fixed;
        right: 20px;
        top: 8px;
        z-index: 2000; 
        background: #ffffff;
        padding: 5px 20px;
        border-radius: 40px;
        border: 2px solid #2c3e50;
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        display: flex;
        align-items: center; 
        gap: 15px;
      }
      .variation-switch-container .btn-group-container-sw {
        display: flex;
        align-items: center;
        gap: 10px;
      }
      .switch-label {
        font-weight: 800;
        color: #2c3e50;
        font-size: 0.85em;
        text-transform: uppercase;
        letter-spacing: 1px;
        margin-bottom: 0;
      }
      /* Custom Picker Styling */
      .variation-switch-container .bootstrap-select .btn {
        background-color: #3498db !important;
        color: white !important;
        border-radius: 20px !important;
        font-weight: bold !important;
        border: none !important;
      }
      
      /* Custom Radio Button Styling */
      .variation-switch-container .btn-default {
        background-color: #f1f3f5 !important;
        color: #495057 !important;
        border: 1px solid #ced4da !important;
        font-weight: 700 !important;
        transition: all 0.2s ease-in-out;
        padding: 5px 12px;
      }
      .variation-switch-container .btn-default.active {
        background-color: #e74c3c !important; /* Bold color for selection */
        color: #ffffff !important;
        border-color: #c0392b !important;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.2) !important;
      }
      .variation-switch-container .btn-default:hover:not(.active) {
        background-color: #e9ecef !important;
      }
      /* Remove potential padding from the widget itself */
      .variation-switch-container .form-group {
        margin-bottom: 0 !important;
      }
    ")),
    # Fixed positioned container for the switch
    div(class = "variation-switch-container",
        span(class = "switch-label", "Subtype"),
        pickerInput(
          inputId = "global_subtype",
          label = NULL,
          choices = sort(na.omit(unique(metadata_global$Group))),
          selected = sort(na.omit(unique(metadata_global$Group)))[1],
          width = "120px"
        ),
        div(style = "width: 1px; height: 25px; background: #ced4da;"),
        div(class = "btn-group-container-sw",
            span(class = "switch-label", "Mode"),
            radioGroupButtons(
              inputId = "variation_type",
              label = NULL,
              choices = c("AA", "NT"),
              selected = "AA",
              status = "default",
              size = "sm"
            )
        )
    )
  ),
  
  footer = tags$footer(
    style = "text-align: center; padding: 15px; background-color: #f8f9fa; border-top: 1px solid #e7e7e7; color: #6c757d; margin-top: 30px; width: 100%;",
    HTML(paste0("&copy; ", format(Sys.Date(), "%Y"), " FLU Amino Acid Divergence Explorer. All rights reserved."))
  ),
  
  # ---------------------------------------------------------
  # TAB 0: HOME
  # ---------------------------------------------------------
  tabPanel("Home",
           fluidPage(
             div(class = "jumbotron",
                 h1("Welcome to the FLU Amino Acid Divergence Explorer", style = "color: #2c3e50; font-weight: bold;"),
                 p("A high-resolution visualization tool for analyzing Influenza Virus (H1N1 and H3N2) genetic diversity across lineages.", style = "font-size: 1.2em; color: #7f8c8d;"),
                 hr(),
                 div(style = "text-align: center; margin-top: 20px; margin-bottom: 30px;",
                     img(src = "welcome_banner.png", style = "max-width: 100%; height: auto; border-radius: 10px; box-shadow: 0 4px 8px rgba(0,0,0,0.2);")
                 ),
                 h3("How to Use This App:", style = "color: #2980b9;"),
                 fluidRow(
                   column(4, h4(icon("chart-bar"), " Single Position Explorer"), p("Dive deep into the amino acid or nucleotide distribution of any specific position within an Influenza gene (HA, NA, etc.).")),
                   column(4, h4(icon("not-equal"), " Pairwise Comparison"), p("Instantly identify robust, fixed differences between any two Influenza clades across all genes.")),
                   column(4, h4(icon("globe"), " Gene-Wide Landscapes"), p("Explore whole-gene visualizations including Entropy conservation plots, Lollipop mutation trackers, and Consensus Alignments."))
                 )
             )
           )
  ),
  
  # TAB 1: DATASET STATS (World Map + Static Plots)
  tabPanel("Dataset Insights",
           fluidPage(
             fluidRow(
               column(4, wellPanel(style = "text-align: center; background-color: #f8f9fa; border: 1px solid #ddd;", 
                                   h4("Total Sequences", style="color: #2c3e50;"), h2(textOutput("total_seqs"), style="color: #2980b9; font-weight: bold;"))),
               column(4, wellPanel(style = "text-align: center; background-color: #f8f9fa; border: 1px solid #ddd;", 
                                   h4("Countries Represented", style="color: #2c3e50;"), h2(textOutput("total_countries"), style="color: #2980b9; font-weight: bold;"))),
               column(4, wellPanel(style = "text-align: center; background-color: #f8f9fa; border: 1px solid #ddd;", 
                                   h4("Time Span", style="color: #2c3e50;"), h2(textOutput("time_range"), style="color: #2980b9; font-weight: bold;")))
             ),
             
             wellPanel(
               fluidRow(
                 column(2, selectInput("map_geo_level", "Grouping:", choices = c("Region", "Country"))),
                 column(2, selectInput("map_clade_type", "Pie Data:", choices = c("HA-Clade" = "clade", "NA-Clade" = "G_clade"))),
                 column(3, selectInput("map_year", "Select Year:", choices = c("All", sort(na.omit(unique(metadata_global$Year)), decreasing = TRUE)))),
                 column(5, style = "margin-top: 25px;", helpText("Global FLU Distribution Analysis. Subtype is controlled globally from the top right."))
               )
             ),
             
             fluidRow(
               column(12, 
                      h4("Global Clade Distribution", style="font-weight: bold; color: #2c3e50;"),
                      leafletOutput("world_map", height = "500px"),
                      hr())
             ),
             
             fluidRow(
               column(6, 
                      h4("Sequencing Over Time (Seasonality)", style="font-weight: bold; margin-top: 20px;"),
                      plotlyOutput("stats_time_plot", height = "400px")
               ),
               column(6, 
                      h4("Regional Breakdown", style="font-weight: bold; margin-top: 20px;"),
                      plotlyOutput("stats_geo_plot", height = "400px")
               )
             ),
             
             fluidRow(
               column(12, hr()),
               column(4, selectInput("clade_plot_x", "Primary Category (X-axis):", 
                                     choices = c("HA-Clade" = "clade", "NA-Clade" = "G_clade", "Year" = "Year", "Region" = "region", "Country" = "country"),
                                     selected = "Year")),
               column(4, selectInput("clade_plot_fill", "Sub-Category (Color):", 
                                     choices = c("HA-Clade" = "clade", "NA-Clade" = "G_clade", "Year" = "Year", "Region" = "region", "Country" = "country"),
                                     selected = "clade")),
               column(12, 
                      h4("Custom Dataset Breakdown", style="font-weight: bold; margin-top: 10px;"),
                      plotlyOutput("stats_clade_plot", height = "500px")
               )
             )
           )
  ),
  
  # ---------------------------------------------------------
  # MACRO-LEVEL DROPDOWN MENU
  # ---------------------------------------------------------
  navbarMenu("Gene-Wide Landscapes",
             
             tabPanel("Conservation (Entropy)",
                      fluidPage(
                        wellPanel(
                          fluidRow(
                            column(4, helpText("Calculates Shannon Entropy to identify highly conserved valleys and hypervariable peaks across the entire gene. Subtype is controlled globally.")),
                            column(4, selectInput("ent_gene", "Gene:", choices = NULL)),
                            column(4, selectInput("ent_clade", "Clade:", choices = NULL))
                          ),
                          fluidRow(
                            column(3, sliderInput("ent_min_seqs", "Min Sequences:", min = 0, max = 1000, value = 10, step = 10)),
                            column(3, sliderInput("ent_font_size", "Plot Font Size:", min = 10, max = 24, value = 14, step = 1)),
                            column(3, radioButtons("ent_plot_format", "Format:", choices = c("PNG", "PDF"), inline = TRUE)),
                            column(3, downloadButton("downloadEntPlot", "Download Plot", class = "btn-info", style="margin-top: 25px; width: 100%;"))
                          )
                        ),
                        h3(textOutput("ent_plot_title")),
                        plotlyOutput("ent_plot", height = "450px") 
                      )
             ),
             
             tabPanel("Mutation Tracker (Lollipop)",
                      fluidPage(
                        wellPanel(
                          fluidRow(
                            column(4, helpText("Visualize fixed amino acid mutations in a Target Clade compared to a Reference Clade. Subtype is controlled globally.")),
                            column(4, selectInput("lol_gene", "Gene:", choices = NULL)),
                            column(4, numericInput("lol_min_freq", "Min Dominant Freq (%):", value = 90.0, min = 50.0, max = 100.0))
                          ),
                          fluidRow(
                            column(3, selectInput("lol_ref_clade", "Reference Clade:", choices = NULL)),
                            column(3, selectInput("lol_tar_clade", "Target Clade:", choices = NULL)),
                            column(2, sliderInput("lol_font_size", "Font Size:", min = 10, max = 24, value = 14, step = 1)),
                            column(2, radioButtons("lol_plot_format", "Format:", choices = c("PNG", "PDF"), inline = TRUE)),
                            column(2, downloadButton("downloadLolPlot", "Download Plot", class = "btn-info", style="margin-top: 25px; width: 100%;"))
                          )
                        ),
                        h3(textOutput("lol_plot_title")),
                        plotlyOutput("lol_plot", height = "550px") 
                      )
             ),
             
             # TAB 5: CONSENSUS MSA (FULL-WIDTH LAYOUT)
             tabPanel("Consensus MSA Map",
                      fluidPage(
                        wellPanel(
                          fluidRow(
                            column(4, helpText("Interactive Multiple Sequence Alignment. Subtype is controlled globally.")),
                            column(4, selectInput("heat_gene", "Gene:", choices = NULL)),
                            column(4, div(style = "margin-top: 25px;", checkboxInput("show_mut_only", "Show Mutations Only", value = FALSE)))
                          )
                        ),
                        fluidRow(
                          column(12,
                                 h3(textOutput("heat_plot_title")),
                                 uiOutput("msa_dynamic_container") 
                          )
                        )
                      )
             )
  ),
  
  # ---------------------------------------------------------
  # TAB 1: SINGLE POSITION
  # ---------------------------------------------------------
  tabPanel("Single Position Explorer",
           sidebarLayout(
             sidebarPanel(
               h5("Setting", style="font-weight: bold; color: #2980b9;"),
               selectInput("sp_group_by", "Group by:", choices = c("Clade", "Year", "Year-Month" = "Year_Month")),
               checkboxInput("sp_show_counts", "Show raw counts instead of percentage", value = FALSE),
               hr(),
               div(id = "sp_quick_access_section",
                 h5("Quick Access", style="font-weight: bold; color: #2980b9;"),
                 selectInput("sp_quick_visit", "Jump to Key Position:", 
                             choices = if(nrow(important_pos_df) > 0) c("Manual Selection" = "None", setNames(1:nrow(important_pos_df), important_pos_df$label)) else c("Manual Selection" = "None")),
                 hr()
               ),
               h5("Precise Access", style="font-weight: bold; color: #2980b9;"),
               selectInput("sp_gene", "Gene:", choices = NULL),
               
               uiOutput("sp_range_label"),
               div(style = "display: flex; align-items: center; gap: 5px; margin-bottom: 15px;",
                   actionButton("sp_pos_minus", "-", class = "btn-primary", style = "padding: 6px 12px; font-weight: bold; height: 34px;"),
                   div(style = "width: 80px;margin-top: 15px", 
                       tags$style(HTML("#sp_position { margin-bottom: 0px !important; height: 34px; text-align: center; }")),
                       numericInput("sp_position", label = NULL, value = 1, min = 1, max = 1000)
                   ),
                   actionButton("sp_pos_plus", "+", class = "btn-primary", style = "padding: 6px 12px; font-weight: bold; height: 34px;")
               ),
               sliderInput("sp_min_seqs", "Min Seqs:", min = 1, max = 500, value = 5),
               sliderInput("sp_font_size", "Plot Font Size:", min = 10, max = 24, value = 14),
               hr(),
               radioButtons("sp_plot_format", "Download Format:", choices = c("PDF", "PNG"), inline = TRUE),
               downloadButton("downloadSpPlot", "Download Plot", class = "btn-info", style="width: 100%;")
             ),
             mainPanel(
               uiOutput("sp_position_details"),
               plotlyOutput("sp_aa_plot", height = "500px"),
               DTOutput("sp_aa_table")
             )
           )
  ),
  
  # ---------------------------------------------------------
  # TAB 2: PAIRWISE COMPARISON
  # ---------------------------------------------------------
  tabPanel("Pairwise Comparison",
           sidebarLayout(
             sidebarPanel(
               selectInput("pw_clade1", "Clade 1:", choices = NULL),
               selectInput("pw_clade2", "Clade 2:", choices = NULL),
               numericInput("pw_min_freq", "Minimum Dominant Frequency (%):", value = 90.0, min = 50.0, max = 100.0, step = 1.0),
               hr(),
               downloadButton("downloadPairwiseCSV", "Download Table (CSV)", class = "btn-primary", style="margin-bottom: 5px; width: 100%;"),
               downloadButton("downloadPairwiseExcel", "Download Excel Matrices", class = "btn-success", style="width: 100%;")
             ),
             mainPanel(
               h3("Cross-Gene Pairwise Differences"),
               p("Click on any highlighted Position to view the full amino acid distribution for that specific site."),
               DTOutput("pw_diff_table")
             )
           )
  )
)