library(shiny)
# library(shinydashboard)
library(htmltools)

library(dplyr)
library(tidyr)
library(readr)
# library(stringr)

library(igraph)
library(tidygraph)
library(ggraph)

library(sf)
library(rgdal)
library(leaflet)
# library(leaflet.esri)
library(units)
library(ggplot2)


# UI

shinyUI(
    
    # 1 Menu? Tabs?
    # menuItem("Species Distribution Maps ", tabName = "maps"),
    # menuItem("Spatial Congruence", tabName = "Cs_tab"),
    # menuItem("SCAN Analysis", tabName = "scan"),
    # menuItem("SCAN Viewer", tabName = "SCAN_viewer"),
    # menuItem("about SCAN", tabName = "about_SCAN")

    
    # 2 main body
    # CSS HTML ----

    # 2.1 MAPS ----
    # text
    # htmltools::includeMarkdown('www/scan_eng_text1.Rmd'),
    
    # input map
    fileInput( inputId = "filemap",  label = "Choose shape-files (.shp + .shx + .dbl + .prj).", 
               accept = c('.shp','.dbf','.sbn','.sbx','.shx',".prj"),width = '500px', multiple=TRUE),  
    
    # fix invalid
    column(width = 6, checkboxInput("fix_invalid_shapes", label = "Fix invalid species layers?", value = FALSE)),
    
    # CRS transform
    column(width = 6, checkboxInput("modify_crs", "Change map projection (crs)?", value = FALSE)),
    conditionalPanel(condition = "input.modify_crs == true", 
                     numericInput(inputId = "map_projection", label = "Choose a crs code to project the map", value = 4326),
                     tags$h6("WGS84(degree unit: 4326; metric unit: 3857); SIRGAS2000(4674, UTM: 31985, 31977); SAD69(4618)")
    ),
    
    # sp column
    inputPanel( textInput("colum_sp_map","Which one is the species ID?", value = "sp")),
    
    # LOAD MAP!!!     
    actionButton("get_map", "Load map!",  class = "btn-warning", color = 'black'),
    #Map names in the loaded map
    
    verbatimTextOutput( "map_upload_names" ),
    
    
    # MAP SAMPLE
    tags$h3("Map sample"),
    plotOutput("map_shp"),
    
    # Check MAP
    # tags$h3("Check map's column names, species, and projection"),
    # tags$h4("Map's column names:"),
    textOutput("map_shp_names"),
    # tags$h4("Are there invalid species?"),
    tableOutput("invalid_map_species"),
    # tags$h4("CRS / ESPG specifications"),
    textOutput("map_crs"),
    # tags$h4("All map species"),
    tableOutput("map_species"),
    
    # 2.2 Cs tab ----
    # "Spatial Congruence Index - Cs"
    # # text (write a Md)
    # "The Spatial Congruence Index (aka Hargrove Index - see Gatto & Cohn-Haft 2021) is calculated for each pair of species."
    # "Or you may write an alternative index combining the area variables 'area_overlap', 'area_sp1', 'area_sp2' (see the example in the box)."
    # alternative index
    checkboxInput(inputId = "use_alternative_index", label = "Use an alternative Cs Index", value = FALSE),
    
    conditionalPanel(condition = "input.use_alternative_index == true", 
                     
             textInput(inputId = "cs_similarity_index", label = "Choose a formula index (default = Hargrove Index)", 
                        value = '(area_overlap / area_sp1) * (area_overlap / area_sp2)' ) 
     ),
    
    actionButton("calculate_Cs","Apply Cs index to map", class = "btn-warning"),
    
    
    numericInput( inputId = "filter_Cs", label = "Minimum Cs value", value = 0.1),
    
    # Buffer polygons?
    tags$p("For large datasets, the Cs calculus may take a long time. The use of inner buffers in selected area-size-classes (quartiles) may save precious time,  avoiding marginal spatial overlaps." ),
    tags$p("SCAN algorithm uses original distribution areas but some overlaps may be discarded in the Cs table because of buffers. Be careful to not loose information!"),
    tags$p("Obs. Buffers are better suited for metric-based CRS projections. The map tab allows crs change (e.g. from WGS84:4326 'degree' unit to WGS84:3857 Pseudo-Mercator 'metre'."),

    checkboxInput(inputId = "use_buffer_map", label = "Use an internal buffer to avoid marginal overlaps? (metric-based CRS only)", value = FALSE),
    
    conditionalPanel( condition = "input.use_buffer_map == true", 
                      
            numericInput(inputId = "shrink_factor_buff", label = "Choose an internal buffer", value = 0.01),
                      
            checkboxGroupInput(inputId = "quantiles_to_buffer", label = "Input quartile areas to buffer", choices = c(1,2,3,4), selected = c(), inline = TRUE),
    ),
    
    # No? Then calculate Cs
    conditionalPanel( condition = "input.calculate_Cs > 0",   tags$h6("Do NOT change TAB while running ! This may take a long time depending on the number of species ... ")   
    ),
    
    # or upload Cs table...
    checkboxInput("Cs_upload_csv", "Upload a Cs csv file instead"),
    
    conditionalPanel( condition = "input.Cs_upload_csv == true ",
                      
            fileInput( inputId = "Cs_table", label = "Select a Cs table .csv file", accept = c("text/csv","text/comma-separated-values,text/plain",".csv" ))
    ),
    
    # Download Cs results
    downloadButton("download_Cs", "Download Cs.csv"),
    
    # Check tables
    "Cs table summary and head & tail",
    textOutput("check_Cs_tables"),
    
    # Head and Tail
    column( width = 6,  box(  tableOutput("Cs_head") ) ),
    column( width = 6,  box(  tableOutput("Cs_tail") ) ),
    
    # Check GRAPH
    column( width = 6, tableOutput("graph_nodes") ),
    column( width = 6, tableOutput("graph_edges") ),
    
    # 2.3 SCAN ----
    
    # set parameters
    column( width = 3, 
            numericInput(inputId = "resolution", label = "Resolution (interval between Ct)", value = 0.1) 
    ),
    column( width = 3, 
            numericInput(inputId = "threshold_max", label = "Max value of threshold Ct", value = 0.9, min = 0.2, max = 1) 
    ),
    column( width = 3, 
            numericInput(inputId = "threshold_min", label = "Min value of threshold Ct", value = 0.2, min = 0.05, max = 0.9) 
    ),
    column( width = 3, 
            conditionalPanel(condition = "input.filter_diameter == true",
                     numericInput(inputId = "max_diameter", label = "Choose the maximum (network) diameter", value = 15)
           )
    ),
    
    column( width = 3, 
            checkboxInput(inputId = "overlap", label = "Overlap criterion: Require overlap between all species?", value = TRUE)),
    
    column( width = 3,  
            checkboxInput(inputId = "filter_diameter", label = "Limit the diameter of the chorotypes networks?", value = TRUE)),
    
    column( width =2, 
            actionButton("run_scan", "SCAN!", class = "btn-warning") ),
    
    # Warning to wait scan
    # conditionalPanel( condition = "run_scan > 0", tags$h6("Do NOT change TAB while running ! This also may last a long time... check R console")  ),
    
    # Check parameters
    tableOutput("parameters"),
    
    uiOutput("names_scan_list"),
        
    downloadButton("downloadData", "Download"),
    
    # 'SCAN results - download data preview'
    dataTableOutput('table_download_preview'),
    
    #'All Chorotypes!'
    tableOutput("scan_chorotypes"),
    
    # 2.4 SCAN_viewer ----
    
    # viewer parameters
    
    column(width = 2,
           selectInput(inputId = "palette", label = "choose palette", choices = c("Dark2","Set1","Paired","Accent","Spectral","Greens", "Reds","BrBG" , "RdYlGn","PiYG"))
    ),
    
    column(width = 2,
           selectInput("layout", "graph layout", choices = c("fr", "kk", "dh", "drl", "mds", "gem"))
    ),
    
    column(width = 2,
           numericInput("map_alpha", "map alpha", max = 1, min = 0.01, value = 0.4)
    ),
    
    # update graph?
    column(width = 3, 
            checkboxInput("graph_from_csv", "Upload a graph (nodes.csv + edges.csv)?", FALSE) 
   ),
    
    conditionalPanel(   condition = "input.graph_from_csv == true",
                        
            fileInput(inputId = "graph_nodes", label = "graph NODES in csv file", accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                                
            fileInput(inputId = "graph_edges", label = "graph EDGES in csv",  accept = c("text/csv","text/comma-separated-values,text/plain",".csv"))
    ),
    
    # Ct slider or input?
    # column(width = 2, numericInput(inputId = "threshold", label = "Threshold (0-1)", value = 0.9, max = 1, min = 0.1)),
    
    sliderInput(inputId = "threshold",label = "Select a Congruence Threshold - Ct vaule (0-1)", 
                value = 0.6, max = 1, ticks = 0.05, min = 0, step = 0.01, width = 700 ),
    
    # choose components to show
    uiOutput("original_components"),
    
    # map plot ggplot
    box(width =12, title = "Map & Graph",
        column( width = 6, 
                plotOutput("ggplot_map")),
        column( width = 6, 
                plotOutput("graph_plot2")),
    ),
    # leaflet and graph plots
    box(width = 12, title = "Leaflet plot",  
        
            leafletOutput("map_plot", height = "500px"),
            # plotOutput("graph_plot", height = "500px"), width = 4  ),
    ), 
          
    # "Chorotypes & species"
    dataTableOutput("g_sub_table")
    

)    
    
    