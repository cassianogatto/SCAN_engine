'SCAN_engine_V0.21_24fev2023'#'SCAN_engine_V0.2_02fev2023'

# TO DO
# BIG CHALLENGES
# invent a system to assign chorotypes (communities) names and recognize the synonyms across threshold rounds in SCAN
# fix palette mismatch between leaflet and ggplot2
# show updates and progress of calculations (important)

# not so big
# Settings page to manage formula, palettes, map & graph preferences, etc... 
# filter species by name!
# ggplot2 map and graph - ok
# option to split map and graph in new panels
# download FIGURES !
# show thresholds available to view (in viewer)
# plot sample map according to quantiles (4)


# filter map species
# new tab to filter species based on distribution criteria (e.g., position, area quantiles)
# check box for diameter filter - ok
# write an introduction with instructions
# write the theoretic background ('about SCAN') etc
# 

# HTML('<p> this is p </p>'),
# HTML('<code> this is code </code>'),


library(shiny)
library("shinydashboard")
library(dplyr)
library(igraph)
library(tidygraph)
library(tidyr)
library(ggraph)
library(readr)
library(sf)
library(ggplot2)
library(leaflet)
library(rgdal)
library(units)
library(leaflet.esri)

# UI

shinyUI(

    dashboardPage( #skin = "green",
        
        dashboardHeader(title = "SCAN engine V_0.21"),
        
        dashboardSidebar( width = '230px',
                          
                          sidebarMenu(
                              
                              menuItem("Species Distribution Maps ", tabName = "maps"),
                              
                              menuItem("Spatial Congruence", tabName = "Cs_tab"),
                              
                              menuItem("SCAN Analysis", tabName = "scan"),
                              
                              menuItem("SCAN Viewer", tabName = "SCAN_viewer"),
                              
                              menuItem("about SCAN", tabName = "about_SCAN")
                              
                          ) # menu items
        ),
        
        dashboardBody(
            # CSS HTML configs ----
            tags$head( 
                tags$style( 
                    HTML( 
                        #
                        #.shiny-input-container { color: #474747;   }
                    "   .content {  padding-top: 40px; padding-left: 40px;}
                        .body { background-color: white;   color: black;   }
                        .sidebar {  position:fixed;  overflow:visible; }            
                        .main-header { position: fixed; width = 100%; }
                    "
                    )                            
                )
            ),
            
            tabItems(
                # maps ----
                tabItem("maps",
                        
                    fluidPage(
                            
                            tags$h2("Map of Species' Distributions"),
                            
                            box(width = 12,
                                
                                tags$h4("Upload the shapefile map containing all species' distributions as unique layers with unique IDs and geometries.
                                         All SCAN steps are derived from these maps: Spatial similarities (Cs), Scanning algorithm, and Plotting results."),
                                
                                #HTML("<br>"),
                                
                                tags$h5("Sometimes, individual species' layers may be corrupted, for many reasons. Sometimes a projection transformation 
                                        (CRS of the whole archive) may be enough: just mark the checkbox and choose another EPSG number (see https://epsg.io). 
                                        Sometimes SCAN_eng may be used to apply sf::st_make_valid to fix individual species issues."),
                                
                                tags$h5("Maps must have, at least, one ID column named 'sp', 
                                        and a 'geometry' column. If not, check current names below
                                        and choose the ID column before loading the the map again. 
                                        Information about species and geographical projections are also shown."),
                            ),
                            
                            box(width = 6,    
                                
                                tags$h3("Input Map"),
                                     
                                fileInput( inputId = "filemap",  label = "Choose shape-files (.shp + .shx + .dbl + .prj).", 
                                                
                                                accept = c('.shp','.dbf','.sbn','.sbx','.shx',".prj"),width = '500px', multiple=TRUE),  
                                
                                
                                checkboxInput("fix_invalid_shapes", label = "Fix invalid species layers?", value = FALSE),
                                
                                checkboxInput("modify_crs", "Change map projection (crs)?", value = FALSE),
                                
                                conditionalPanel(condition = "input.modify_crs == true", 
                                                 
                                                 numericInput(inputId = "map_projection", label = "Choose a crs code to project the map", value = 4326),
                                                 
                                                 tags$h6("WGS84(degree unit: 4326; metric unit: 3857); SIRGAS2000(4674, UTM: 31985, 31977); SAD69(4618)")
                                ),
                                
                                
                                HTML('<br>'),
                                
                                tags$h4("Column names"),
                                
                                tags$h5("Shapefile data is encoded in the column 'geometry'.
                                        The column identifying species' IDs, if not 'sp' 
                                        (default), must be written it in the box below 
                                        (current column names are shown below). 
                                        Press 'Load map' again."),
                                
                                inputPanel( textInput("colum_sp_map","Which one is the species ID?", value = "sp")),
                                
                                # LOAD MAP!!!     
                                actionButton("get_map", "Load map!", class = "btn-warning"),
                                
                                tags$h4("Map names in the loaded map:"),
                                
                                verbatimTextOutput( "map_upload_names" ),
                            ),
                            
                            box( width = 6, # map sample ----
                                
                                tags$h3("Map sample"),
                                
                                plotOutput("map_shp")
                            ),
                                
                            box(width = 12,
                                
                                tags$h3("Check map's column names, species, and projection"), 
                                    
                                textOutput("map_shp_names"),
                                
                                HTML('<br>'), HTML('<br>'),
                                
                                tags$h4("Are there invalid species?"),
                                
                                tableOutput("invalid_map_species"),
                                
                                HTML('<br>'), HTML('<br>'),
                                
                                tags$h4("CRS / ESPG specifications"),
                                
                                textOutput("map_crs"),
                                
                                HTML('<br>'), HTML('<br>'),
                                
                                tags$h4("All map species"),
                                
                                tableOutput("map_species"),
                                
                                # checkboxInput('filter_map_species',"Filter species to print? (does not affect the map for analysis)"),
                                # conditionalPanel(condition = "input.filter_map_species == true", textInput(inputId = "filter_species_list", "Use the names as in the map between quotes (e.g. 'Musa paradisiaca', or c('Homo sapiens', 'Homo naledi')" )                  )
                            ),
                        )
                ),
                
                tabItem("Cs_tab",     
              
                    fluidPage(
                        
                        tags$h2("Spatial Congruence Index - Cs"),
                        
                        tags$h3("How much spatially congruent are two species?"),
                        
                        # Cs  calc ----
                        box( width = 6,
                        
                        tags$h4("The Spatial Congruence Index (aka Hargrove Index - see Gatto & Cohn-Haft 2021) is calculated for each pair of species."),
                        
                        tags$h5("Or you may write an alternative index combining the area variables 'area_overlap', 'area_sp1', 'area_sp2' (see the example in the box)."),
                        
                        checkboxInput(inputId = "use_alternative_index", label = "Use an alternative Cs Index", value = FALSE),
                        
                        fluidRow(
                            
                            conditionalPanel(condition = "input.use_alternative_index == true",
                                
                                 textInput(inputId = "cs_similarity_index", label = "Choose a formula index (default = Hargrove Index)", 
                                                
                                           value = '(area_overlap / area_sp1) * (area_overlap / area_sp2)') ),
                            
                            column( width = 4, actionButton("calculate_Cs","Apply Cs index to map", class = "btn-warning") ),
                            
                            column( width = 8, numericInput( inputId = "filter_Cs", label = "Minimum Cs value", value = 0.1))
                        ),
                    ), 
                        
                        # bufer ----    
                        box(width = 6,
                            
                            fluidRow(
                                
                                column(width = 1, ""),
                                
                                column(width = 2, icon ('skull-crossbones')),   
                                
                                column(width = 6, tags$h4("Shrinking Distribution Buffer!")),
                                
                                column(width = 2, icon ('skull-crossbones'))
                            ),
                            
                            tags$h5("For large datasets, the Cs calculus may take a long time. 
                                The use of inner buffers in selected area-size-classes (quartiles) may save precious time, 
                                avoiding marginal spatial overlaps."),
                            
                            tags$h4("SCAN algorithm uses original distribution areas but some overlaps may be discarded in the Cs table because of buffers. 
                                    Be careful to not loose information!"),

                            tags$h5("Obs. Buffers are better suited for metric-based CRS projections. The map tab allows crs change (e.g. from WGS84:4326 'degree' unit 
                                to WGS84:3857 Pseudo-Mercator 'metre'."),
                            
                            checkboxInput(inputId = "use_buffer_map", label = "Use an internal buffer to avoid marginal overlaps? (metric-based CRS only)", value = FALSE),

                            conditionalPanel( condition = "input.use_buffer_map == true",  
                            
                                numericInput(inputId = "shrink_factor_buff", label = "Choose an internal buffer", value = 0.01),
                                
                                checkboxGroupInput(inputId = "quantiles_to_buffer", label = "Input quartile areas to buffer",
                                                   
                                           choices = c(1,2,3,4), selected = c(), inline = TRUE),
                            ),
                            
                            conditionalPanel( condition = "input.calculate_Cs > 0",   tags$h6("Do NOT change TAB while running ! This may take a long time depending on the number of species ... ")   ),         
                        ), 
                        
                        box( width = 6,
                             
                            tags$h4("Already have a Cs.csv table with 'sp1', 'sp2' and 'Cs' columns? Check the box and upload."),
                            
                            checkboxInput("Cs_upload_csv", "Upload a Cs csv file instead"),
                            
                            conditionalPanel( condition = "input.Cs_upload_csv == true ",
                                            
                                   fileInput( inputId = "Cs_table", label = "Select a Cs table .csv file", accept = c("text/csv","text/comma-separated-values,text/plain",".csv" ))
                            ), 
                        
                            # tags$h4("Choose a lower limit to spatial congruence Cs ?");# checkboxInput("apply_filter_Cs", "Apply filter to Cs?", value = FALSE),# conditionalPanel( condition = "input.apply_filter_Cs == true", )
                        ),
                    
                        box(width = 6,
                            
                            tags$h4("You can store your Cs table downloading it from here."),
                            
                            downloadButton("download_Cs", "Download Cs.csv"),
                            
                            tags$h5("Cs table summary and head & tail"),
                            
                            textOutput("check_Cs_tables"),
                            
                            fluidRow(
                                
                                column( width = 6,  box(  tableOutput("Cs_head") ) ),
                                
                                column( width = 6,  box(  tableOutput("Cs_tail") ) )
                            )
                        ), # upload Cs
                        
                        box( width = 6,
                        
                            HTML("<code> Check graph's nodes and edges here!</code>"),    
                               
                               fluidRow(    
                                   
                                   column( width = 6, tableOutput("graph_nodes") ),
                                   
                                   column( width = 6, tableOutput("graph_edges") )
                            )
                        ) # check Cs
                    )
                    
                ),
                
                # scan ----
                tabItem("scan",
                        
                    fluidPage(
                    
                        fluidRow(
                            
                            tags$h2("Spatial Congruence ANalysis"),
                            
                            box( width = 12,
                                 
                                tags$h4("After calculating or uploading a Cs table of pairwise spatial similarities, you may parameterize and run SCAN!"),
                                
                                fluidRow(
                                
                                    column( width = 3, numericInput(inputId = "resolution", label = "Resolution (interval between Ct)", value = 0.1) ),
                                     
                                    column( width = 3, numericInput(inputId = "threshold_max", label = "Max value of threshold Ct", value = 0.9, min = 0.2, max = 1) ),
                                     
                                    column( width = 3, numericInput(inputId = "threshold_min", label = "Min value of threshold Ct", value = 0.2, min = 0.05, max = 0.9) ),
                                    
                                    column( width = 3, conditionalPanel(condition = "input.filter_diameter == true",
                                            
                                            numericInput(inputId = "max_diameter", label = "Choose the maximum (network) diameter", value = 15)
                                        )
                                    )
                                ), 
                                #sliderInput( inputId = "threshold_min_max", label = "Select the threshold range to SCAN", value = c(0.8, 1), min = 0.2, max = 1, step = 0.01 ),
                                fluidRow(
                                    
                                    column( width =2, actionButton("run_scan", "SCAN!", class = "btn-warning") ),
                                    
                                    column( width = 3, 
                                            checkboxInput(inputId = "overlap", label = "Overlap criterion: Require overlap between all species?", value = TRUE)),
                                    
                                    column( width = 3,  
                                            checkboxInput(inputId = "filter_diameter", label = "Limit the diameter of the chorotypes networks?", value = TRUE)),
                                     
                                    
                                ),
                                
                                fluidRow(
                                    
                                    column( width = 8, 
                                 
                                            conditionalPanel( condition = "run_scan > 0",  
                                                  
                                                tags$h6("Do NOT change TAB while running ! This also may last a long time... check R console")  )
                                    )
                                )
                            ),
                            
                            box( width = 12,
                                     
                                 tags$h4("Wait for the results and check the parameters of the SCAN analysis"),
                                 
                                 tableOutput("parameters"),
                                 
                                 fluidRow(
                                     
                                     column( width = 8, uiOutput("names_scan_list") ),
                                     
                                     column( width = 4, downloadButton("downloadData", "Download") )
                                 )
                            )
                        ),
                        
                        fluidRow(
                            
                            box(width = 12,
                            
                                tags$h3('SCAN results - download data preview'),
                                
                                dataTableOutput('table_download_preview'),
                                
                                tags$h4('All Chorotypes!'),
                                
                                tableOutput("scan_chorotypes")
                            
                            )
                        )
                    )
                    
                ),
                
                # scan viewer ----
                tabItem("SCAN_viewer",
                        
                    fluidPage(

                        tags$h2("SCAN Viewer"),# tags$h3("Chorotype maps"),

                        box(width = 12,

                            fluidRow(

                                column(width = 3, checkboxInput("graph_from_csv", "Upload a graph (nodes.csv + edges.csv)?", FALSE)),

                                conditionalPanel(   condition = "input.graph_from_csv == true",

                                        column( width = 3,

                                                fileInput(inputId = "graph_nodes", label = "graph NODES in csv file", accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),

                                                fileInput(inputId = "graph_edges", label = "graph EDGES in csv",  accept = c("text/csv","text/comma-separated-values,text/plain",".csv"))
                                        )
                                ),

                                # column(width = 2, numericInput(inputId = "threshold", label = "Threshold (0-1)", value = 0.9, max = 1, min = 0.1)),

                                column(width = 2,
                                       
                                       selectInput(inputId = "palette", label = "choose palette", choices = c("Dark2","Set1","Paired","Accent","Spectral","Greens", "Reds","BrBG" , "RdYlGn","PiYG"))
                                ),
                                
                                column(width = 2,
                                        
                                       selectInput("layout", "graph layout", choices = c("fr", "kk", "dh", "drl", "mds", "gem"))
                                ),
                                
                                column(width = 2,
                                        
                                       numericInput("map_alpha", "map alpha", max = 1, min = 0.01, value = 0.4)
                                ),
                            ),
                        ),
                        
                        HTML("<br>"),
                        
                        # Ct and chorotypes
                        box( width = 12, title = "Set Ct and choose chorotypes",
                            
                            # Choose Congruence threshold of chorotypes
                            column(width = 12,
                                
                                   sliderInput(inputId = "threshold", 
                                                label = "Select a Congruence Threshold - Ct vaule (0-1)", 
                                                value = 0.6, max = 1, ticks = 0.05,
                                                min = 0, step = 0.01, width = 700
                                )
                            ),
                            
                            # choose components to show
                            uiOutput("original_components"),
                        ),
                        
                        # tags$header("Map & Graph"),
                        box(width =12, title = "Map & Graph",
                                          
                            column( width = 6, plotOutput("ggplot_map")), #tags$h3("Map"), 
                            
                            column( width = 6, plotOutput("graph_plot2")) #tags$h3("Graph"),
                        ),
            
                        # leaflet and graph plots
                        # tags$header("Leaflet plot"),
                        
                        box(width = 12, title = "Leaflet plot",
                            
                            box(  leafletOutput("map_plot", height = "500px"), width = 12  ),#tags$h3("Map"), 
                            
                            # box(  plotOutput("graph_plot", height = "500px"), width = 4  ) # tags$h3("Graph"),  
                        ),
                        
                        tags$header("Selected chorotypes (network components) and their species"),
                        
                        dataTableOutput("g_sub_table")
                    )
                ),
                
                # about ----
                tabItem("about_SCAN",
                        
                    fluidPage(
                        
                        tags$h3("SCAN"),
                        
                        tags$h4("https://github.com/cassianogatto/SCAN_engine_app"),
                            
                        fluidRow(
                            
                            infoBox(width = 12,
                                title = "Chorotypes", 
                                value = " Chorotypes are unique combinations of species with spatial congruences 'Cs' higher between themselves than to any species of other such groups.
                                In SCAN species groupings are relative to (and determined by) thresholds of congruence Ct.
                                Each chorotype is a 'community' (in network terminology), as represented in the graph: links are Cs values.
                                The map depicts the actual spatial distribution of each component species of a chorotype.
                                Chorotype may 'evolve' as thresholds get lower, grouping more species, until a criterion of spatial overlap is violated.
                                Some groups exist only at higher Ct; others only at low Ct - it depends on the ecology and history of species and environments.
                                see Gatto & Cohn-Haft 2021 - PlosOne https://doi.org/10.1371/journal.pone.0245818" ),
                            imageOutput("photo1"),
                            box(
                                tags$h4("Abstract (from Gatto & Cohn-Haft 2021)"),
                                tags$h6("Species with congruent geographical distributions, potentially caused by common historical and 
                                ecological spatial processes, constitute biogeographical units called chorotypes. Nevertheless, the
                                degree of spatial range congruence characterizing these groups of species is rarely used as an explicit 
                                parameter. Methods conceived for the identification of patterns of shared ranges often suffer from scale 
                                bias associated with the use of grids, or the incapacity to describe the full complexity of patterns, 
                                from core areas of high spatial congruence, to long gradients of range distributions expanding from 
                                these core areas. Here, we propose a simple analytical method, Spatial Congruence Analysis (SCAN), 
                                which identifies chorotypes by mapping direct and indirect spatial relationships among species. 
                                Assessments are made under a referential value of congruence as an explicit numerical parameter. 
                                A one-layered network connects species (vertices) using pairwise spatial congruence estimates (edges). 
                                This network is then analyzed for each species, separately, by an algorithm which searches for spatial 
                                relationships to the reference species. The method was applied to two datasets: a simulated gradient of 
                                ranges and real distributions of birds. The simulated dataset showed that SCAN can describe gradients 
                                of distribution with a high level of detail. The bird dataset showed that only a small portion of 
                                range overlaps is biogeographically meaningful, and that there is a large variation in types of patterns 
                                that can be found with real distributions. Species analyzed separately may converge on similar or 
                                identical groups, may be nested in larger chorotypes, or may even generate overlapped patterns with no 
                                species in common. Chorotypes can vary from simple ones, composed by few highly congruent species, to 
                                complex, with numerous alternative component species and spatial configurations, which offer insights 
                                about possible processes driving these patterns in distinct degrees of spatial congruence. Metrics 
                                such as congruence, depth, richness, and ratio between common and total areas can be used to describe 
                                chorotypes in detail, allowing comparisons between patterns across regions and taxa."),
                            )                    
                        )
                    )
                    
                )
                
            )   
        )
    )
)

infoBox(title, value = NULL, subtitle = NULL,
        icon = shiny::icon("bar-chart"), color = "aqua", width = 12,
        href = NULL, fill = FALSE)
