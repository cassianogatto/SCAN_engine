library(shiny)
library(shinydashboard)
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

# source (SCAN functions)
{ # source
    
    # fix igraph's 'to_subgraph' function
    to_subgraph <- function(graph, ..., subset_by = NULL, delete.vertices = TRUE) {
        if (is.null(subset_by)) {
            subset_by <- active(graph)
            message('SCAN network subsetting by ', subset_by, ' Ct = ', threshold)
        }
        ind <- as_tibble(graph, active = subset_by)
        ind <- mutate(ind, .tidygraph_index = seq_len(n()))
        ind <- filter(ind, ...)
        ind <- ind$.tidygraph_index
        subset <- switch(
            subset_by,
            nodes = induced_subgraph(graph, ind),
            edges = subgraph.edges(graph, ind, delete.vertices = delete.vertices)
        )
        list(subgraph = as_tbl_graph(subset))
    }
    
    SCAN_lite =  function(graph = C, max_Ct = max(graph %>% activate(edges) %>% as_tibble %>% .$Cs),
                          min_Ct = 0.76, Ct_resolution = -0.02, max_diameter = 10, mark_overlap = TRUE, 
                          filter_overlap = FALSE, filter_diameter = FALSE, filter_out_spp = c()    ) {
        #  setup
        chorotypes = list()
        g_spp_all = tibble()
        g_summary_all = tibble()
        Ct_resolution = ifelse(Ct_resolution > 0, Ct_resolution * (-1), Ct_resolution)
        if(isTRUE(filter_overlap)) {mark_overlap = TRUE}
        if(isTRUE(mark_overlap)) {graph = graph %>% activate(nodes) %>% mutate(no_overlap = NA)}
        
        # MAIN LOOP --
        for(threshold in seq(max_Ct,min_Ct,Ct_resolution)){
            
            # any species to be filtered out? (1)
            if(length(filter_out_spp) != 0) {
                graph = graph %>% morph(to_subgraph, subset_by = "nodes",  name %in% filter_out_spp,
                                        remove_multiples = TRUE, delete.vertices= TRUE) %>% mutate(filter = 1) %>%
                    unmorph() }
            
            # get the communities (components) using criteria of Cs, overlap and diameter (2)
            graph = partial_components(graph = graph, threshold = threshold, filter_diameter = filter_diameter,
                                       filter_overlap = filter_overlap)
            
            # get statistics by component (3)
            g = map_by_component(graph = graph, threshold = threshold)
            
            # filter by diameter (4)
            g = g %>% activate(nodes) %>%
                mutate("filter" = ifelse(get(paste0("diameter",threshold)) > max_diameter, threshold, NA))
            
            g = g %>% activate(nodes) %>% mutate("betweenness{threshold}" := round(betweenness(g),1))
            
            # species' list
            g_spp = g %>% activate(nodes) %>% as_tibble() %>%
                group_by(name, get(paste0("components",threshold)), get(paste0("diameter",threshold)),
                         get(paste0("order",threshold)), get(paste0("centrality",threshold))) %>%
                summarize(Ct = threshold, betweenness = get(paste0("betweenness",threshold))) %>%
                select(1,Ct, components = 2, diameter = 3, order = 4, centrality = 5, betweenness) %>%
                arrange(name)
            
            # communities
            g_summary = g %>% activate(nodes) %>% as_tibble %>% group_by(get(paste0("components",threshold)),
                                                                         get(paste0("order",threshold))) %>%
                summarize(Ct = threshold, chorotype_spp = paste(name, collapse = ", "), richness_spp = n(),
                          diameter = max(get(paste0("diameter",threshold))),
                          max_centrality = max(get(paste0("centrality",threshold))),
                          max_betweenness = max(get(paste0("betweenness", threshold)))) %>%
                select(component = 1, Ct, chorotype_spp, richness_spp, diameter, max_centrality, max_betweenness)
            
            # check overlap
            if(isTRUE(mark_overlap)){
                
                g_spp = g_spp %>% mutate(no_overlap = NA)
                g_summary = g_summary %>% mutate(no_overlap = NA)
                
                are_connected = data.frame()
                for(comp in g_summary$component){
                    spp = g_summary %>% filter(component == comp) %>% pull(.,"chorotype_spp") %>% strsplit(.,", ") %>% .[[1]]
                    
                    for (sp1 in spp){
                        for(sp2 in spp[which(spp != sp1)]){
                            conn = tibble(species1 = sp1, species2 = sp2,
                                          connected = igraph::are.connected(graph, sp1,sp2))
                            are_connected = rbind(are_connected, conn)
                        }       }        }
                
                # if all species in a component which sp1 belongs are connected -> TRUE
                connected_nodes_in_components =  are_connected %>%
                    group_by(species1) %>% summarize(all_connected = ifelse(all(connected), TRUE, FALSE)) %>%
                    left_join(g_spp, by = c("species1" = "name")) %>% select(component = 4, name = 1,2) %>%
                    arrange(component, name)
                
                # identify and remove communities in which not all components are connected (overlapped)
                all_connected_components = connected_nodes_in_components %>% group_by(component) %>%
                    summarize(all_connected = ifelse(all(all_connected), TRUE, FALSE))
                
                not_connected_components = all_connected_components %>% filter(all_connected == FALSE) %>% pull(.,"component")
                spp_in_not_connected_components = g_spp %>% filter(components %in% not_connected_components) %>%
                    pull(.,'name')
                
                
                if(length(not_connected_components) > 0 & isTRUE(mark_overlap)) {
                    
                    # write non-all-overlapped chorotypes an species
                    g_spp = g_spp %>% mutate(no_overlap = replace(NA, name %in% spp_in_not_connected_components, threshold ) )
                    
                    g_summary = g_summary %>% mutate(no_overlap = replace(NA, component %in% not_connected_components, threshold))
                    
                    # write those non all-overlapped components to non-overlap column in graph
                    if(isTRUE(mark_overlap)) {
                        graph = graph %>% morph(to_subgraph, subset_by = "nodes",
                                                # criteria to write non-overlap for the first (and only) time in graph
                                                is.na(no_overlap) &                                    # [4]
                                                    name %in% spp_in_not_connected_components,
                                                remove_multiples = TRUE, delete.vertices= TRUE) %>%
                            mutate(no_overlap = threshold) %>%
                            unmorph()
                    }
                }
            }
            
            # update tables
            
            g_spp_all = rbind(g_spp_all, g_spp)
            
            g_summary_all = rbind(g_summary_all, g_summary)
            
        }  # main loop ends
        
        # summarize results and return list of objects
        
        if(isTRUE(mark_overlap)){
            chorotypes[['chorotypes']] = g_summary_all %>% group_by(chorotype_spp, richness_spp, diameter) %>%
                summarise(Ct_max = max(Ct), Ct_min = min(Ct), max_centrality = max(max_centrality),
                          max_betweenness = max(max_betweenness), no_overlap = max(no_overlap)) %>%
                arrange(chorotype_spp, desc(Ct_max))
            
            chorotypes[['all_spp_summary']] = g_spp_all %>% group_by(name, components, order) %>%
                summarise(max_Ct = max(Ct),
                          min_Ct = min(Ct), max_diam = max(diameter), min_diam = min(diameter),
                          max_between = max(betweenness), no_overlap = max(no_overlap))
        } else {
            chorotypes[['chorotypes']] = g_summary_all %>% group_by(chorotype_spp, richness_spp, diameter) %>%
                summarise(Ct_max = max(Ct), Ct_min = min(Ct), max_centrality = max(max_centrality),
                          max_betweenness = max(max_betweenness)) %>%
                arrange(chorotype_spp, desc(Ct_max))
            
            chorotypes[['all_spp_summary']] = g_spp_all %>% group_by(name, components, order) %>% summarise(max_Ct = max(Ct),
                                                                                                            min_Ct = min(Ct), max_diam = max(diameter), min_diam = min(diameter),
                                                                                                            max_between = max(betweenness))
        }
        
        chorotypes[['all_spp']] = g_spp_all
        
        chorotypes[['graph']] = graph
        
        chorotypes[["parameters"]] = tibble(max_diameter = max_diameter, max_Ct = max_Ct, min_Ct = min_Ct,
                                            Ct_resolution = Ct_resolution, mark_overlap = mark_overlap,
                                            filter_overlap = filter_overlap, filter_diameter = filter_diameter)
        
        return(chorotypes)
        
        
    }
    
    partial_components = function (graph = graph, threshold = threshold, filter_diameter = FALSE, filter_depth = FALSE, filter_overlap = FALSE, ...){
        
        print("using 'igraph::group_components' - see more options of community structurig in '?group_components'")
        
        if(isTRUE(filter_overlap) & isTRUE(filter_diameter)) { graph %>% morph(to_subgraph, subset_by = "edges",
                                                                               
                                                                               (Cs >= threshold & is.na(.N()$no_overlap[from]) & is.na(.N()$filter[from])), # check the node respective to the 'from' edge table
                                                                               
                                                                               remove_multiples = TRUE, delete.vertices= TRUE) %>%
                
                activate(edges) %>% mutate("Ct{threshold}" := TRUE) %>%
                
                activate(nodes) %>% mutate("Ct{threshold}" := TRUE) %>%
                
                mutate("components{threshold}" := group_components("weak")) %>% # identify connected elements in COMPONENTS (simpler commuity definition)
                
                unmorph()
            
        } else {
            if(isTRUE(filter_overlap)) {
                graph %>% morph(to_subgraph, subset_by = "edges",
                                (Cs >= threshold & is.na(.N()$no_overlap[from])),
                                remove_multiples = TRUE, delete.vertices= TRUE) %>%
                    activate(edges) %>% mutate("Ct{threshold}" := TRUE) %>%
                    activate(nodes) %>% mutate("Ct{threshold}" := TRUE) %>%
                    mutate("components{threshold}" := group_components("weak")) %>%
                    unmorph()
            } else {
                
                if(isTRUE(filter_diameter)) {
                    graph %>% morph(to_subgraph, subset_by = "edges",
                                    (Cs >= threshold & is.na(.N()$filter[from])),
                                    remove_multiples = TRUE, delete.vertices= TRUE) %>%
                        activate(edges) %>% mutate("Ct{threshold}" := TRUE) %>%
                        activate(nodes) %>% mutate("Ct{threshold}" := TRUE) %>%
                        mutate("components{threshold}" := group_components("weak")) %>% # identify connected elements
                        unmorph()
                } else{
                    graph %>% morph(to_subgraph, subset_by = "edges", Cs >= threshold,
                                    remove_multiples = TRUE, delete.vertices= TRUE) %>%
                        activate(edges) %>% mutate("Ct{threshold}" := TRUE) %>%
                        activate(nodes) %>% mutate("Ct{threshold}" := TRUE) %>%
                        mutate("components{threshold}" := group_components("weak")) %>%
                        unmorph()
                }       }       }       }
    
    # calculate diameter, order, and centrality for each component separately
    map_by_component = function(graph = graph, threshold = threshold){ #, filter_diameter = filter_diameter
        
        # 'I'll suspend this filtering here by now'
        # if(filter_depth) {graph = graph %>% activate(nodes) %>% filter(is.na(depth_filter))}
        
        graph %>% activate(edges) %>% filter(!is.na(get(paste0("Ct",threshold)))) %>%
            activate(nodes) %>% filter(!is.na(get(paste0("Ct",threshold)))) %>%
            # split by components
            morph(to_split, group_by = get(paste0("components",threshold)), subset_by = "nodes") %>%
            # diameter = depth in SCAN
            mutate("diameter{threshold}" := graph_diameter(unconnected = TRUE),
                   # order = richness of species
                   "order{threshold}" := graph_order(),
                   # centrality by component
                   "centrality{threshold}" := centrality_degree()) %>%
            # # betweenness of each node - (cannot calculate betweenness like the above parameters... don't know why...)
            # mutate("betweenness" = betweenness()) %>%
            
            unmorph()
    }
    
    # function overlapping species by quantile
    overlapping_species <- function( map = map(), shrink_factor = 0.01, quartiles_to_buffer = c(3,4)) {
        
        shrink_factor <- ifelse(shrink_factor < 0, shrink_factor, (-1)*shrink_factor)
        
        buffer <- ifelse(map$quantile_area %in% quartiles_to_buffer, sqrt(map$area) * shrink_factor, 0)
        
        map <- st_buffer(map, dist = buffer) #(11s)
        
        overlapping <- st_intersects(map, map, sparse = F) |> as_tibble() |> 
            
            setNames(map$sp) |> mutate(sp1 = map$sp) |> select(sp1, everything()) |> 
            
            pivot_longer(cols = !1, names_to = "sp2") |> filter(value) |> filter(sp1 != sp2) |> 
            
            filter(!duplicated(paste0(pmax(sp1, sp2), pmin(sp1, sp2))))
        
        overlapping
    }
    
    # rowwise mutate vector approach (sp1,sp2, area_sp1, area_sp2, area_overlap)
    overlap_areas_function <- function(map = map, over = overlapping ){
        
        over <- over |> left_join(map, by = c('sp1' = 'sp')) |> select(sp1,sp2, value, geom_sp1 = geometry)
        
        over <- over |> left_join(map, by = c('sp2' = 'sp')) |> select(sp1,sp2, value, geom_sp1, geom_sp2 = geometry)
        
        over <- over |> rowwise() |> mutate(area_overlap = st_area(st_intersection(geom_sp1, geom_sp2)))
        
        over <- over |> rowwise() |> mutate(area_sp1 = st_area(geom_sp1), area_sp2 = st_area(geom_sp2))
        
        over <- over |> select(sp1, sp2, area_sp1, area_sp2, area_overlap)
        
        over 
        
    }
    
}

ui <- dashboardPage( #skin = "green",
        
        dashboardHeader(title = "SCAN engine V_0.21"),
        
        dashboardSidebar( width = '230px',
                          
            # menu items
            sidebarMenu(
              
                menuItem("Species Distribution Maps ", tabName = "maps"),
                
                menuItem("Spatial Congruence", tabName = "Cs_tab"),
                
                menuItem("SCAN Analysis", tabName = "scan"),
                
                menuItem("SCAN Viewer", tabName = "SCAN_viewer"),
                
                menuItem("about SCAN", tabName = "about_SCAN")
            )
        ),
        
        dashboardBody(
            
            # CSS HTML configs ----
            tags$head( 
                
                tags$style( 
                
                    HTML( 
                        
                        #.shiny-input-container { color: #474747;   }
                        ".content {  padding-top: 40px; padding-left: 40px;}
                        .body { background-color: white;   color: black;   }
                        .sidebar {  position:fixed;  overflow:visible; }            
                        .main-header { position: fixed; width = 100%; } "
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
                                
                                tags$h4("Map's column names:"),
                                
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



server <- function(input, output) {

    options(shiny.maxRequestSize=1000*1024^2) # this is required for uploading large datasets
    
    # maps and Cs----
    
    map_upload <- eventReactive(input$get_map,{
        
        shpdf <- input$filemap
        
        if(is.null(shpdf)){    return()    }
        
        previouswd <- getwd()
        
        uploaddirectory <- dirname(shpdf$datapath[1])
        
        setwd(uploaddirectory)
        
        for(i in 1:nrow(shpdf)){   file.rename(shpdf$datapath[i], shpdf$name[i])    }
        
        setwd(previouswd)
        
        map1 <- readOGR(paste(uploaddirectory, shpdf$name[grep(pattern="*.shp$", shpdf$name)], sep="/"))
        
        map1 <- map1 |> st_as_sf(.)
        
        map1
    })
    
    map1 <- reactive({
        
        map_col <- which(map_upload() |> names() == input$colum_sp_map)
        
        if(is.null(map_col)){ return() }
        
        map <- map_upload() |> select(sp = map_col, geometry)
        
        map <- map |> group_by(sp) |> summarise()
        
        return( map )
    })
    
    # the question is if there is an exaggeration in memory requirement after map multiplication (map1(), map2(), map())...
    # check and change projection
    map2 <- reactive({
        
        if (isTRUE(input$modify_crs)) {
            
            if(st_crs(map1())$epsg != as.numeric(input$map_projection)) { map <- map1() |> st_transform(crs = input$map_projection)
            
            } else { map <- map1() }
            
        } else { map <- map1() }
        
        return(map)
    })
    
    # any invalid shapefile?
    invalid_map <-  reactive ({ 
        
        invalid <- which( ! st_is_valid( map2() ) ) 
        
        if( length( invalid ) > 0) {return( invalid )} else { return( c() ) }
    })
    
    # render invalid species output
    output$invalid_map_species <- renderTable( 
        
        map2()[invalid_map(), ] |>  st_drop_geometry() |> select(sp) |> summarise(invalid_species = paste(sp, collapse = ', '))  
    )
    
    # check and fix invalid shapefiles and finally get map()
    map <- reactive({
        
        # Check that the input is valid (could not make it work)
        # validate(  need(is.null(map2()), "Please enter a valid map.") )
        
        if( isTRUE(input$fix_invalid_shapes) & ( length(invalid_map() ) > 0 ) ) {           # map[invalid_map(),] <-  ifelse(length(invalid_map()) > 0, st_make_valid(map[invalid_map(), ]), map) }
            
            map <- map2()
            
            map[invalid_map(),] <- st_make_valid(map[invalid_map(), ]) 
            
        } else { map <- map2() }
        
        return(map)
        
    })
    
    # which pairs of species really overlap in a buffered map?
    overlapping <- reactive({ # overlapping <- overlapping_species(map = map(), shrink_factor = input$shrink_factor_buff, quartiles_to_buffer = input$quantiles_to_buffer) })
        
        # add an area and a quantile_area column ----
        map_areas = st_area(map())
        map <- map() |> mutate(area = map_areas)
        # map <- within(map, quantile_area <- as.integer(cut(map_areas, quantile(map_areas))))
        map <- map() |> mutate(quantile_area = as.integer(cut(map_areas, quantile(map_areas))))
        
        map <- map |> mutate(quantile_area = ifelse(is.na(quantile_area), 1, quantile_area))
        
        # apply buffer
        shrink_factor <- ifelse(input$shrink_factor_buff < 0, input$shrink_factor_buff, (-1) * input$shrink_factor_buff)
        
        buffer <- ifelse( map$quantile_area %in% input$quantiles_to_buffer, sqrt(map_areas) * shrink_factor, 0 )
        
        map <- st_buffer(map, dist = buffer)
        
        # which buffered species-pairs still overlap?
        overlapping <- st_intersects(map, map, sparse = F) |> as_tibble() |> 
            
            setNames(map$sp) |> mutate(sp1 = map$sp) |> select(sp1, everything()) |> 
            
            pivot_longer(cols = !1, names_to = "sp2") |> filter(value) |> filter(sp1 != sp2) |> 
            
            filter(!duplicated(paste0(pmax(sp1, sp2), pmin(sp1, sp2))))
        
        overlapping
    })
    
    Cs_calc <- eventReactive(input$calculate_Cs, {
        
        req( map() )
        
        sf_use_s2(FALSE)
        
        # areas of species and overlaps -> back to the original map (not buffered) (sp1,sp2, area_sp1, area_sp2, area_overlap)
        areas_df <- overlapping() |> left_join( map(), by = c('sp1' = 'sp')) |> select(sp1,sp2, value, geom_sp1 = geometry)
        
        areas_df <- areas_df |> left_join( map(), by = c('sp2' = 'sp')) |> select(sp1,sp2, value, geom_sp1, geom_sp2 = geometry)
        
        areas_df <- areas_df |> rowwise() |> mutate(area_overlap = st_area(st_intersection(geom_sp1, geom_sp2)))
        
        areas_df <- areas_df |> rowwise() |> mutate(area_sp1 = st_area(geom_sp1), area_sp2 = st_area(geom_sp2))
        
        areas_df <- areas_df |> select(sp1, sp2, area_sp1, area_sp2, area_overlap)
        
        areas_df <- areas_df |> units::drop_units()
        
        # apply Cs Index# Hargrove is the default: (area_overlap / area_sp1) * (area_overlap / area_sp2)
        Cs_calc <- areas_df |> mutate( Cs = eval( parse( text = input$cs_similarity_index ) ) ) # instead parse(text = ) may use eval(str2lang(ind)) or str2expression
        
        Cs_calc <- Cs_calc |> filter(Cs >= input$filter_Cs) |>  select(sp1,sp2, Cs) |> arrange(desc(Cs)) |> mutate(Cs = round(Cs, 3))
        
        #if(!is.null(filter_Cs)){  Cs_calc <- Cs_calc |> filter(Cs > input$filter_Cs)}
        
        Cs_calc
        
    })
    
    Cs_up <- reactive({
        
        req(input$Cs_table)
        
        Cs_file <- input$Cs_table
        
        Cs_upload <- read.csv(Cs_file$datapath, header = TRUE) |> as_tibble() |> select(sp1,sp2,Cs)
        
        #if(!is.null(filter_Cs)){  Cs_upload <- Cs_upload |> filter(Cs > input$filter_Cs)}
        
        return(Cs_upload)
    })
    
    Cs <- reactive( {
        
        if(isTRUE(input$Cs_upload_csv)) {  return(  Cs_up()  |> filter(Cs > input$filter_Cs )  )
            
        } else {
            
            return(  Cs_calc() |> filter(Cs > input$filter_Cs )  )  }
        
    })
    
    graph <- reactive({
        
        graph <- Cs() |> as_tbl_graph( from = sp1, to = sp2, directed = FALSE )
        
        graph <- graph %>% igraph::simplify(remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb="first")
        
        graph <- graph %>% as_tbl_graph(directed = FALSE)
        
        return(graph)
    })
    
    # scan
    SCANlist <- eventReactive( input$run_scan, {
        
        req(graph())
        
        thres_max <- input$threshold_max
        
        while(thres_max >= max_thres_graph()) {  thres_max <- thres_max - input$resolution  }
        
        SCANlist <- SCAN_lite(
            graph = graph(),
            max_Ct =  thres_max,
            min_Ct =  input$threshold_min,
            Ct_resolution =  input$resolution,
            max_diameter = input$max_diameter,
            mark_overlap = input$overlap,
            filter_overlap = input$overlap
        )
        
        SCANlist[['graph_nodes']] <- SCANlist[['graph']] |> activate(nodes) |> as_tibble()
        
        SCANlist[['graph_edges']] <- SCANlist[['graph']] |> activate(edges) |> as_tibble()
        
        return(SCANlist)
    })
    
    dataset_SCAN_ouput <- reactive({
        
        switch(input$scan_data_to_download,
               "chorotypes" = SCANlist()[['chorotypes']],
               "all_spp_summary" = SCANlist()[['all_spp_summary']],
               "all_spp" = SCANlist()[['all_spp']],
               "parameters" = SCANlist()[['parameters']],
               "graph_nodes" = SCANlist()[['graph_nodes']],
               "graph_edges" = SCANlist()[['graph_edges']]
        )
    }) # Reactive selected dataset ----
    
    max_thres_graph <- reactive({ graph() |> activate(edges) |> as_tibble() |> summarise(max(Cs)) |> unlist() |> round(2)  })
    
    # viewer v1,1
    threshold <- reactive({   input$threshold   })
    
    g_full <- reactive({
        
        if(!isTRUE(input$graph_from_csv)){
            
            g <- SCANlist()[['graph']]
            
            g <- g  |>  activate(edges) %>% select(from, to, Cs) %>%  filter(Cs >= threshold()) %>%
                activate(nodes) %>% filter(!is.na(get(paste0("components",threshold())))) %>%
                select( name, comps = paste0('components', threshold() ) ) %>% arrange(comps, name)
            # get an object from environment
            # get(input$g) %>% .[['graph']] %>% activate(edges) %>% select(from, to, Cs) %>%  filter(Cs >= threshold()) %>%
            #     activate(nodes) %>% filter(!is.na(get(paste0("components",threshold())))) %>%
            #     select(name, comps = paste0('components',threshold())) %>% arrange(comps, name)
            return(g)
            
        } else {     
            
            node_file <- input$graph_nodes
            nodes <- read.csv(node_file$datapath, header = TRUE)
            
            edge_file <- input$graph_edges
            edges <- read.csv(edge_file$datapath, header = TRUE)
            
            g <- tbl_graph(nodes = nodes, edges = edges, directed = F)
            
            g <- g  |>  activate(edges) %>% select(from, to, Cs) %>%  filter(Cs >= threshold()) %>%
                activate(nodes) %>% filter(!is.na(get(paste0("components",threshold())))) %>%
                select( name, comps = paste0('components', threshold() ) ) %>% arrange(comps, name)
            
            return(g)
        }
    })
    
    original_components <- reactive({
        
        g_full() %>% activate(nodes) %>%  select(comps) %>%
            
            arrange(comps) %>% pull() %>% unique()
    })
    
    g_sub <- reactive({
        
        g_sub <- g_full() %>% activate(nodes) %>%
            
            filter(comps %in% input$selected_components)
        
        g_sub
    })
    
    g_map <- reactive({
        
        g_spp <- g_sub() |> activate(nodes) |> as_tibble()
        
        g_map1 <- right_join( map(), g_spp, by = c('sp' = 'name')) %>% select(comps, everything())
        
        # g_map1 |> st_transform(crs = 4326)
    })
    
    # OUTPUTS
    # maps and Cs
    
    output$map_upload_names <- renderText( paste(map_upload() |> names(), collapse = ', ') )
    
    output$check_Cs_tables <- renderText(
        
        paste("Cs.csv table with ", ncol(Cs()), "columns:" , names(Cs())[1], ",",  names(Cs())[2], ",", names(Cs())[3], ", with", nrow(Cs()), "rows" )
        
    )
    
    output$Cs_head <- renderTable(   Cs() |> head()  )
    
    output$Cs_tail <- renderTable(   Cs() |> tail()  )
    
    output$download_Cs <- downloadHandler(
        
        filename = function() { paste("Cs_table", ".csv", sep = "")  },
        
        content = function(file) {  write.csv(Cs(), file, row.names = FALSE)    }
    )
    
    output$map_species <- renderTable(   map() |> st_drop_geometry() |> select(sp) |> summarise(species = paste(sp,collapse = ', '))  )
    
    output$map_shp_names <- renderText(   paste(map() |> names(), collapse = '-   /   -')  )
    
    output$map_crs <- renderText(  paste( st_crs(map()) ) )
    
    output$map_shp <- renderPlot(
        
        map_shp <- if(nrow(map()) > 250) { 
            
            map()[1:200, "sp" ] |>  plot(col = sf.colors(categorical = TRUE, alpha = 0.5)) 
            
        } else {  
            
            map()[, "sp"] |> plot(col = sf.colors(categorical = TRUE, alpha = 0.5)) 
        }  
    )
    
    output$graph_nodes <- renderTable(  graph() |> activate(nodes) |> as_tibble() |> head() )
    
    output$graph_edges <- renderTable(  graph() |> activate(edges) |> as_tibble() |> head() )
    
    output$test_graph_head <- renderPlot(    (graph() |> activate(edges) |> filter(Cs > 0.5) |> create_layout(layout = "kk") |> ggraph()) + theme_bw() )
    
    output$scan_chorotypes <- renderTable( SCANlist()[['chorotypes']] |> arrange(desc(Ct_max)))
    
    output$parameters <- renderTable({   SCANlist()[['parameters']]  })
    
    # DOWNLOAD
    output$names_scan_list <- renderUI({   
        
        names <- names(SCANlist())
        
        selectInput(inputId = "scan_data_to_download", label = "Choose a SCAN dataset preview (below) or download", choices = names[names !="graph"] )
    }) # Download csv ----# https://shiny.rstudio.com/articles/download.html
    
    output$downloadData <- downloadHandler(
        
        filename = function() { paste(input$scan_data_to_download, ".csv", sep = "")  },
        
        content = function(file) {  write.csv(dataset_SCAN_ouput(), file, row.names = FALSE)    }
    )
    
    output$table_download_preview <- renderDataTable({ 
        
        if(is.null(SCANlist())) { return()  } else {  return( dataset_SCAN_ouput() )   }
        
    })# |> head() 
    
    # viewer
    pal <- reactive({ colorFactor(  palette = input$palette, domain = input$selected_components)  })#original_components()) 
    
    output$original_components <- renderUI({
        
        checkboxGroupInput(inputId = "selected_components", 
                           label = paste("Choose the chorotypes at Ct =", threshold() ) , 
                           choices = original_components(), 
                           inline = TRUE, 
                           selected = NULL) #ifelse(input$select_all_components == TRUE,components, NULL )) #components ) #try ifelse later  ifelse(input$select_all_components == 1, components, NULL)  to select all - but did not work
    })
    
    output$g_sub_table <- renderDataTable({
        
        g_sub() |> activate(nodes) |> as_tibble() |>
            
            group_by(comps) |> summarise(n_spp = n(), species = paste0(name, collapse = ', '))
    })
    
    output$map_plot <- renderLeaflet({  # pal <- colorFactor(input$palette, domain = original_components() )#colorBin(input$palette, domain = original_components() ) #input$selected_components)# colorNumeric # pal <- colorBin(input$palette, domain = original_components(), bins = 7)# labels <- sprintf("%s %s", g_map$comps, g_map$sp) %>% lapply(htmltools::HTML)   # %s use the first 'string'
        
        g_map() |> leaflet() |> 
            
            #addTiles() |>
            
            addProviderTiles("Esri.WorldTopoMap") %>%
            # addEsriTiledMapLayer( url = "https://services.arcgisonline.com/ArcGIS/rest/services/USA_Topo_Maps/MapServer") |> 
            
            addPolygons(  weight = 1,  fillColor = ~ pal()(comps), color = "black", dashArray = "1", fillOpacity = input$map_alpha  ) # highlightOptions = highlightOptions(color = "white", weight = 2, bringToFront = TRUE),# %>% leaflet::addLegend(  pal = pal, values = ~comps,  opacity = 0.7, title = "Chorotypes" )
    })
    
    output$ggplot_map <- renderPlot({
        
        # sa <- #st_as_sfc(scan("C:/Users/cassiano/hubic/Amazon_birds/R_Amazon_birds/SCAN_Viewer/SCAN_engine/SCAN_engine/www/sa.txt"))
        # st_read("www/World_Continents.shp")
        # map <- st_union(g_map(), sa)
        
        gmap_points <- cbind(g_map() |> st_drop_geometry(), st_coordinates(st_centroid(g_map()$geometry)))
        
        ggplot(data = g_map() %>% arrange(comps)) +
            
            geom_sf( aes(fill = comps), # an IFELSE TURNS FILL TO CONTINUOUS... use distiller, otherwise scale_fill_brewer to discrete palette
                     alpha = input$map_alpha, 
                     color = 'black', 
                     show.legend = F) +
            # geom_sf(data = sa, fill = NA, color = 'black') +
            
            scale_fill_distiller( direction = 1, 
                                  palette =   input$palette, 
                                  na.value = "transparent", 
                                  aesthetics = "fill") + #start = 0.2, end = 0.8, #Diverging  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
            # scale_fill_continuous(values = palette(g_map$comps)) +
            
            ggtitle( paste0(" SCAN - Chorotypes at Ct = ", threshold() ) ) + 
            # these functions are from ggspatial (not loaded)
            # annotation_scale(location = "bl", width_hint = 0.5) +
            # annotation_north_arrow(location = "bl", which_north = "true",  pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
            
            geom_text(data= gmap_points, aes(x = X, y = Y, label = sp),
                      color = "black", size = 4, fontface = "italic", check_overlap = TRUE)  +
            
            theme_bw()#theme_classic() #theme_minimal() 
    })
    
    output$graph_plot <- renderPlot({
        
        lay <- create_layout(g_sub(), layout = input$layout)      # lou <- cluster_louvain(g_sub) hmmm there are other options to graph building...
        
        ggraph(lay) +
            
            geom_edge_link( aes( alpha = (Cs+0.75)) , width = 1.25 , show.legend = FALSE) +
            
            geom_node_point( aes( fill = comps),  size = (degree( g_sub(), mode="all") + 15) / 4, shape = 21, show.legend = FALSE) +  # how to synchronize with palette used in MAP? ~pal ??
            
            # scale_fill_distiller( direction = 1, palette = input$palette, na.value = "transparent", aesthetics = "fill") +
            # scale_fill_manual(values = pal(as.factor(comps))) +
            # scale_fill_brewer (values = palette()(comps) ) +
            
            geom_node_text( aes( label = name), size = 3, col = "black", repel= TRUE) +
            
            labs( subtitle = paste0("Ct = ", threshold() )) +
            
            theme_graph()
        
        # tring to implementate hulls showing groups... not succesfull
        # basic_graph2 <-  basic_graph1 + geom_mark_hull(aes(x, y, group = comps), label = comps, label.fontsize = 15, fill = "transparent", lty = "dotted", concavity = 1, expand = unit(3, "mm"), alpha = 0.05) + theme(legend.position = "none")
    })
    
    output$graph_plot2 <- renderPlot({
        
        lay <- create_layout(g_sub(), layout = input$layout)
        
        ggraph(lay) +
            
            geom_edge_link(aes(alpha = (Cs+0.75)) , width = 1.25 , show.legend = FALSE) +
            
            geom_node_point(aes(fill = comps), # how to synchronize with palette used in MAP? ~pal ??
                            
                            size =  (degree(g_sub(), mode="all") + 20) / 4, shape = 21, show.legend = FALSE) +
            
            scale_fill_distiller( direction = 1, palette = input$palette, na.value = "transparent", aesthetics = "fill") +
            
            geom_node_text(aes(label = name), size = 4, col = "black", repel=TRUE) +
            
            labs( subtitle = paste0("Ct = ", threshold() )) +
            
            theme_graph()
        
    })
    
    output$photo <- renderImage({ list( src = file.path("www", "journal.pone.0245818.g004.PNG") ,width = 500, height = 650 #contentType = "image/jpeg",
    )}, deleteFile = FALSE)
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
