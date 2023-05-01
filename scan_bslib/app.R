

library(shiny)
source("scan_source_functions.R") 

ui <-
    navbarPage("SCAN_engine", id="scan_nav",
               
                tabPanel("Interactive map",
                        
                         # check this class
                        div(class="outer",
                            
                            tags$head(
                                # Include our custom CSS
                                includeCSS("styles.css"),
                                
                                #includeScript("gomap.js")
                            ),
                            
                            # If not using custom CSS, set height of leafletOutput to a number instead of percent
                            leafletOutput("map_shp"),
                            
                            absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                                          draggable = TRUE, top = 60, left = 20, right = "auto", bottom = "auto",
                                          width = 530, height = "auto",
                                          
                                          h2("Load Species' Maps"),
                                          
                                          # input map
                                          fileInput( inputId = "filemap",  label = "Choose shape-files (.shp + .shx + .dbl + .prj).", 
                                                     accept = c('.shp','.dbf','.sbn','.sbx','.shx',".prj"),width = '500px', multiple=TRUE),  
                                          
                                          # sp column
                                          inputPanel( textInput("colum_sp_map","Which one is the species ID?", value = "sp")),
                                          
                                          # GET MAP!!!     
                                          actionButton("get_map", "Load map!",  class = "btn-warning", color = 'black'),
                                          
                              ),
                            
                            
                            
                            

# Define server logic required to draw a histogram
server <- function(input, output) {

    
}

# Run the application 
shinyApp(ui = ui, server = server)
