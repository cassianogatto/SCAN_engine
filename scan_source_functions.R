# scan functions source
    
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
    
