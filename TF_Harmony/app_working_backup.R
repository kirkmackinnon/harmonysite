library(data.table)
library(ggplot2)
library(gprofiler2)
library(DT)
library(shiny)
library(plotly)
library(motifStack)
library(ggraph)
library(igraph)
library(tidygraph)
library(visNetwork)
library(viridis)
library(ggdendro)
library(dendextend)
library(universalmotif)
library(shinycssloaders)
library(phylogram)
library(seqinr)
library(ape)

## Spinner style
options(spinner.type = 5, spinner.color = "darkblue")

## Load and preprocess all data
source("data_loading.R", local = TRUE)

## UI definition
source("ui.R", local = TRUE)

## Helper functions (Fisher test, DEG lookup, pairwise matching)
source("helpers.R", local = TRUE)

# Define server logic 
## Backend computation
  server <- function(input, output, session) {
  
  ## Multi-TF DEG lookup. selectize input limited to TF ids
  ## findfile was refactored to use the for loop and cycle through for as many tfs as needed. 
  ## Uses NarPV to calculate overlapping DEGs on the fly.
  ## Previous results are cached so only new TF selections trigger recalculation.
  prev_tarfs <- reactiveVal(list())
  tarfs <- reactive({
    result <- findfile(input$TFs, prev_tarfs())
    prev_tarfs(result)
    result
  })

  prev_matched <- reactiveVal(data.table())
  matched <- reactive({
    validate(need(length(input$TFs) >= 2, "Select at least two TFs for pairwise analysis."))
    result <- matchtidy(tarfs(), prev_matched())
    prev_matched(result)
    result
  })
  
  ## Restricts Harmony to only use the input TFs
  ## Uses newdt - the harmony data with fewer columns 
  subgraph <- reactive({
    validate(need(length(input$TFs) >= 2, "Select at least two TFs for pairwise analysis."))
    newdt[(TF1 %in% input$TFs)][TF2 %in% input$TFs]
    })


  ## This is the motif plot from Pairwise Analysis
  ## Subsets the pfms dataset and then grabs the pre-existing SVG (changed from PNG) file and displays it 
  ## Users can select what type of visual - tree vs stack vs no dendrogram
  ## So it compares motifs based on pfms and then renders images from pre-computed motifs
  
  output$motifplot <- renderImage({
    req(input$pairwise == "pairwisemotifs")
    req(length(input$TFs) >= 1)

    pattern <- paste0("^(", paste0(input$TFs, collapse = "|"), ")_")

    subpfms <- pfms[grepl(pattern, names(pfms))]
    validate(need(length(subpfms) >= 1, "No motifs found for the selected TFs."))

    outfile <- tempfile(fileext = ".svg")
    result <- tryCatch({
      svg(outfile, width = 10, height = 5)
      grid::grid.newpage()
      motifStack(subpfms, layout = input$motiftreestyle)
      dev.off()
      TRUE
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      FALSE
    })
    validate(need(result && file.exists(outfile) && file.info(outfile)$size > 500,
                  "Motif plot failed to render. Please switch tabs and try again."))

    list(src = outfile,
         contentType = 'image/svg+xml',
         width = "100%",
         height = "auto",
         alt = "Motif stack plot")
    
  }, deleteFile = TRUE)
  ## Showing subset harmony table in an interactive DataTable with scroll bars
  ## This is under *pairwise harmony*
  ## Only shows harmony for subsetted TF list based on user input
    output$harmonyTable <- DT::renderDataTable(
      {
        setDT(subgraph())
      },
      selection = 'single',
      extensions = 'Buttons',
      options = list(
        scrollY = '400px',
        paging = FALSE,
        dom = 'Bfrtip',
        buttons = list(list(extend = 'csv', filename = 'pairwise_harmony'))
      )
    )  
    
    ## Reactive function to subset harmony data on Global Analyses based on user input
    ## Most future Global Analysis harmony references use subdt**
    ## Subsetted in such a way that you can specify TF or TF Family
    subdt <- reactive({
      validate(need(length(input$family1) >= 1, "Select at least one TF or TF Family."))
        dt[(((TF1_Family %in% input$family1) &
               (TF2_Family %in% input$family1)
        ) |
          ((TF1 %in% input$family1) &
             (TF2 %in% input$family1))), .(
               TF1,
               TF2,
               Concordant_Intersect,
               Concordant_PValue,
               Concordant_Harmony,
               Discordant_Intersect,
               Discordant_PValue,
               Discordant_Harmony,
               TF1_Family,
               TF2_Family
             )][order(-Concordant_Harmony)]
    })
    
    ## Showing subset harmony table in an interactive DataTable with scroll bars
    ## This is under Global Analyses harmony
    ## Only shows harmony for subsetted TF list based on user input
    output$fullharmonydt <- DT::renderDataTable(
      {
        setDT(subdt())
      },
      selection = 'single',
      extensions = 'Buttons',
      options = list(
        scrollY = '800px',
        paging = FALSE,
        dom = 'Bfrtip',
        buttons = list(list(extend = 'csv', filename = 'global_harmony'))
      )
    )  
    
    ## This is the cell type specificity plot under Pairwise Analyses
    ## Subsets the cell type data to the shared DEGs from the TFs specified by the user
    ## Calculates mean cell type expression and scales it
    ## Then creates a matrix of the data, calculates distance, and then clusters
    ## then generates a ggplot that I wrap in ggplotly for interactivity 
    
    output$ctplot <- renderPlotly({
      subcte <- cte[`Gene ID` %in% matched()$rn]
      subcte[, Mean := mean(TPM), by = .(`Gene ID`, Group)][, Scaled := scale(Mean), by = `Gene ID`]
      subcte <- unique(subcte[, .(Tissue = `short name`, ID = `Gene ID`, Scaled = Scaled, Name = `long name`)])
      
      ## Clustering:
      # convert narrow data to matrix using scaled values
      submat <- as.matrix(dcast.data.table(subcte, ID ~ Tissue, value.var = "Scaled"), rownames = "ID")
      # only use rows without NAs
      submat <- submat[complete.cases(submat),]
      # Calculate euclidean distance from matrix
      d <- dist(submat)
      # use complete clustering on distance matrix
      fit <- hclust(d, method = "complete")
      # assign hierarchy to initial arabidopsis ids based on cluster order
      subcte$ID <- factor(subcte$ID, levels = fit$labels[fit$order], ordered = T)
      
      hm <- ggplot(subcte, aes(x = Tissue, y = ID, fill = Scaled, text = Name)) + 
        geom_raster() + 
        theme_bw(base_size = 20) + 
        theme(
          panel.grid = element_blank(), 
          axis.text.x = element_text(angle = -90, vjust = 0.5),
          axis.text.y = element_blank()
        ) +
        scale_fill_viridis_c() 
      
      ggplotly(hm)
    })
    
    
    ## NxTime overlap plot from Pairwise Analyses - Root data
    ## Subsets NxTime based on shared DEGs
    ## No clustering - instead it orders based Just-In-Time category
    
    output$nxtplot1 <- renderPlotly({
      subjitr <- jitr[GeneID %in% matched()$rn]
      
      subjitr$GeneID <- factor(subjitr$GeneID, levels = subjitr[order(JIT), GeneID])
      
      hm <- ggplot(subjitr, aes(x = JIT, y = GeneID, fill = pvalue)) + 
        geom_raster() + 
        theme_bw(base_size = 20) + 
        theme(
          panel.grid = element_blank(), 
          axis.text.x = element_text(angle = -90, vjust = 0.5),
          axis.text.y = element_blank()
        ) +
        scale_fill_viridis_c() + 
        labs(title = "Roots")

      ggplotly(hm)
    })
    
    ## Same thing as nxtplot1 except its doing it for Shoot data now
    output$nxtplot2 <- renderPlotly({
      subjits <- jits[GeneID %in% matched()$rn]
      
      subjits$GeneID <- factor(subjits$GeneID, levels = subjits[order(JIT), GeneID])
      
      hs <- ggplot(subjits, aes(x = JIT, y = GeneID, fill = pvalue)) +
        geom_raster() +
        theme_bw(base_size = 20) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = -90, vjust = 0.5),
          axis.text.y = element_blank()
        ) +
        scale_fill_viridis_c() +
        labs(title = "Shoots")
      
      ggplotly(hs)
    })
    
    ## GO Term submission
    ## Pulls out shared DEG ids then submits them to GOST server
    ## from g:Profiler package https://biit.cs.ut.ee/gprofiler/gost
    goout <- reactive({
      targs <- split(matched()$rn, f = matched()$Inter)
      gost(targs, organism = "athaliana", multi_query = T)
    }) %>% bindCache(sort(input$TFs))
    
    
    ## Rendering plot output of GO terms
    output$goplot <- renderPlotly({
      
      ggplotly(gostplot(goout(), interactive = F), height = input$goplotheight)
      
    })
    
    ## Rendering GO term table from gost results as a scrollable DataTable
    output$gotable <- DT::renderDT({
      res <- as.data.frame(goout()$result)
      ## Select and order columns to match the original gosttable layout
      ## p_values is a list column from multi_query — take the minimum per term
      res$p_values <- sapply(res$p_values, min)
      DT::datatable(
        res[, c("source", "term_id", "term_name", "term_size", "p_values")],
        extensions = 'Buttons',
        options = list(
          scrollX = TRUE,
          paging = TRUE,
          pageLength = 25,
          autoWidth = FALSE,
          dom = 'Bfrtip',
          buttons = list(list(extend = 'csv', filename = 'GO_terms')),
          columnDefs = list(
            list(width = '30px', targets = c(0, 1)),
            list(width = '22px', targets = 2),
            list(width = '60px', targets = c(4, 5)),
            list(width = '200px', targets = 3),
            list(targets = 5, render = DT::JS(
              "function(data, type, row) {",
              "  if (type === 'display') { return parseFloat(data).toExponential(2); }",
              "  return data;",
              "}"
            ))
          )
        ),
        class = 'compact cell-border'
      ) %>% DT::formatStyle(columns = 1:5, fontSize = '85%')
    })

    ## GO terms are now only computed when the user visits the tab
    ## (previously forced on every load via suspendWhenHidden = FALSE)
    
    
    ## User input - list of TFs to upload
    ## App watches for tfs then reads in the data
    ## Then updates the selectInput "TFs" based on the new data
    ## Pairwise Analyses
    observeEvent(input$tfupload, {
      file <- input$tfupload
      fdt <- fread(file$datapath, header = F)
      
      updateSelectInput(
        session = session,
        inputId = "TFs",
        choices = idoptions,
        selected = fdt[V1 %in% idoptions, V1]
      )
    })
    
    ## Populate allgenes selectize server-side for performance with large gene list
    updateSelectizeInput(session, "allgenes", choices = allgeneids, selected = allgeneids[1], server = TRUE)

    ## User option to upload a list of targets
    ## Target regulation
    ## Updates selectInput to match upload
    observeEvent(input$targetupload, {
      file <- input$targetupload
      fdt <- fread(file$datapath, header = F)

      updateSelectizeInput(
        session = session,
        inputId = "allgenes",
        choices = allgeneids,
        selected = fdt[V1 %in% allgeneids, V1],
        server = TRUE
      )
    })
    
    ## For the harmony cutoff user input
    ## The slider is updated automatically to use the table minimum and maximum as the extreme values
    observeEvent(input$TFs, {
      if (length(input$TFs) < 2) return()
      meltdt <- melt.data.table(subgraph(), id.vars = c("TF1", "TF2"), measure.vars = c("Concordant", "Discordant"))
      meltdt <- meltdt[!(is.na(value) | is.infinite(value))]
      req(nrow(meltdt) > 0)
      updateSliderInput(session, inputId = "HarmonyRange", value = 0, min = 0, max = max(abs(meltdt$value)))
    })
    
    ## This should be the DEG Network plot from Pairwise Analyses
    ## Shows all the overlapping degs
    output$tfnetworkplot <- renderVisNetwork({
      
      ## Getting input from the user
      padjcutoff <- as.numeric(input$tfregpcutoff)
      l2fccutoff <- as.numeric(input$tfreglfccutoff)
      degreecutoff <- as.numeric(input$tfdegreecutoff)
      
      ## Subsetting DEGs to match user cutoffs
      subnar <- narpv[TF %in% input$TFs][padj < padjcutoff][abs(log2FoldChange) > l2fccutoff]
      
      ## Cleaning up the datatable - adding labels based on directionality
      ## Removing self referential rows - no need to show TF over expression here
      ## ordering based activity/color/log2foldchange
      subnar[, group := ifelse(sign(log2FoldChange) == 1, "Upregulated", "Downregulated")]
      subnar <- subnar[!(TF_ID == rn)]
      subnar <- subnar[order(-Activity, Color, log2FoldChange)]

      ## Creating the graph data structure from the DEG subset
      mygraph <- graph_from_data_frame(
        d = subnar[, .(TF_ID, rn, label = TF, value = log2FoldChange, padj, group)], 
        vertices = unique(vertices[(V1 %in% subnar$TF) | (V1 %in% subnar$rn) | (V1 %in% subnar$TF_ID), .(V1, label = TF, group = Activity, value = TedStrength, shape = shape)])
        )
      
      
      ## Calculating network degree
      V(mygraph)$Degree <- igraph::degree(mygraph)
      ## Pruning the network to only show vertices above the user cutoff for degree minimum
      mygraph <- delete.vertices(mygraph, V(mygraph)$Degree < degreecutoff)
      
      ## Setting the colors based on the color scheme established in the setup section - before UI
      V(mygraph)$color <- testcolors[V(mygraph)$group]
      E(mygraph)$color <- testcolors[E(mygraph)$group]

      ## Showing the network graph and defining a handful of options
      visIgraph(mygraph, layout = input$degnetworkstyle) %>% 
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), collapse = T ) %>% 
        visNodes(color = list(border = "black", highlight = "yellow")) %>%
        visGroups(groupname = "Activator", color = "red", shape = "triangle") %>%
        visGroups(groupname = "Minimally Active", color = "gray", shape = "triangle") %>%
        visGroups(groupname = "Repressor", color = "blue", shape = "triangle") %>%
        visGroups(groupname = "Unknown", color = "gray") %>%
        visLegend()
        
    })
    
    ## This is a different network graph showing harmony - TF Regulation tab under Pairwise Analyses
    ## Generates a network graph from subgraph - which is the subsetted harmony datatable
    ## The graph itself shows concordant vs. discordant harmony by colored arrows assuming there is significant harmony between two tfs
    
    output$networkplot <- renderVisNetwork({
      ## Melting datatable, subsetting, adding color labels
      melted <- melt.data.table(subgraph(), id.vars = c("TF1", "TF2"), measure.vars = c("Concordant", "Discordant"), variable.name = "group")
      melted <- melted[!(is.na(value) | is.infinite(value)) & abs(value) > input$HarmonyRange ]
      melted[, color := ifelse(group == "Concordant", "red", "blue")]
      
      
      ## Checks to see if DEGs is checked
      ## Grabs all the TF2s -> TFs that have harmony with TF1 from the melted datatable
      ## Goes through all the TFs in the input and double filters - must be a TF and must be an AtID in tarffs
      ## Then subsets melted to only be the tfs in tf2ids
      if (input$networkdegs){
        
        tf2tfs <- melted$TF2
        tf2ids <- c()
        for (item in tarfs()){
          tf2ids <- c(tf2ids, ids[TF %in% tf2tfs][(Ids %in% item$rn)])
        }
        melted <- melted[TF2 %in% tf2ids$TF]
        
      }
      
      
      ## Creating graph data structure, ensuring it's directional to match harmony 
      mygraph <- as_tbl_graph(as.data.table(melted), directed = T)
      ## calculating network degree
      V(mygraph)$Degree <- igraph::degree(mygraph)
      
      ## Displaying the network graph
      visIgraph(mygraph, layout = input$tfregnetworkstyle) %>%
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), collapse = T, width = "100%") %>%
        visNodes(color = list(border = "black", highlight = "yellow")) %>%
        visGroups(groupname = "Concordant", color = "red") %>%
        visGroups(groupname = "Discordant", color = "blue") %>%
        visLegend() 
      
    })
    
    ## This is the harmony graph on the first page of Pairwise Analyses
    ## creates a ggplot comparing log2FC for shared DEGs between two TFs
    ## points are colored by whether the two TFs agree on up/down regulation of shared TF
    ## plot is wrapped by Plotly for interactivity
    ## Pre-computed correlation exists but isn't displayed - had trouble making it work with Plotly 
    ## Linear modeling is added via geom_smooth
    
    output$harmonyplotly <- renderPlotly({
      setDT(matched())
  
      p <- ggplot(matched(), aes(x = `log2FC.x`, y = `log2FC.y`, color = Harmony, shape = Harmony)) +
        
        geom_point(size = 2, aes(text = rn)) +
        
        theme_bw(base_size = 30) +
        
        theme(
          legend.position = "bottom",
          panel.grid.minor = element_blank()
        ) +
        
        {if(input$harmonyinterfacet)facet_wrap(~Inter, scales = "free", ncol = 1)} + 
        
        scale_color_viridis_d(direction = -1, end = 0.75) +
        
        geom_smooth(aes(group = interaction(Inter, Harmony)), method = "lm", se = F, linewidth = 1) +
        
        geom_vline(xintercept = 0) +
        
        geom_hline(yintercept = 0) + 
        
        geom_line(aes(group = rn), color = "black", alpha = 0.25)
        
       ggplotly(p, height = input$harmonyplotheight)
    })
    
    ## This is the second tab under Global Analyses
    ## Family Harmony
    ## The idea was to show harmony relationships but color coded by family to look for patterns
    ## starts by subsetting the harmony datatable based on user tf family inputs
    ## removes rows with NA harmony and subsets further based on user cutoffs
    ## 
    
    output$relationships <- renderPlotly({
      subdt <- subdt()[, .(TF1, TF2, Concordant_Harmony, Discordant_Harmony, TF1_Family, TF2_Family)]
      subdt <- subdt[!(is.na(Discordant_Harmony)) &
                       (abs(Discordant_Harmony) >= as.numeric(input$familyharmonycutoff))]
      subdt <- subdt[!(is.na(Concordant_Harmony)) &
                       (abs(Concordant_Harmony) >= as.numeric(input$familyharmonycutoff))]

      ## Defining a new column called interaction which is used to color the points
      ## This is a unique identifier for which families are interacting 
      ## this step is used so TF1-MYB <-> TF2-bZIP is treated the same as TF1-bZIP <-> TF2-MYB 
      subdt[, Interaction := paste(pmin(TF1_Family, TF2_Family), pmax(TF1_Family, TF2_Family), sep = ".")]
      
      ## Scatter plot of Harmony values - discordant shown as negatives
      ## colored by family interaction
      ## Text overlay in plotly of the TF ids
      p <- ggplot(subdt, 
                  aes(
                    x = Concordant_Harmony, 
                    y = Discordant_Harmony * (-1), 
                    text = interaction(TF1, TF2), 
                    color = factor(Interaction))) + 
        geom_point(size = 2) + 
        theme(
          legend.position = "none"
        )
      ggplotly(p)
    })
    
    ## Network graph calculating proportionality of one family interacting with another family
    ## Starts with subsetting and user cutoffs
    ## Calculates the mean concordant/discordant harmony of one family with another
    ## Sums the total harmony by family
    ## Then divides the mean by the total to calculate proportion for each family interaction
    ## Then starts creating the graph structure, manually defining edges and ensuring they're directional
    output$familyproportion <- renderVisNetwork({
      subdt <- subdt()[, .(TF1, TF2, Concordant_Harmony, Discordant_Harmony, TF1_Family, TF2_Family)]
      subdt <- subdt[!(is.na(Discordant_Harmony)) &
                       (abs(Discordant_Harmony) >= as.numeric(input$familyharmonycutoff))]
      subdt <- subdt[!(is.na(Concordant_Harmony)) &
                       (abs(Concordant_Harmony) >= as.numeric(input$familyharmonycutoff))]
      subdt[, Interaction := paste(pmin(TF1_Family, TF2_Family),
                                   pmax(TF1_Family, TF2_Family),
                                   sep = ".")]
      
      family_harmony <- subdt[, .(
        mean_H_concordant = mean(Concordant_Harmony, na.rm = TRUE),
        mean_H_discordant = mean(Discordant_Harmony, na.rm = TRUE),
        n_interactions = .N
      ), by = .(TF1_Family, TF2_Family)]
      
      # Normalize concordant harmony
      family_harmony[, total_H_concordant := sum(mean_H_concordant, na.rm = TRUE), by = TF1_Family]
      family_harmony[, prop_H_concordant := mean_H_concordant / total_H_concordant]
      
      
      # Normalize discordant harmony (optional)
      family_harmony[, total_H_discordant := sum(mean_H_discordant, na.rm = TRUE), by = TF1_Family]
      family_harmony[, prop_H_discordant := mean_H_discordant / total_H_discordant]
      
      # -------------------------
      # 1. Define Harmony edges
      # -------------------------
      # Concordant edges: Family1 → Family2
      edges_concordant <- family_harmony[mean_H_concordant > 0, .(
        from = TF1_Family,
        to = TF2_Family,
        value = mean_H_concordant * 10,
        arrows = "to",
        color = "blue",
        title = paste0("Concordant Harmony: ", round(mean_H_concordant, 3))
      )]
      
      # Discordant edges: Family2 → Family1 (reversed direction)
      edges_discordant <- family_harmony[mean_H_discordant > 0, .(
        from = TF2_Family,
        to = TF1_Family,
        value = mean_H_discordant * 10,
        arrows = "to",
        color = "red",
        title = paste0("Discordant Harmony: ", round(mean_H_discordant, 3))
      )]
      
      # Combine both directions
      edges <- rbind(edges_concordant, edges_discordant)
      
  
      # -------------------------
      # 3. Build nodes table with coordinates
      # -------------------------
      # All unique family names
      all_families <- unique(c(edges$from, edges$to))
      
      # Build node table
      nodes <- data.table(
        id = all_families,
        label = all_families
      )
      
      
      # Create igraph object
      g <- graph_from_data_frame(edges[, .(from, to)], directed = TRUE)

      # Sugiyama layout with vertical stretch
      sugiyama_layout <- layout_with_sugiyama(g)$layout

      nodes[, x := sugiyama_layout[, 1] * 70]
      nodes[, y := -sugiyama_layout[, 2] * 1500]

      # -------------------------
      # 4. Plot in visNetwork
      # -------------------------
      visNetwork(nodes, edges) %>%
        visNodes(shape = "dot", size = 25, fixed = FALSE) %>%
        visEdges(smooth = TRUE) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visPhysics(enabled = FALSE)
      
      
    })
    
    ## Reactive objects to calculate harmony clustering for heatmaps and tanglegrams.
    ## Done separately for both X and Y axes as harmony is directional.
    ## Uses harmony_hclust() from helpers.R (was previously ~120 lines of duplicated code).
    conhclustx <- reactive({
      harmony_hclust(subdt(), "Concordant", transpose = FALSE,
        input$familyharmonycutoff, input$familyintersectcutoff, input$familypvcutoff,
        input$distmeth, input$clustmeth)
    })

    conhclusty <- reactive({
      harmony_hclust(subdt(), "Concordant", transpose = TRUE,
        input$familyharmonycutoff, input$familyintersectcutoff, input$familypvcutoff,
        input$distmeth, input$clustmeth)
    })

    dishclustx <- reactive({
      harmony_hclust(subdt(), "Discordant", transpose = FALSE,
        input$familyharmonycutoff, input$familyintersectcutoff, input$familypvcutoff,
        input$distmeth, input$clustmeth)
    })

    dishclusty <- reactive({
      harmony_hclust(subdt(), "Discordant", transpose = TRUE,
        input$familyharmonycutoff, input$familyintersectcutoff, input$familypvcutoff,
        input$distmeth, input$clustmeth)
    })
    
    ## This is the harmony heatmap on the first tab of Global Analyses
    ## Clustering is done is separate reactive objects
    ## Done as a ggplot so that I can make it interactive with Plotly
    ## Harmony is first subsetted and filtered based on user input
    ## Then ordered based on clustered object
    ## Then displayed as a raster graph
    output$conHarmonyHeatmap <- renderPlotly({
      
      subdtcon <- subdt()[, .(TF1, TF2, Concordant_Intersect, Concordant_PValue, Concordant_Harmony, TF1_Family, TF2_Family)]
      
      
      subdtcon <- subdtcon[!(is.na(Concordant_Harmony))][
        (abs(Concordant_Harmony) > as.numeric(input$familyharmonycutoff))][
          Concordant_Intersect > as.numeric(input$familyintersectcutoff)][
            Concordant_PValue < as.numeric(input$familypvcutoff)][
              !(is.infinite(Concordant_Harmony))
            ]
      
      
      conorderx <- conhclustx()$labels[conhclustx()$order]
      conordery <- conhclusty()$labels[conhclusty()$order]

       subdtcon$TF1 <- factor(subdtcon$TF1, levels = conorderx)
       subdtcon$TF2 <- factor(subdtcon$TF2, levels = conordery)
       setkey(subdtcon, TF1)

      hm <- ggplot() +
        
        geom_raster(
          data = subdtcon,
          aes(
            x = TF1,
            y = TF2,
            fill = Concordant_Harmony,

          ),

        ) +

        theme_bw(base_size = 15) + 
        
        theme(
          
          axis.text.x = element_text(
            angle = -90, 
            vjust = 0.5
          ),

          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(), 
          legend.position = "bottom"
        ) + 

        scale_fill_viridis_c(option = "A") +

        ggtitle("Concordant Harmony")
      
      ggplotly(hm) 
      
    })
     
    
    ## Same this as previous graph but for discordant harmony
    ## Heatmap of discordant harmony under Global Analyses
    ## Clustering calculated in separate reactive object
    ## Done as ggplot to allow interactivity with Plotly
    output$disHarmonyHeatmap <- renderPlotly({
      subdtdis <- subdt()[, .(TF1, TF2, Discordant_Intersect, Discordant_PValue, Discordant_Harmony)]
      
      subdtdis <- subdtdis[!(is.na(Discordant_Harmony)) & (abs(Discordant_Harmony) > as.numeric(input$familyharmonycutoff))][
        Discordant_Intersect > as.numeric(input$familyintersectcutoff)][
          Discordant_PValue < as.numeric(input$familypvcutoff)][
            !(is.infinite(Discordant_Harmony))
          ]
      
      
      disorderx <- dishclustx()$labels[dishclustx()$order]
      disordery <- dishclusty()$labels[dishclusty()$order]
      
      subdtdis$TF1 <- factor(subdtdis$TF1, levels = disorderx)
      subdtdis$TF2 <- factor(subdtdis$TF2, levels = disordery)
      setkey(subdtdis, TF1)
      
      hm2 <- ggplot() + 
        
        geom_raster(
          data = subdtdis,
          aes(

            x = TF1, 
            y = TF2,
            fill = Discordant_Harmony * -1,

          ),

        ) +
        

      theme_bw(base_size = 15) + 
        
        theme(
          axis.text.x = element_text(
            angle = -90, 
            vjust = 0.5
            
          ),
          
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          legend.title = element_blank()
        ) + 

        scale_fill_viridis_c(option = "E") +

        labs(title = "Discordant Harmony")
        
      
     ggplotly(hm2)
      })
    
    
  ## Tanglegram with user-selectable left/right dendrograms
  ## Replaces the five individual tanglegram outputs (conmotif, phylocon, phylodis,
  ## dismotif, condis) that were orphaned — not referenced in the UI.
  output$tanglegram <- renderPlot({
    ## Uses pre-computed phylo_dend from data_loading.R instead of recomputing from MSA each render
    choices <- c("Concordant", "Discordant", "Phylogeny", "Up Motif similarity", "Down Motif similarity")

    ph <- phylo_dend
    selids <- unique(subdt()[, TF1])

    subids <- idoptions[!(idoptions %in% selids)]

    muc <- motifupclust()
    muc$labels <- sub("_up", "", muc$labels)
    muc$labels <- as.character(muc$labels)

    mdc <- motifdownclust()
    mdc$labels <- sub("_down", "", mdc$labels)
    mdc$labels <- as.character(mdc$labels)
    
    if (length(subids) > 0) {
      ph <- phylogram::prune(ph, pattern = subids)
    } 
    
    dd <- dendlist(as.dendrogram(conhclustx()), as.dendrogram(dishclustx()), ph, as.dendrogram(muc), as.dendrogram(mdc))
    
    kb <- as.integer(input$kbreaks)
    tanglegram(dd, sub = paste(input$tanglechoice1, "x", input$tanglechoice2), k_branches = kb, k_labels = kb, sort = T, which = c(which(choices == input$tanglechoice1), which(choices == input$tanglechoice2)))
  })
  
  ## This creates the matrix of motif similarity 
  ## Subsets the pwms object based on user input then calls compare_motifs
  motifmat <- reactive({
    subids <- unique(subdt()[, TF1])

    matches <- unlist(lapply(subids, function(tf) {
      which(startsWith(names(pwms), paste0(tf, "_")))
    }))
    
    subpwms <- pwms[unique(matches)]

    compare_motifs(subpwms, method = input$motifcompmethod)
  }) %>% bindCache(sort(input$family1), input$motifcompmethod)
  
  ## Calculates clustering object for motifs
  motifhclust <- reactive({
    motifdist <- dist(motifmat(), method = input$distmeth)
    hclust(motifdist, method = input$clustmeth)
  })
  
  ## Specifically calculating the motif similarity but for motifs calculated on only UP Regulated targets
  ## I don't remember where this is used
  motifupclust <- reactive({
    keep <- grep("_up", rownames(motifmat()), value = TRUE)
    
    submat <- motifmat()[keep, keep]
    motifdist <- dist(submat, method = input$distmeth)
    hclust(motifdist, method = input$clustmeth)
  })

  ## Calculating clustered object but for motifs calculated only on DOWN Regulated targets
  motifdownclust <- reactive({
    keep <- grep("_down", rownames(motifmat()), value = TRUE)
    submat <- motifmat()[keep, keep]

    motifdist <- dist(submat, method = input$distmeth)
    hclust(motifdist, method = input$clustmeth)
  })
  
  
  ## Heatmap of motif similarity from precomputed clustered objects
  
  output$motifheatmap <- renderPlotly({
    
    
    motiforder <- motifhclust()$labels[motifhclust()$order]
      
    motifdt <- as.data.table(motifmat(), keep.rownames = "TF1")
    motifnar <- melt.data.table(motifdt, id.vars = "TF1", variable.name = "TF2")
    
    motifnar$TF1 <- factor(motifnar$TF1, levels = motiforder)
    motifnar$TF2 <- factor(motifnar$TF2, levels = motiforder)

    ggplot(motifnar, aes(x = TF1, y = TF2, fill = value)) + geom_raster()
  })
  
  
  ## Motif dendrogram from Global Analyses
  ## Returns an interactive dendrogram
  output$motifdend <- renderPlotly({
    dd <- as.dendrogram(motifhclust())
    ddplot <- ggdendrogram(dd)
    ggplotly(ddplot)
  })
  
  ## This is a datatable over in Target Regulation - third tab
  ## Returns a datatable subsetted and in narrow format from the DEGs data
  ## Shows all TFs that target any of the provided Targets in the allgenes input
  output$tfregdt <- DT::renderDT({
    DT::datatable(
      setDT(targetsubnar()),
      extensions = 'Buttons',
      options = list(
        dom = 'Bfrtip',
        buttons = list(list(extend = 'csv', filename = 'target_DEGs'))
      )
    )
  })

  ## This is a datatable in Target Regulation - third tab
  ## Takes the target subsetted/filtered data from DEGs data
  ## and then filters the vertices data structure based on DEG target regulation
  ## Returns the datatable showing TF, Family, Transcriptional Enhancer Strength, and shape - is a tf or not
  output$targetregdt <- DT::renderDT({
    DT::datatable(
      setDT(subvs()),
      extensions = 'Buttons',
      options = list(
        dom = 'Bfrtip',
        buttons = list(list(extend = 'csv', filename = 'target_regulators'))
      )
    )
  })
  
  ## reactive function to create targetsubnar datatable
  ## uses User supplied cutoffs to filter and subset DEG data
  ## Restricts based on a list of Targets
  targetsubnar <- reactive({
    validate(need(length(input$allgenes) >= 1, "Select at least one target gene."))
    padjcutoff <- as.numeric(input$regpcutoff)
    l2fccutoff <- as.numeric(input$reglfccutoff)
    
    subnar <- narpv[rn %in% input$allgenes][padj < padjcutoff][abs(log2FoldChange) > l2fccutoff]
  })
  
  ## Reactive function that subsets and filters vertice datastructure
  ## Restricts vertices to only the TFs regulating a given set of targets
  ## Relabels resulting output columns for readability
  subvs <- reactive({
    subnar <- targetsubnar()
    unique(vertices[(V1 %in% subnar$TF) | (V1 %in% subnar$rn) | (V1 %in% subnar$TF_ID), .(V1, label = TF, group = Family, value = TedStrength, shape = shape)])
  })
  
  
  ## Generating the actual network graph for the Target Regulation tab
  ## Takes the data from DEG data and vertice data (targetsubnar() and subvs() just created)
  ## Creates a directed graph then makes it interactive with VisNetwork
  output$regulators <- renderVisNetwork({
    
    subnar <- targetsubnar()

    mygraph <- graph_from_data_frame(
      
      subnar[,.(TF_ID, rn, padj, value = log2FoldChange, color = Color, label = TF)], directed = T, 
      vertices = subvs()
      
    )
    
    
    V(mygraph)$Degree <- igraph::degree(mygraph)
    visIgraph(
      mygraph, 
      layout = input$targetnetworkstyle
      
              ) %>% visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), collapse = T ) %>%
              visLegend()
    
  })

}

# Run the application 
shinyApp(ui = ui, server = server)