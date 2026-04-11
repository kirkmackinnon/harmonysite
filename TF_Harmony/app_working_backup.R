library(data.table)
library(ggplot2)
library(gprofiler2)
library(DT)
library(shiny)
library(plotly)
library(motifStack)
library(igraph)
library(visNetwork)
library(viridis)
library(ggdendro)
library(dendextend)
library(universalmotif)
library(shinycssloaders)
library(phylogram)
library(seqinr)
library(ape)
library(ggrepel)

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
      svg(outfile, width = 10, height = 7)
      grid::grid.newpage()
      motifStack(subpfms, layout = "tree")
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

  ## Pairwise motif similarity table
  ## Uses compare_motifs (ALLR_LL) on the selected TFs' PWMs
  ## Displays as a long-format table of TF1 vs TF2 similarity scores
  output$motifSimTable <- DT::renderDT({
    req(input$pairwise == "pairwisemotifs")
    req(length(input$TFs) >= 2)

    pattern <- paste0("^(", paste0(input$TFs, collapse = "|"), ")_")
    matches <- which(grepl(pattern, names(pwms)))
    req(length(matches) >= 2)

    subpwms <- pwms[matches]
    sim_mat <- compare_motifs(subpwms, method = "ALLR_LL")

    sim_dt <- as.data.table(sim_mat, keep.rownames = "Motif1")
    sim_long <- melt.data.table(sim_dt, id.vars = "Motif1", variable.name = "Motif2", value.name = "Similarity")
    sim_long <- sim_long[as.character(Motif1) < as.character(Motif2)][order(-Similarity)]

    DT::datatable(sim_long,
      options = list(scrollY = '600px', paging = FALSE, dom = 't'),
      rownames = FALSE
    ) %>% DT::formatRound("Similarity", 3)
  })

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
      
      ## Clustering rows (genes) and columns (tissues):
      # convert narrow data to matrix using scaled values
      submat <- as.matrix(dcast.data.table(subcte, ID ~ Tissue, value.var = "Scaled"), rownames = "ID")
      # only use rows without NAs
      submat <- submat[complete.cases(submat),]
      # Cluster rows (genes)
      row_fit <- hclust(dist(submat), method = "complete")
      subcte$ID <- factor(subcte$ID, levels = row_fit$labels[row_fit$order], ordered = T)
      # Cluster columns (tissues)
      col_fit <- hclust(dist(t(submat)), method = "complete")
      subcte$Tissue <- factor(subcte$Tissue, levels = col_fit$labels[col_fit$order], ordered = T)
      
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
      
      subjitr$GeneID <- factor(subjitr$GeneID, levels = subjitr[order(pvalue), GeneID])
      
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
      
      subjits$GeneID <- factor(subjits$GeneID, levels = subjits[order(pvalue), GeneID])
      
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
      
      ggplotly(gostplot(goout(), interactive = F), height = input$harmonyplotheight)
      
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
    
    ## DEG Networks - shows TF->target edges for selected TFs
    ## TFs as module-colored triangles with synonyms, targets as boxes
    ## Louvain module detection on harmony between the selected TFs
    output$tfnetworkplot <- renderVisNetwork({
      req(length(input$TFs) >= 1)

      padjcutoff <- as.numeric(input$tfregpcutoff)
      l2fccutoff <- as.numeric(input$tfreglfccutoff)
      degreecutoff <- as.numeric(input$tfdegreecutoff)

      subnar <- narpv[TF %in% input$TFs][padj < padjcutoff][abs(log2FoldChange) > l2fccutoff]
      subnar <- subnar[!(TF_ID == rn)]
      req(nrow(subnar) > 0)

      ## Module detection on selected TFs via harmony
      hit_tf_ids <- unique(subnar$TF_ID)
      tf_lk <- unique(narpv[, .(TF_ID, TF)])

      if (length(hit_tf_ids) >= 2) {
        h <- copy(dt)
        h <- merge(h, tf_lk, by.x = "TF1", by.y = "TF", all.x = TRUE)
        setnames(h, "TF_ID", "TF1_ID")
        h <- merge(h, tf_lk, by.x = "TF2", by.y = "TF", all.x = TRUE)
        setnames(h, "TF_ID", "TF2_ID")
        h[, abs_harmony := abs(Concordant_Harmony - Discordant_Harmony)]

        h_sub <- h[TF1_ID %in% hit_tf_ids & TF2_ID %in% hit_tf_ids & TF1_ID != TF2_ID,
                   .(TF1_ID, TF2_ID, abs_harmony)]
        mods <- compute_louvain_modules(h_sub, weight_col = "abs_harmony", cutoff = 0.02)
        if (is.null(mods)) mods <- data.table(TF_ID = hit_tf_ids, module = seq_along(hit_tf_ids))
      } else {
        mods <- data.table(TF_ID = hit_tf_ids, module = 1L)
      }

      module_colors <- c(
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
        "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62"
      )

      tf_info <- unique(subnar[, .(TF_ID, TF)])
      tf_info <- merge(tf_info, mods, by = "TF_ID", all.x = TRUE)
      tf_info[is.na(module), module := 0L]

      tf_nodes <- tf_info[, .(
        id = TF_ID,
        label = TF,
        group = "TF",
        shape = "triangle",
        value = 10,
        color = module_colors[((module - 1) %% length(module_colors)) + 1],
        title = paste0("<b>", TF, "</b> (", TF_ID, ")<br>Module: ", module)
      )]

      ## Degree filter on targets
      target_deg <- subnar[, .(degree = uniqueN(TF_ID)), by = rn]
      keep_targets <- target_deg[degree >= degreecutoff, rn]
      subnar <- subnar[rn %in% keep_targets]
      req(nrow(subnar) > 0)

      gene_ids <- setdiff(unique(subnar$rn), tf_nodes$id)
      gene_nodes <- data.table(
        id = gene_ids, label = gene_ids,
        group = "Gene", shape = "box", value = 1, color = "#d9d9d9", title = gene_ids
      )

      nodes <- rbind(tf_nodes, gene_nodes, fill = TRUE)

      edges <- subnar[, .(
        from = TF_ID, to = rn, arrows = "to",
        width = pmax(1, 2 * abs(log2FoldChange)),
        title = paste0(TF, " -> ", rn, "<br>log2FC: ", round(log2FoldChange, 3), "<br>padj: ", signif(padj, 3)),
        color = ifelse(log2FoldChange > 0, "darkred", "darkblue")
      )]

      visNetwork(nodes, edges) %>%
        visGroups(groupname = "TF", shape = "triangle") %>%
        visGroups(groupname = "Gene", shape = "box") %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                   nodesIdSelection = TRUE, collapse = TRUE) %>%
        visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
        visPhysics(solver = "forceAtlas2Based",
                   stabilization = list(enabled = TRUE, iterations = 800, fit = TRUE)) %>%
        visEvents(stabilizationIterationsDone = "function () { this.setOptions({physics: false}); }")
    })

    ## TF Regulation - shows harmony edges between selected TFs
    ## TFs as module-colored triangles with synonyms
    ## Concordant edges in red, discordant in blue, filtered by HarmonyRange slider
    output$networkplot <- renderVisNetwork({
      req(length(input$TFs) >= 2)

      sg <- subgraph()
      melted <- melt.data.table(sg, id.vars = c("TF1", "TF2"),
                                measure.vars = c("Concordant", "Discordant"), variable.name = "group")
      melted <- melted[!(is.na(value) | is.infinite(value)) & abs(value) > input$HarmonyRange]
      req(nrow(melted) > 0)

      ## Module detection on selected TFs
      tf_lk <- unique(narpv[, .(TF_ID, TF)])
      all_tfs <- unique(c(melted$TF1, melted$TF2))
      tf_info <- tf_lk[TF %in% all_tfs]

      if (nrow(tf_info) >= 2) {
        h <- copy(dt)
        h <- merge(h, tf_lk, by.x = "TF1", by.y = "TF", all.x = TRUE)
        setnames(h, "TF_ID", "TF1_ID")
        h <- merge(h, tf_lk, by.x = "TF2", by.y = "TF", all.x = TRUE)
        setnames(h, "TF_ID", "TF2_ID")
        h[, abs_harmony := abs(Concordant_Harmony - Discordant_Harmony)]

        h_sub <- h[TF1_ID %in% tf_info$TF_ID & TF2_ID %in% tf_info$TF_ID & TF1_ID != TF2_ID,
                   .(TF1_ID, TF2_ID, abs_harmony)]
        mods <- compute_louvain_modules(h_sub, weight_col = "abs_harmony", cutoff = 0.02)
        if (is.null(mods)) mods <- data.table(TF_ID = tf_info$TF_ID, module = seq_along(tf_info$TF_ID))
      } else {
        mods <- data.table(TF_ID = tf_info$TF_ID, module = 1L)
      }

      module_colors <- c(
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
        "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62"
      )

      tf_info <- merge(tf_info, mods, by = "TF_ID", all.x = TRUE)
      tf_info[is.na(module), module := 0L]

      nodes <- tf_info[, .(
        id = TF,
        label = TF,
        shape = "triangle",
        value = 10,
        color = module_colors[((module - 1) %% length(module_colors)) + 1],
        title = paste0("<b>", TF, "</b> (", TF_ID, ")<br>Module: ", module)
      )]

      edges <- melted[, .(
        from = TF1, to = TF2, arrows = "to",
        width = pmax(0.5, abs(value)),
        title = paste0(TF1, " -> ", TF2, "<br>", group, ": ", round(value, 3)),
        color = ifelse(group == "Concordant", "red", "blue"),
        dashes = ifelse(group == "Discordant", TRUE, FALSE)
      )]

      visNetwork(nodes, edges) %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                   nodesIdSelection = TRUE, collapse = TRUE) %>%
        visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
        visPhysics(solver = "forceAtlas2Based",
                   stabilization = list(enabled = TRUE, iterations = 800, fit = TRUE)) %>%
        visEvents(stabilizationIterationsDone = "function () { this.setOptions({physics: false}); }")
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

    ## Color axis labels red (up) / blue (down)
    label_colors_x <- ifelse(grepl("_up$", levels(motifnar$TF1)), "red", "blue")
    label_colors_y <- ifelse(grepl("_up$", levels(motifnar$TF2)), "red", "blue")

    hm <- ggplot(motifnar, aes(x = TF1, y = TF2, fill = value)) +
      geom_raster() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = -90, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7)
      )

    p <- ggplotly(hm)

    ## Plotly doesn't support per-tick colors natively.
    ## Workaround: replace tick labels with colored HTML spans.
    colored_x <- mapply(function(lab, col) {
      paste0("<span style='color:", col, "'>", lab, "</span>")
    }, levels(motifnar$TF1), label_colors_x, USE.NAMES = FALSE)

    colored_y <- mapply(function(lab, col) {
      paste0("<span style='color:", col, "'>", lab, "</span>")
    }, levels(motifnar$TF2), label_colors_y, USE.NAMES = FALSE)

    p <- layout(p,
      xaxis = list(tickvals = seq_along(levels(motifnar$TF1)) - 1,
                   ticktext = colored_x),
      yaxis = list(tickvals = seq_along(levels(motifnar$TF2)) - 1,
                   ticktext = colored_y)
    )
    p
  })


  ## Motif dendrogram from Global Analyses
  ## Labels colored red (up) / blue (down), branches colored by up/down majority
  output$motifdend <- renderPlot({
    dd <- as.dendrogram(motifhclust())
    labels <- labels(dd)
    label_cols <- ifelse(grepl("_up$", labels), "red", "blue")
    dd <- dendextend::set(dd, "labels_col", label_cols)
    dd <- dendextend::set(dd, "labels_cex", 0.7)
    plot(dd, main = "Motif Similarity Dendrogram")
  })
  
  ## reactive function to create targetsubnar datatable
  ## uses User supplied cutoffs to filter and subset DEG data
  ## Restricts based on a list of Targets
  targetsubnar <- reactive({
    validate(need(length(input$allgenes) >= 1, "Select at least one target gene."))
    padjcutoff <- as.numeric(input$regpcutoff)
    l2fccutoff <- as.numeric(input$reglfccutoff)

    narpv[rn %in% input$allgenes][padj < padjcutoff][abs(log2FoldChange) > l2fccutoff]
  })

  ## TF ID to TF name lookup from narpv
  tf_lookup <- reactive({
    unique(narpv[, .(TF_ID, TF)])
  })

  ## Harmony table enriched with TF IDs and derived weight columns for module detection
  target_harmony_dt <- reactive({
    lk <- tf_lookup()
    h <- copy(dt)

    h <- merge(h, lk, by.x = "TF1", by.y = "TF", all.x = TRUE)
    setnames(h, "TF_ID", "TF1_ID")
    h <- merge(h, lk, by.x = "TF2", by.y = "TF", all.x = TRUE)
    setnames(h, "TF_ID", "TF2_ID")

    h[, Harmony := Concordant_Harmony - Discordant_Harmony]
    h[, abs_harmony := abs(Harmony)]
    h[, abs_harmony_sq := abs(Harmony)^2]
    h[, concordant_weight := fifelse(is.na(Concordant_Harmony), 0, Concordant_Harmony)]
    h[, discordant_weight := fifelse(is.na(Discordant_Harmony), 0, Discordant_Harmony)]

    h[]
  })

  ## Per-TF summary for the selected gene set
  tf_reg_summary <- reactive({
    subnar <- targetsubnar()
    req(nrow(subnar) > 0)

    out <- subnar[, .(
      TF = first(TF),
      Family = first(Family),
      n_targets = uniqueN(rn),
      n_up_targets = uniqueN(rn[log2FoldChange > 0]),
      n_down_targets = uniqueN(rn[log2FoldChange < 0]),
      mean_abs_l2fc = mean(abs(log2FoldChange), na.rm = TRUE),
      mean_l2fc = mean(log2FoldChange, na.rm = TRUE),
      mean_up_l2fc = mean(log2FoldChange[log2FoldChange > 0], na.rm = TRUE),
      mean_down_l2fc = mean(log2FoldChange[log2FoldChange < 0], na.rm = TRUE),
      sum_abs_l2fc = sum(abs(log2FoldChange), na.rm = TRUE),
      sum_up_l2fc = sum(log2FoldChange[log2FoldChange > 0], na.rm = TRUE),
      sum_down_l2fc = sum(log2FoldChange[log2FoldChange < 0], na.rm = TRUE),
      best_padj = min(padj, na.rm = TRUE)
    ), by = TF_ID]

    for (j in c("mean_up_l2fc", "mean_down_l2fc")) {
      out[is.nan(get(j)), (j) := 0]
    }
    out[]
  })

  ## Mean Harmony of each hit TF to the other hit TFs
  tf_mean_harmony <- reactive({
    reg <- tf_reg_summary()
    hit <- unique(reg$TF_ID)
    req(length(hit) >= 1)

    hdt <- target_harmony_dt()
    h <- hdt[
      TF1_ID %in% hit & TF2_ID %in% hit & TF1_ID != TF2_ID,
      .(TF1_ID, TF2_ID, Concordant_Harmony, Discordant_Harmony, Harmony)
    ]

    out_sum <- h[, .(
      sum_concordant_harmony_out = sum(Concordant_Harmony[!is.na(Concordant_Harmony)], na.rm = TRUE),
      sum_discordant_harmony_out = sum(Discordant_Harmony[!is.na(Discordant_Harmony)], na.rm = TRUE),
      mean_concordant_harmony_out = mean(Concordant_Harmony[!is.na(Concordant_Harmony) & Concordant_Harmony != 0], na.rm = TRUE),
      mean_discordant_harmony_out = mean(Discordant_Harmony[!is.na(Discordant_Harmony) & Discordant_Harmony != 0], na.rm = TRUE),
      n_harmony_connections_out = uniqueN(TF2_ID[!is.na(Harmony) & Harmony != 0])
    ), by = .(TF_ID = TF1_ID)]

    in_sum <- h[, .(
      sum_concordant_harmony_in = sum(Concordant_Harmony[!is.na(Concordant_Harmony)], na.rm = TRUE),
      sum_discordant_harmony_in = sum(Discordant_Harmony[!is.na(Discordant_Harmony)], na.rm = TRUE),
      mean_concordant_harmony_in = mean(Concordant_Harmony[!is.na(Concordant_Harmony) & Concordant_Harmony != 0], na.rm = TRUE),
      mean_discordant_harmony_in = mean(Discordant_Harmony[!is.na(Discordant_Harmony) & Discordant_Harmony != 0], na.rm = TRUE),
      n_harmony_connections_in = uniqueN(TF1_ID[!is.na(Harmony) & Harmony != 0])
    ), by = .(TF_ID = TF2_ID)]

    out <- merge(out_sum, in_sum, by = "TF_ID", all = TRUE)
    num_cols <- setdiff(names(out), "TF_ID")
    for (j in num_cols) {
      out[is.na(get(j)) | is.nan(get(j)), (j) := 0]
    }
    out[]
  })

  ## Louvain module detection on the hit TFs
  tf_modules <- reactive({
    reg <- tf_reg_summary()
    hit <- unique(reg$TF_ID)
    req(length(hit) >= 1)

    hdt <- target_harmony_dt()
    weight_col <- switch(
      input$louvain_weight_mode,
      "abs_harmony" = "abs_harmony",
      "abs_harmony_sq" = "abs_harmony_sq",
      "concordant" = "concordant_weight",
      "discordant" = "discordant_weight",
      "abs_harmony"
    )

    h <- hdt[
      TF1_ID %in% hit & TF2_ID %in% hit & TF1_ID != TF2_ID,
      .(TF1_ID, TF2_ID, abs_harmony, abs_harmony_sq, concordant_weight, discordant_weight)
    ]

    mods <- compute_louvain_modules(h, weight_col = weight_col, cutoff = 0.02)

    if (is.null(mods)) {
      return(data.table(TF_ID = hit, module = seq_along(hit)))
    }
    mods[]
  })

  ## Combined per-TF data for scatter plot and network
  tf_scatter_dt <- reactive({
    reg  <- tf_reg_summary()
    mh   <- tf_mean_harmony()
    mods <- tf_modules()

    out <- merge(reg, mh, by = "TF_ID", all.x = TRUE)
    out <- merge(out, mods, by = "TF_ID", all.x = TRUE)

    fill_zero_cols <- c(
      "sum_concordant_harmony_out", "sum_discordant_harmony_out",
      "mean_concordant_harmony_out", "mean_discordant_harmony_out",
      "n_harmony_connections_out",
      "sum_concordant_harmony_in", "sum_discordant_harmony_in",
      "mean_concordant_harmony_in", "mean_discordant_harmony_in",
      "n_harmony_connections_in"
    )
    for (j in fill_zero_cols) {
      if (j %in% names(out)) {
        out[is.na(get(j)) | is.nan(get(j)), (j) := 0]
      }
    }

    out[, sum_concordant_harmony := sum_concordant_harmony_out + sum_concordant_harmony_in]
    out[, sum_discordant_harmony := sum_discordant_harmony_out + sum_discordant_harmony_in]
    out[, mean_concordant_harmony := mean_concordant_harmony_out + mean_concordant_harmony_in]
    out[, mean_discordant_harmony := mean_discordant_harmony_out + mean_discordant_harmony_in]
    out[, n_harmony_connections := n_harmony_connections_out + n_harmony_connections_in]
    out[, sum_harmony := sum_concordant_harmony - sum_discordant_harmony]
    out[, mean_harmony := mean_concordant_harmony - mean_discordant_harmony]
    out[is.na(module), module := 0L]
    out[, label := fifelse(is.na(TF) | TF == "", TF_ID, TF)]

    mod_sizes <- out[, .(module_size = .N), by = module]
    out <- merge(out, mod_sizes, by = "module", all.x = TRUE)
    out[]
  })

  ## Update module selector when scatter data changes
  observe({
    d <- tf_scatter_dt()
    req(nrow(d) > 0)
    mods <- sort(unique(d$module))
    updateSelectInput(session, "module_select",
      choices = c("All", as.character(mods)), selected = "All")
  })

  ## Scatter plot: # targets vs outgoing Harmony balance, colored by module
  output$tfHarmonyScatter <- renderPlot({
    d <- tf_scatter_dt()
    req(nrow(d) > 0)

    d[, harmony_y := sum_concordant_harmony_out - sum_discordant_harmony_out]

    ggplot(d, aes(x = n_targets, y = harmony_y, color = factor(module))) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(size = n_targets), alpha = 0.85) +
      ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
      labs(x = "# regulated target genes", y = "Outgoing Harmony balance",
           color = "Module", size = "# targets") +
      theme_bw()
  })

  ## Heatmap of TF x Target regulation direction, faceted by module
  heatmap_dt <- reactive({
    subnar <- targetsubnar()
    reg    <- tf_scatter_dt()
    d <- merge(subnar, reg[, .(TF_ID, label, module)], by = "TF_ID", all.x = TRUE)
    d[]
  })

  output$targetRegHeatmap <- renderPlot({
    d <- heatmap_dt()
    req(nrow(d) > 0)

    mat_dt <- d[, .(value = mean(log2FoldChange, na.rm = TRUE)), by = .(label, rn, module)]
    mat_wide <- data.table::dcast(mat_dt, label + module ~ rn, value.var = "value", fill = NA)

    row_meta <- mat_wide[, .(label, module)]
    mat <- as.matrix(mat_wide[, !c("label", "module")])
    rownames(mat) <- mat_wide$label

    mat_clust <- mat
    mat_clust[is.na(mat_clust)] <- 0

    if (ncol(mat_clust) > 1) {
      col_cor <- cor(mat_clust, use = "pairwise.complete.obs")
      col_dist <- as.dist(1 - col_cor)
      col_hc <- hclust(col_dist, method = "complete")
      gene_order <- colnames(mat_clust)[col_hc$order]
    } else {
      gene_order <- colnames(mat_clust)
    }

    plot_dt <- data.table::melt(
      data.table(label = rownames(mat), mat, check.names = FALSE),
      id.vars = "label", variable.name = "rn", value.name = "log2FoldChange"
    )
    plot_dt <- merge(plot_dt, row_meta, by = "label", all.x = TRUE)
    plot_dt[, label := factor(label, levels = unique(plot_dt[order(module, label)]$label))]
    plot_dt[, rn := factor(rn, levels = gene_order)]
    plot_dt[, reg_dir := sign(log2FoldChange)]

    ggplot(plot_dt, aes(x = rn, y = label, fill = factor(reg_dir))) +
      geom_tile(color = NA) +
      facet_grid(module ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(
        values = c("-1" = "blue", "0" = "white", "1" = "red"),
        name = "Regulation"
      ) +
      labs(x = "Target genes (clustered)", y = "TFs (grouped by module)", fill = "log2FC") +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0),
            axis.text.x = element_text(angle = 45, hjust = 1))
  })

  ## Network graph for Target Regulation
  ## Shows TF->Gene edges (colored by direction) and TF->TF Harmony edges (dashed)
  ## Module-colored TF nodes with rich hover tooltips
  output$regulators <- renderVisNetwork({
    reg_all <- tf_scatter_dt()
    req(nrow(reg_all) > 0)

    reg_all[, label := fifelse(is.na(TF) | TF == "", TF_ID, TF)]

    reg <- copy(reg_all)
    if (!is.null(input$module_select) && input$module_select != "All") {
      reg <- reg[module == as.integer(input$module_select)]
    }
    req(nrow(reg) > 0)

    hit <- unique(reg$TF_ID)
    subnar <- copy(targetsubnar())
    subnar <- subnar[TF_ID %in% hit]
    req(nrow(subnar) > 0)

    ## Nodes
    module_colors <- c(
      "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
      "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
      "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"
    )

    tf_nodes <- reg[, .(
      id = TF_ID, label = label, group = "TF", shape = "triangle",
      value = n_targets,
      color = module_colors[((module - 1) %% length(module_colors)) + 1],
      title = paste0(
        "<b>", label, "</b><br>",
        "TF_ID: ", TF_ID, "<br>",
        "Family: ", Family, "<br>",
        "Module: ", module, "<br>",
        "Module size: ", module_size, "<br>",
        "Module definition: ", input$louvain_weight_mode, "<br><br>",
        "<b>Target regulation</b><br>",
        "Total targets: ", n_targets, "<br>",
        "Up targets: ", n_up_targets, "<br>",
        "Down targets: ", n_down_targets, "<br>",
        "Mean |log2FC|: ", round(mean_abs_l2fc, 3), "<br>",
        "Mean up log2FC: ", round(mean_up_l2fc, 3), "<br>",
        "Mean down log2FC: ", round(mean_down_l2fc, 3), "<br>",
        "Sum up log2FC: ", round(sum_up_l2fc, 3), "<br>",
        "Sum down log2FC: ", round(sum_down_l2fc, 3), "<br><br>",
        "<b>Harmony</b><br>",
        "Outgoing Harmony connections: ", n_harmony_connections_out, "<br>",
        "Incoming Harmony connections: ", n_harmony_connections_in, "<br>",
        "Sum Concordant Harmony (out): ", round(sum_concordant_harmony_out, 3), "<br>",
        "Sum Discordant Harmony (out): ", round(sum_discordant_harmony_out, 3), "<br>",
        "Mean Concordant Harmony (out): ", round(mean_concordant_harmony_out, 3), "<br>",
        "Mean Discordant Harmony (out): ", round(mean_discordant_harmony_out, 3), "<br>",
        "Sum Concordant Harmony (in): ", round(sum_concordant_harmony_in, 3), "<br>",
        "Sum Discordant Harmony (in): ", round(sum_discordant_harmony_in, 3), "<br>",
        "Mean Concordant Harmony (in): ", round(mean_concordant_harmony_in, 3), "<br>",
        "Mean Discordant Harmony (in): ", round(mean_discordant_harmony_in, 3), "<br><br>",
        "Best padj: ", signif(best_padj, 3)
      )
    )]

    gene_ids <- setdiff(unique(subnar$rn), tf_nodes$id)
    gene_nodes <- data.table(
      id = gene_ids, label = gene_ids,
      group = "Gene", shape = "box", value = 1, color = "#d9d9d9"
    )

    nodes <- rbind(tf_nodes, gene_nodes, fill = TRUE)

    ## TF -> Gene edges
    tf_gene_edges <- subnar[, .(
      from = TF_ID, to = rn, arrows = "to",
      value = abs(log2FoldChange),
      width = pmax(1, 3 * abs(log2FoldChange)),
      title = paste0(
        "<b>", fifelse(is.na(TF) | TF == "", TF_ID, TF), " -> ", rn, "</b><br>",
        "Direction: ", ifelse(log2FoldChange > 0, "Positive", "Negative"), "<br>",
        "log2FC: ", round(log2FoldChange, 3), "<br>",
        "padj: ", signif(padj, 3)
      ),
      edge_col = ifelse(
        log2FoldChange > 0,
        paste0("rgba(139,0,0,", input$logtrans, ")"),
        paste0("rgba(0,0,139,", input$logtrans, ")")
      )
    )]
    tf_gene_edges[, color := lapply(edge_col, function(cc) list(color = cc, highlight = cc, hover = cc, inherit = FALSE))]
    tf_gene_edges[, edge_col := NULL]

    ## TF -> TF Harmony edges (directed)
    harmony_cutoff <- as.numeric(input$networkhcutoff)
    k <- 5L
    hdt <- target_harmony_dt()

    h <- hdt[
      TF1_ID %in% hit & TF2_ID %in% hit & TF1_ID != TF2_ID &
        !is.na(Harmony) & abs(Harmony) >= harmony_cutoff,
      .(TF1_ID, TF2_ID, Harmony, Concordant_Harmony, Discordant_Harmony)
    ]

    if (nrow(h) > 0) {
      h <- h[order(TF1_ID, -abs(Harmony))]
      h <- h[, head(.SD, k), by = TF1_ID]

      tf_tf_edges <- h[, .(
        from = TF1_ID, to = TF2_ID, arrows = "to",
        value = abs(Harmony),
        width = pmax(0.1, 0.8 * abs(Harmony)),
        dashes = TRUE,
        smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.15),
        title = paste0(
          "<b>", TF1_ID, " -> ", TF2_ID, "</b><br>",
          "Harmony: ", round(Harmony, 3), "<br>",
          "Concordant Harmony: ", round(Concordant_Harmony, 3), "<br>",
          "Discordant Harmony: ", round(Discordant_Harmony, 3), "<br>",
          "Type: ", ifelse(Harmony >= 0, "Concordant-biased", "Discordant-biased")
        ),
        edge_col = ifelse(
          Harmony >= 0,
          paste0("rgba(217,95,2,", input$harmtrans, ")"),
          paste0("rgba(117,112,179,", input$harmtrans, ")")
        )
      )]
      tf_tf_edges[, color := lapply(edge_col, function(cc) list(color = cc, highlight = cc, hover = cc, inherit = FALSE))]
      tf_tf_edges[, edge_col := NULL]
    } else {
      tf_tf_edges <- data.table(
        from = character(), to = character(), arrows = character(),
        value = numeric(), width = numeric(), dashes = logical(),
        smooth = list(), title = character(), color = list()
      )
    }

    edges <- rbind(tf_gene_edges, tf_tf_edges, fill = TRUE)

    visNetwork(nodes, edges, width = "100%", height = "900px") %>%
      visGroups(groupname = "TF", shape = "triangle") %>%
      visGroups(groupname = "Gene", shape = "box") %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = TRUE, collapse = TRUE
      ) %>%
      visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        stabilization = list(enabled = TRUE, iterations = 800, fit = TRUE)
      ) %>%
      visEvents(stabilizationIterationsDone = "function () { this.setOptions({physics: false}); }")
  })

}

# Run the application 
shinyApp(ui = ui, server = server)