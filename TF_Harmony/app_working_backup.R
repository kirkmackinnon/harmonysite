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
library(gprofiler2)
library(visNetwork)
library(viridis)
library(ggdendro)
library(dendextend)
library(universalmotif)
library(Cairo)
library(shinycssloaders)
library(phylogram)
library(seqinr)
library(ape)


reactiveConsole(TRUE)


### Make narpv from degs: ###
# degs <- fread("../../nebtarget/full_degs_nobatch.tsv")
# tfswithfamilies <- fread("tfidswithfamilies.tsv")
# degs <- merge.data.table(degs, tfswithfamilies, by.x = "TF", by.y = "Name", all.x = T)
# setnames(degs, c("GeneID.x", "GeneID.y"), c("rn", "TF_ID"))
# degs <- merge.data.table(degs, allteds, by.x = "TF_ID", by.y = "Locus", all.x = T)
# fwrite(degs, "newnarpv_nebs_nobatch.tsv", sep = "\t")
#####
options(spinner.type = 5, spinner.color = "darkblue")

testcolors <- c("red", "gray", "blue", "gray", "darkred", "darkblue")
names(testcolors) <- c("Activator", "Minimally Active", "Repressor", "Unknown", "Upregulated", "Downregulated")

htmsa <- Biostrings::readAAMultipleAlignment("Data/msa.fa")

allteds <- fread("Data/allteds.tsv")

#dt <- fread("harmonytable_justnobatch.tsv")
dt <- fread("Data/harmonytable_nobatch_nebs.tsv")

# Rename columns - only for nobatch_nebs.tsv
setnames(dt, c("TF", "Intersect_Concordant", "Intersect_Discordant", "PValue_Concordant", "PValue_Discordant", "Correlation_Concordant", "Correlation_Discordant", "Harmony_Concordant", "Harmony_Discordant"), c("TF1", "Concordant_Intersect", "Discordant_Intersect", "Concordant_PValue", "Discordant_PValue", "Concordant_Correlation", "Discordant_Correlation", "Concordant_Harmony", "Discordant_Harmony"))


dt <- dt[!(TF1 == TF2)]
dt[, TF1 := sub("-.*", "", TF1)]
dt[, TF2 := sub("-.*", "", TF2)]
tfswithfamilies <- fread("Data/tfidswithfamilies.tsv")

setkey(dt, TF1)
setkey(tfswithfamilies, Name)
dt <- dt[tfswithfamilies[, .(Name, TF1_Family = Family)]]
setkey(dt, TF2)
dt <- dt[tfswithfamilies[, .(Name, TF2_Family = Family)]]

pwms <- readRDS("Data/pwms.RDS")
#names(pwms) <- sub(".*?_", "", names(pwms))
pfms <- convert_motifs(pwms, class = "motifStack-pfm")

#narpv <- fread("nar_agg_no-batch_wTED.tsv")
narpv <- fread("Data/newnarpv_nebs_nobatch.tsv")
narpv[, Color := ifelse(sign(log2FoldChange) < 0, "Blue", "Red")]
narpv[TF == "HSF.A4A", TF_ID := "AT4G18880"]

exp17nogo <- narpv[EXP == "1-17", unique(TF)]

narpv <- narpv[EXP != "1-17"]

allgeneids <- unique(narpv$rn)

gis <- narpv[, .(unique(rn))]
vertices <- merge(gis, allteds, by.x = "V1", by.y = "Locus", all.x = T)
vertices <- merge(vertices, unique(narpv[, .(TF_ID, TF)]), by.x = "V1", by.y = "TF_ID", all.x = T)
#vertices[!(is.na(TF)), V1 := TF]

cte <- fread("Data/cellTypeExpression.tsv", key = "Gene ID")
jitr <- fread("Data/JITGenes_root.csv", skip = 1, select = c("Gene", "FDR adjusted p-value", "First Response (Just-in-time bin)"), col.names = c("GeneID", "pvalue", "JIT"))
setkey(jitr, GeneID)


jits <- fread("Data/JITGenes_shoot.csv", skip = 1, select = c("AtID", "FDR adjusted p-value", "First Response (Just-in-time bin)"), col.names = c("GeneID", "pvalue", "JIT"))
setkey(jits, GeneID)

dt <- dt[!(TF1 %in% exp17nogo)]
dt <- dt[!(TF2 %in% exp17nogo)]

#idoptions <- dt[, unique(sub("AT.*?_", "", V1))]
#idoptions <- sub("_DES...*.csv", "", list.files("./DEGS", pattern = "^A.*"))

### modified when nobatch nebs were added because they were TSVs instead
idoptions <- sub("_DES...*.tsv", "", list.files("Data/DEGS", pattern = "^A.*"))

ids <- as.data.table(idoptions)
ids[, c("Ids", "TF") := tstrsplit(idoptions, "_")]
idoptions <- sub("^A.*?_", "", idoptions)
idoptions <- narpv[order(TF), unique(TF)]
newdt <-
  dt[, .(
    TF1 = sub("AT.*?_", "", TF1),
    TF2 = sub("AT.*?_", "", TF2),
    Concordant = Concordant_Harmony,
    Discordant = Discordant_Harmony
  )][order(-Concordant)]

fulltfids <- fread("Data/all_ath_tf_ids.tsv")

vertices <- merge(vertices, fulltfids, by.x = "V1", by.y = "Gene_ID", all.x = T)
vertices[!(is.na(Family)), shape := "triangle"]
vertices[(is.na(Family)), shape := "circle"]
vertices[!(is.na(TF)), shape := "triangle"]
#vertices[!(is.na(TF)), V1 := TF]
vertices <- unique(vertices)
# Define UI for application that draws a histogram
ui <- navbarPage("Landscape of TF Harmony",
                 
                 
          tabPanel("Global Analyses",

                     inputPanel(
                       
                       textInput(inputId = "familyharmonycutoff", label = "Harmony Minimum:", value = 0),
                       textInput(inputId = "familypvcutoff", label = "Fisher P-Value Maximum:", value = 0.05),
                       textInput(inputId = "familyintersectcutoff", label = "Intersect Minimum:", value = 0),
                       
                     selectInput(
                       inputId = "distmeth", 
                       label = "Distance Method", choices = c("euclidean", "maximum", "manhattan", "canberra"),
                       selected = "maximum", 
                       multiple = F
                     ),
                     selectInput(
                       inputId = "clustmeth", 
                       label = "Clustering Method", 
                       choices = c("complete", "average", "ward.D", "single", "mcquitty", "median", "centroid", "ward.D2"),
                       selected = "ward.D2", 
                       multiple = F
                     ),
                     selectInput(inputId = "motifcompmethod", label = "Motif Comparison Method:", choices = c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT", "HELL", "SEUCL", "MAN", "ALLR_LL", "WEUCL", "WPCC"), selected = "ALLR_LL", multiple = F, selectize = F)
                     
                   
                     
                   ),
                   selectInput(
                     inputId = "family1", 
                     label = "Search by TF or TF Family:", 
                     choices = unique(c(tfswithfamilies$Family, idoptions)),
                     selected = unique(tfswithfamilies$Family),
                     multiple = T,
                     width = "100%",  
                     selectize = T
                   ),
                   
                   tabsetPanel(
                     
                     tabPanel("Harmony Heatmap",
                              tabsetPanel(
                                tabPanel("Concordant Harmony",
                              plotlyOutput(
                                outputId = "conHarmonyHeatmap",
                                
                                height = 2000
                              ) %>% withSpinner(), ),
                              tabPanel("Discordant Harmony",
                              plotlyOutput(
                                outputId = "disHarmonyHeatmap",
                                
                                height = 2000
                              )%>% withSpinner(),) 
                              ),
                     ),
                     tabPanel("Family Harmony", 
                              splitLayout(
                                plotlyOutput( outputId= "relationships", 
                                              height = 1000) %>% withSpinner(),
                                visNetworkOutput( outputId= "familyproportion",
                                                  height = 1000) %>% withSpinner()
                                
                              ),
                     ),
                     tabPanel("Tanglegram",
                              inputPanel(
                                selectInput(
                                  inputId = "tanglechoice1", 
                                  label = "Left side:", 
                                  choices = c("Concordant", "Discordant", "Phylogeny", "Up Motif similarity", "Down Motif similarity"),
                                  selected = "Concordant", 
                                  multiple = F
                                  
                                  ), 
                                selectInput(
                                  inputId = "tanglechoice2", 
                                  label = "Right side:", 
                                  choices = c("Concordant", "Discordant", "Phylogeny", "Up Motif similarity", "Down Motif similarity"),
                                  selected = "Discordant", 
                                  multiple = F
                                ),
                                selectInput(
                                  inputId = "kbreaks",
                                  label = "Number of breaks:",
                                  choices = c(0,1,2,3,4,5,6,7,8,9),
                                  selected = 1, 
                                  multiple = F
                                )
                              ), 
                              plotOutput(outputId = "tanglegram", height = 1000) %>% withSpinner()

                              

                     ),
                     tabPanel("Harmony Table", 
                                DT::dataTableOutput(outputId = "fullharmonydt")
                              ),
                     
                     tabPanel("Motif Heatmap",
                            
                              plotlyOutput(
                                outputId = "motifheatmap", height = 1000
                                )%>% withSpinner()
                              ),
                     tabPanel("Motif Dendrogram",
                              plotlyOutput(
                                outputId = "motifdend", height = 1000
                              )%>% withSpinner()
                              ),
                     # tabPanel("Browse Motifs",
                     #          plotOutput(
                     #            outputId = "motifbrowser",
                     #            height = 1000
                     #          ) %>% withSpinner()
                     #          )
                              
                   )
          ),
          
          tabPanel("Pairwise Analyses",
                   
                     inputPanel(
                       selectizeInput(
                         inputId = "TFs",
                         label = "TFs",
                         choices = idoptions,
                         multiple = T,
                         selected = idoptions[1:2],
                         options = list(maxItems = 10)
                       ), 
                       fileInput(inputId = "tfupload", label = "Upload list of TFs"),

                       
                       sliderInput(
                         inputId = "harmonyplotheight",
                         min = 500,
                         max = 10000,
                         step = 100,
                         value = 1000,
                         label = "Plot height:"

                       ),
                       checkboxInput(inputId = "harmonyinterfacet", label = "Facet by interaction:", value = T),
                     ),
                   tabsetPanel(id = "pairwise",
                   tabPanel("TF x TF",
                   
                   splitLayout(
                        plotlyOutput(outputId = "harmonyplotly", height = "100%") %>% withSpinner(),
                        DT::dataTableOutput("harmonyTable")
                   )

                   
          ),
          tabPanel("Motif Comparison", value = "pairwisemotifs",
                   selectInput("motiftreestyle", label = "Tree Style:", choices = c("stack", "tree", "radialPhylog"), multiple = F),
                   #uiOutput(outputId = "tfmotif"),
                   imageOutput(
                     outputId = "motifplot", 
                     width = "80%", 
                     height = 800
                   )%>% withSpinner()
                   
                   
          ),
          tabPanel("Cell Type Specificity",
                   plotlyOutput(
                     outputId = "ctplot",  
                     height = 1000
                   )
          ),
          tabPanel("NxTime Overlay", 
                   splitLayout(
                     cellWidths = c("50%", "50%"), 
                     plotlyOutput(
                       outputId = "nxtplot1", 
                       height = 1000
                     ), 
                     plotlyOutput(
                       outputId = "nxtplot2", 
                       height = 1000
                     )
                   )
          ), 
          tabPanel("GO Terms",
                   sliderInput(inputId = "goplotheight", min = 500, max = 10000, step = 100, value = 1000, label = "Plot height:"),
                   splitLayout(
                   plotlyOutput(
                     outputId = "goplot", 
                     width = "100%",
                     height = 1000
                   ),
                   plotOutput(
                     outputId = "gotable",
                     width = "100%", 
                     height = 1000
                   ), style = "overflow:auto;" 
                   )
                   
          ),
          tabPanel("TF Regulation",
                   inputPanel(
                   sliderInput(
                     inputId = "HarmonyRange", 
                     label = "Harmony Cutoff:", 
                     min = 0, 
                     max = 215000, 
                     value = 0, 
                     round = T, 
                     animate = T, 
                     width = 800
                     ),
                   selectInput(
                     inputId = "tfregnetworkstyle",
                     label = "Network Style:",
                     choices = c(
                       "layout_with_sugiyama",
                       "layout_on_sphere",
                       "layout_with_fr",
                       "layout_with_gem",
                       "layout_with_kk",
                       "layout_with_lgl",
                       "layout_with_mds",
                       "layout_in_circle",
                       "layout_nicely",
                       "layout_as_star",
                       "layout_as_tree",
                       "layout_on_grid",
                       "layout_with_dh"
                     ),
                     selected = "layout_with_sugiyama",
                     multiple = F
                   ),
                   
                   checkboxInput(
                     inputId = "networkdegs", 
                     label = "DEGs only:", 
                     value = TRUE
                     ), 
                   ),
            visNetworkOutput(
              outputId = "networkplot", 
              height = 1000
              )%>% withSpinner()
          ), 
          tabPanel("DEG Networks",
                   inputPanel(
                     textInput(inputId = "tfregpcutoff", label = "Adjusted P-value Cutoff:", value = "0.05"),
                     textInput(inputId = "tfreglfccutoff", label = "Log2FoldChange Cutoff:", value = "0"),
                     textInput(inputId = "tfdegreecutoff", label = "Degree Cutoff:", value = "0"),
                     selectInput(
                       inputId = "degnetworkstyle",
                       label = "Network Style:",
                       choices = c(
                         "layout_with_sugiyama",
                         "layout_on_sphere",
                         "layout_with_fr",
                         "layout_with_gem",
                         "layout_with_kk",
                         "layout_with_lgl",
                         "layout_with_mds",
                         "layout_in_circle",
                         "layout_nicely",
                         "layout_as_star",
                         "layout_as_tree",
                         "layout_on_grid",
                         "layout_with_dh"
                       ),
                       selected = "layout_with_sugiyama",
                       multiple = F
                     )
                   
                   ),
                   visNetworkOutput(
                     outputId = "tfnetworkplot", 
                     height = 1000
                   )%>% withSpinner()
            
          ),),),
          
          tabPanel("Target Regulation",
                   inputPanel(
                     fileInput(inputId = "targetupload", label = "Upload list of Targets"),
                     textInput(inputId = "regpcutoff", label = "Adjusted P-value Cutoff:", value = "0.05"),
                     textInput(inputId = "reglfccutoff", label = "Log2FoldChange Cutoff:", value = "0"),
                     selectInput(
                       inputId = "targetnetworkstyle",
                       label = "Network Style:",
                       choices = c(
                         "layout_with_sugiyama",
                         "layout_on_sphere",
                         "layout_with_fr",
                         "layout_with_gem",
                         "layout_with_kk",
                         "layout_with_lgl",
                         "layout_with_mds",
                         "layout_in_circle",
                         "layout_nicely",
                         "layout_as_star",
                         "layout_as_tree",
                         "layout_on_grid",
                         "layout_with_dh"
                       ),
                       selected = "layout_with_sugiyama",
                       multiple = F
                     )
                   ),
                   selectInput(
                     inputId = "allgenes",
                     label = "Targets: ",
                     choices = allgeneids,
                     multiple = T, 
                     width = "100%",
                     selected = allgeneids[1]
                   ),
                   splitLayout(
                   visNetworkOutput(
                    outputId = "regulators",
                    width = "50%",
                    height = 1000
                   )%>% withSpinner(),
                   verticalLayout(
                   DT::DTOutput(outputId = "tfregdt", height = 500, width = "50%"),
                   DT::DTOutput(outputId = "targetregdt", height = 500, width = "50%")
                   )
                   ),
          includeCSS("~/Desktop/github/HTT_144TFs/TF_Harmony/www/style.css")

          )
          
        )
        


ftestfun <- function(shared, Tar1, Tar2){
  totalgenes <- 32031
  inNewOnly <- nrow(Tar1) - shared
  inOldOnly <- nrow(Tar2) - shared
  notInEither <- totalgenes - (shared + inNewOnly + inOldOnly)
  
  contTable <- data.table(
    notA = c(notInEither, inNewOnly),
    inA = c(inOldOnly, shared))

  #print(contTable)
  
  test <- fisher.test(contTable, alternative = "greater" )
  #print(test)
  return(test)
}

findfile <- function(pat, oldll = list()){
  lllength <- length(pat) + length(oldll)
  
  ll <- vector("list", length(pat))
  

  for (item in seq_along(pat)) {
    if (pat[item] %in% names(oldll)) next
    ll[[item]] <- narpv[TF == pat[item]][padj < 0.05, .(rn, log2FC = log2FoldChange, padj)]
    names(ll)[item] <- pat[item]
  }
    
  ll <- c(ll, oldll)
  ll <- ll[names(ll) %in% pat]
  
  
  return(ll)
}


matchtidy <- function(Tar1, oldmatch = data.table()){
  #matched <- Reduce(function(...) merge.data.table(..., by = "rn", suffixes = paste0("_", names(Tar1))), Tar1)
  matched <- data.table()
  combos <- c()
  listlength <- choose(length(Tar1), 2)
  tempagg <- vector("list", listlength)
  count <- 1
  for (i in 1:length(Tar1)){
    for (j in i:length(Tar1)){
      if (i == j) next
      inter <- paste0(names(Tar1)[i], "_", names(Tar1)[j])
      combos <- c(combos, inter)
      if (inter %in% oldmatch$Inter) next
      temp <- merge.data.table(Tar1[[i]], Tar1[[j]], by = "rn")
      temp[, Inter := inter]
      tempagg[[count]] <- temp
      count <- count + 1
      #matched <- rbind(matched, temp)
    }
    
  }
  #matched <- do.call(rbind, tempagg)
  matched <- rbindlist(tempagg)
  matched[, Harmony := ifelse(sign(`log2FC.x`) == sign(`log2FC.y`), "Concordant", "Discordant")]
  matched <- rbind(matched, oldmatch)
  matched <- matched[Inter %in% combos]
  return(matched)
}

motifloader <- function(fn){
  return(tags$figure(
    class = "motif", 
    tag$img(
      src = fn,
      height = 250,
      alt = fn
    ),
    tags$figcaption(fn)
  ))
}

# Define server logic required to draw a histogram
  server <- function(input, output, session) {
  
  Tar1 <- reactive({findfile(input$TF1)})
  Tar2 <- reactive({findfile(input$TF2)})
  
  #print(input$TFs)
  tarfs <- reactive({
    if (exists("tarfs()")) findfile(input$TFs, tarfs())
    else findfile(input$TFs)
    })
  matched <- reactive({
    if (exists("matched()")) matchtidy(tarfs(), matched())
    else matchtidy(tarfs())
    })
  
  subgraph <- reactive({
    #newdt[(TF1 == input$TF1 | TF2 == input$TF2 | TF1 == input$TF2 | TF2 == input$TF1) & !(is.na(Concordant) & is.na(Discordant))]
    newdt[(TF1 %in% input$TFs)][TF2 %in% input$TFs]#[!(is.na(Concordant) & is.na(Discordant))]
    })
  
  
  
  shared <- nrow(matched())
  
  # concordantshared <- reactive({nrow(matched()[Harmony == "Concordant"])})
  # discordantshared <- reactive({nrow(matched()[Harmony == "Discordant"])})
  # 
  # concorrelation <- reactive({cor(matched()[Harmony == "Concordant", log2FC.x], matched()[Harmony == "Concordant", log2FC.y])})
  # discorrelation <- reactive({cor(matched()[Harmony == "Discordant", log2FC.x], matched()[Harmony == "Discordant", log2FC.y])})
  # 
  # conftest <- reactive({ftestfun(concordantshared(), Tar1(), Tar2())})
  # disftest <- reactive({ftestfun(discordantshared(), Tar1(), Tar2())})
  
  trigger_motif_redraw <- reactiveVal(0)
  
  observeEvent(input$tabs, {
    if (input$tabs == "pairwise") {
      trigger_motif_redraw(trigger_motif_redraw() + 1)
    }
  })

  # shinyjs::delay(50, {
  output$motifplot <- renderImage({
    req(input$pairwise == "pairwisemotifs")
    
    pattern <- paste0("^(", paste0(input$TFs, collapse = "|"), ")_")
    
    subpfms <- pfms[grepl(pattern, names(pfms))]
    
    grid::grid.newpage()
    #subpfms <- pfms[sub("_.*", "", names(pfms)) %in% input$TFs]
    tryCatch({
      
      outfile <- tempfile(fileext = ".png")
      
      svg(outfile, width = 10, height = 5)
      motifStack(subpfms, layout = input$motiftreestyle)
      dev.off()
      list(src = outfile,
           contentType = 'image/svg+xml',
           width = "100%",  # SVG scales naturally
           height = "auto",
           alt = "Motif stack plot")
    
       
    }, error = function(e) {
      print(paste("Motif plot error:", e$message))
    })
    
  }, deleteFile = TRUE)
  # })
  output$motifbrowser <- renderPlot({
    pattern <- paste0("^(", paste0(input$TFs, collapse = "|"), ")_")
    
    subpfms <- pfms[grepl(pattern, colnames(motifmat()))]
    #subpfms <- pfms[names(pfms) %in% colnames(motifmat())]
    #phylog <- ade4::hclust2phylog(motifhclust())
    #browseMotifs(subpfms, phylog = phylog)
    motifStack(subpfms, reorder = F)
    # subpwms <- pwms[names(pwms) %in% colnames(motifmat())]
    # motif_tree(subpwms, layout = "rectangular")
  })
  
    output$harmonyTable <- DT::renderDataTable(
      {
        setDT(subgraph())
      }, 
      selection = 'single',
      options = list(
        scrollY = '400px',
        paging = FALSE
        
      )
    )  
    
    subdt <- reactive({
      subdt <-
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
    
    output$fullharmonydt <- DT::renderDataTable(
      {
        
        setDT(subdt())
      }, 
      selection = 'single',
      options = list(
        scrollY = '800px',
        paging = FALSE
        
      )
    )  
    
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
      #subcte$ID <- factor(subcte$ID, levels = row.names(submat)[fit$order], ordered = T)
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
    
    goout <- reactive({
      targs <- split(matched()$rn, f = matched()$Inter)
      gost(targs, organism = "athaliana", multi_query = T)
    })
    
    output$goplot <- renderPlotly({
      
      ggplotly(gostplot(goout(), interactive = F), height = input$goplotheight)
      
    })
    
    output$gotable <- renderPlot({
      #print(head(goout()))
      publish_gosttable(goout(), ggplot = T)
    })
    
    outputOptions(output, "goplot", suspendWhenHidden = FALSE)
    
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
    observeEvent(input$targetupload, {
      file <- input$targetupload
      fdt <- fread(file$datapath, header = F)
      
      updateSelectInput(
        session = session,
        inputId = "allgenes",
        choices = allgeneids,
        selected = fdt[V1 %in% allgeneids, V1]
      )
    })
    # 
    # observeEvent(input$harmonyTable_rows_selected, {
    #   updateSelectInput(
    #     session = session,
    #     inputId = "TF2",
    #     choices = idoptions,
    #     selected = newdt[input$harmonyTable_rows_selected, TF2]
    #   )
    # })
    
    proxy <- dataTableProxy(outputId = "harmonyTable")
    
    observeEvent(input$TFs, {
      meltdt <- melt.data.table(subgraph())
      meltdt <- meltdt[!(is.na(value) | is.infinite(value))]
      proxy %>% selectRows(newdt[TF1 %in% input$TFs & TF2 %in% input$TFs, which = T],)
      updateSliderInput(session, inputId = "HarmonyRange", value = 0, min = 0, max = max(abs(meltdt$value)))
    })
    
    observeEvent(input$TF2, {
      meltdt <- melt.data.table(subgraph())
      meltdt <- meltdt[!(is.na(value) | is.infinite(value))]
      proxy %>% selectRows(newdt[TF1 == input$TF1 & TF2 == input$TF2, which = T],)
      updateSliderInput(session, inputId = "HarmonyRange", value = 0, min = 0, max = max(meltdt$value))
    })
    
    output$tfnetworkplot <- renderVisNetwork({
      padjcutoff <- as.numeric(input$tfregpcutoff)
      l2fccutoff <- as.numeric(input$tfreglfccutoff)
      degreecutoff <- as.numeric(input$tfdegreecutoff)
      
      subnar <- narpv[TF %in% input$TFs][padj < padjcutoff][abs(log2FoldChange) > l2fccutoff]
      
      
      subnar[, group := ifelse(sign(log2FoldChange) == 1, "Upregulated", "Downregulated")]
      subnar <- subnar[!(TF_ID == rn)]
      subnar <- subnar[order(-Activity, Color, log2FoldChange)]
      #mygraph <- as_tbl_graph(subnar[,.(TF, rn, padj, value = log2FoldChange, group)], directed = T)
      
      #print(head(subnar))
      
      mygraph <- graph_from_data_frame(
        d = subnar[, .(TF_ID, rn, label = TF, value = log2FoldChange, padj, group)], 
        vertices = unique(vertices[(V1 %in% subnar$TF) | (V1 %in% subnar$rn) | (V1 %in% subnar$TF_ID), .(V1, label = TF, group = Activity, value = TedStrength, shape = shape)])
        )
      
      
      
      V(mygraph)$Degree <- igraph::degree(mygraph)
      
      mygraph <- delete.vertices(mygraph, V(mygraph)$Degree < degreecutoff)
      
      #V(mygraph)$group <- ifelse(names(V(mygraph)) %in% narpv$TF, narpv[TF %in% names(V(mygraph)),  Activity], "Unknown")
      #V(mygraph)$value <- ifelse(names(V(mygraph)) %in% narpv$TF, narpv[TF %in% names(V(mygraph)),  TedStrength], 0)
      
      V(mygraph)$color <- testcolors[V(mygraph)$group]
      E(mygraph)$color <- testcolors[E(mygraph)$group]
      
      #print(head(mygraph))
      
      visIgraph(mygraph, layout = input$degnetworkstyle) %>% 
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), collapse = T ) %>% 
        visNodes(color = list(border = "black", highlight = "yellow")) %>%
        visGroups(groupname = "Activator", color = "red", shape = "triangle") %>%
        visGroups(groupname = "Minimally Active", color = "gray", shape = "triangle") %>%
        visGroups(groupname = "Repressor", color = "blue", shape = "triangle") %>%
        visGroups(groupname = "Unknown", color = "gray") %>%
        # visGroups(groupname = "Upregulated", color = "darkred", shape =) %>%
        # visGroups(groupname = "Downregulated", color = "darkblue") %>%
        visLegend() 
        
    })
    
  
    output$networkplot <- renderVisNetwork({
      
      melted <- melt.data.table(subgraph(), variable.name = "group")
      melted <- melted[!(is.na(value) | is.infinite(value)) & abs(value) > input$HarmonyRange ]
      melted[, color := ifelse(group == "Concordant", "red", "blue")]
      
      if (input$networkdegs){
        
        tf2tfs <- melted$TF2
        tf2ids <- c()
        #t2test <- c()
        for (item in tarfs()){
          tf2ids <- c(tf2ids, ids[TF %in% tf2tfs][(Ids %in% item$rn)])
          #t2test <- c(t2test, paste0(rep(item$rn), ids[TF %in% tf2tfs][(Ids %in% item$rn)], sep = "_"))
        }
        #tf2ids <- ids[TF %in% tf2tfs][(Ids %in% tarfs())]
        #print(t2test)
        melted <- melted[TF2 %in% tf2ids$TF]
        
      }
      
      
      
      mygraph <- as_tbl_graph(as.data.table(melted), directed = T)
      
      V(mygraph)$Degree <- igraph::degree(mygraph)
      
      #print(mygraph)
      
      visIgraph(mygraph, layout = input$tfregnetworkstyle) %>% 
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), collapse = T, width = "100%") %>%
        #visEdges(arrows = "to" ) %>%
        visNodes(color = list(border = "black", highlight = "yellow")) %>%
        visGroups(groupname = "Concordant", color = "red") %>%
        
        visGroups(groupname = "Discordant", color = "blue") %>%
        
     
        visLegend() 
      
      
      # ggraph(mygraph, "dendrogram", circular = T) + 
      #   geom_conn_bundle(data = get_con(from = from, to = to), alpha = 0.1) 
      
     # ggraph(mygraph, "stress") +
     # 
     #    geom_edge_parallel2(
     #      aes(
     #        color = variable,
     #        fill = variable,
     #        width = abs(value),
     #      ),
     #      arrow = arrow(
     #        #angle = ifelse(sign(log2FC) == 1, 30, 90),
     #        type = "closed"
     #      ),
     #      
     #      start_cap = circle(5, 'mm'),
     #      end_cap = circle(5, 'mm'),
     #      #alpha = 0.6,
     #      strength = 0.2
     #    ) +
     #    # geom_edge_density(
     #    #   aes(fill = variable)
     #    # )+
     #    geom_node_text(
     #      aes(label = name),
     #      size = 10,
     #      repel = T,
     #      #nudge_x = -0.1
     #    ) +
     # 
     #    geom_node_point(
     #      aes(size = Degree)
     #    ) +
     # 
     #    theme_void(base_size = 20) +
     # 
     #    theme(
     #    #  legend.position = "left",
     #      legend.position = "bottom",
     #      legend.direction = "vertical",
     # 
     #    ) +
     #    # 
     #    # scale_edge_alpha(range = c(0.3,1)) +  
     #    # 
     #    scale_edge_width(range = c(1, 7)) +
     #    # 
     #    scale_edge_color_manual(
     #      values = viridis::viridis(2, end = 0.75, direction = -1) #,
     #      #breaks = c("Negative", "Positive")
     #    ) #+
     #    # 
     #    # labs(
     #    #   width = "Harmony Strength",
     #    #   size = "Degree Connectivity",
     #    #   color = "\nHarmony Direction"
     #    # )  
      
      # simpleNetwork(
      #   melted,
      #   zoom = T
      #   )
    })
    
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
        
        geom_smooth(aes(group = interaction(Inter, Harmony)), method = "lm", se = F, size = 1) +
        
        geom_vline(xintercept = 0) +
        
        geom_hline(yintercept = 0) + 
        
        geom_line(aes(group = rn), color = "black", alpha = 0.25)
        
       ggplotly(p, height = input$harmonyplotheight)
    })
    
    
    output$relationships <- renderPlotly({
      subdt <- dt[(((TF1_Family %in% input$family1) & (TF2_Family %in% input$family1)) | ((TF1 %in% input$family1) & (TF2 %in% input$family1))), .(TF1, TF2, Concordant_Harmony, Discordant_Harmony, TF1_Family, TF2_Family)]
      subdt <- subdt[!(is.na(Discordant_Harmony)) &
                       (abs(Discordant_Harmony) >= as.numeric(input$familyharmonycutoff))]
      subdt <- subdt[!(is.na(Concordant_Harmony)) &
                       (abs(Concordant_Harmony) >= as.numeric(input$familyharmonycutoff))]
      # [
      #   Discordant_Intersect > as.numeric(input$familyintersectcutoff)][
      #     Discordant_PValue < as.numeric(input$familypvcutoff)][
      #       !(is.infinite(Discordant_Harmony))
      #     ]
      subdt[, Interaction := paste(pmin(TF1_Family, TF2_Family), pmax(TF1_Family, TF2_Family), sep = ".")]
      
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
    
    output$familyproportion <- renderVisNetwork({
      subdt <- dt[(((TF1_Family %in% input$family1) &
                      (TF2_Family %in% input$family1)
      ) |
        ((TF1 %in% input$family1) &
           (TF2 %in% input$family1))), .(TF1,
                                         TF2,
                                         Concordant_Harmony,
                                         Discordant_Harmony,
                                         TF1_Family,
                                         TF2_Family)]
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
      
      # Sugiyama layout
      sugiyama_layout <- layout_with_sugiyama(g)$layout
      
      # Add layout to nodes (scale/flip for visNetwork)
      nodes[, x := sugiyama_layout[, 1] * 100]
      nodes[, y := -sugiyama_layout[, 2] * 100]
      
      # -------------------------
      # 4. Plot in visNetwork
      # -------------------------
      visNetwork(nodes, edges) %>%
        visNodes(shape = "dot", size = 25, fixed = FALSE) %>%
        visEdges(smooth = TRUE) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visPhysics(enabled = FALSE)  # Respect igraph layout
      
      
      # # Create edges for Concordant Harmony
      # edges_concordant <- family_harmony[mean_H_concordant > 0.5, .(
      #   from = TF1_Family,
      #   to = TF2_Family,
      #   value = mean_H_concordant * 10,
      #   title = paste0("Concordant Harmony: ", round(mean_H_concordant, 3)),
      #   color = "blue",
      #   arrows = "to"
      # )]
      # 
      # 
      # # Create edges for Discordant Harmony (reverse direction)
      # edges_discordant <- family_harmony[mean_H_discordant > 0.5, .(
      #   from = TF2_Family,  # Reverse direction
      #   to = TF1_Family,
      #   value = mean_H_discordant * 10,
      #   title = paste0("Discordant Harmony: ", round(mean_H_discordant, 3)),
      #   color = "red",
      #   arrows = "to"
      # )]
      # 
      # # Combine edges
      # edges <- rbind(edges_concordant, edges_discordant)
      # 
      # 
      # 
      # 
      # # Prepare nodes
      # nodes <- unique(c(edges$from, edges$to))
      # nodes <- data.table(id = nodes, label = nodes)
      # 
      # 
      # 
      # tryCatch({
      # visNetwork(nodes, edges) %>%
      #   visEdges(smooth = TRUE) %>%  # smoother layout for overlapping edges
      #   visNodes(shape = "dot", size = 25) %>%
      #   visHierarchicalLayout(direction = "LR", sortMethod = "directed") %>%
      #   #visLayout(randomSeed = 42) %>%
      #   visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      #   visPhysics(stabilization = TRUE) #%>%
      #   #visInteraction(navigationButtons = TRUE) %>%
      #   # visLegend(addEdges = list(
      #   #   list(label = "Concordant Harmony", color = "blue"),
      #   #   list(label = "Discordant Harmony", color = "red")
      #   # ))
      # }, error = function(e){
      #   message("Error rendering visNetwork:")
      #   message(e)
      #   return(NULL)
      # }
      # )
      
      # # Prepare edges (only show those above a threshold to reduce clutter)
      # edges <- family_harmony[prop_H_concordant > 0.01, .(
      #   from = TF1_Family,
      #   to = TF2_Family,
      #   value = prop_H_concordant * 10,  # scale for visualization
      #   title = paste0("Harmony: ", round(prop_H_concordant, 3))
      # )]
      # 
      # # Create unique node list from both `from` and `to`
      # nodes <- unique(c(edges$from, edges$to))
      # nodes <- data.table(id = nodes, label = nodes)
      # 
      # visNetwork(nodes, edges) %>%
      #   visEdges(arrows = "to", smooth = FALSE) %>%
      #   visNodes(shape = "dot", size = 20) %>%
      #   visLayout(randomSeed = 42) %>%
      #   visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      #   visPhysics(stabilization = TRUE) %>%
      #   visInteraction(navigationButtons = TRUE) %>%
      #   visLegend()
    })
    
    conhclustx <- reactive({
      subdtcon <- dt[(((TF1_Family %in% input$family1) & (TF2_Family %in% input$family1)) | ((TF1 %in% input$family1) & (TF2 %in% input$family1))), .(TF1, TF2, Concordant_Intersect, Concordant_PValue, Concordant_Harmony, TF1_Family, TF2_Family)]

      subdtcon <- subdtcon[!(is.na(Concordant_Harmony))][
        (abs(Concordant_Harmony) > as.numeric(input$familyharmonycutoff))][
          Concordant_Intersect > as.numeric(input$familyintersectcutoff)][
            Concordant_PValue < as.numeric(input$familypvcutoff)][
              !(is.infinite(Concordant_Harmony))
            ]
      

      subdtconwdt <- dcast.data.table(subdtcon[, .(TF1, TF2, Concordant_Harmony)], TF1 ~ TF2, value.var = "Concordant_Harmony", fill = 0)

      conmat <- as.matrix(subdtconwdt, rownames = "TF1")

      
      
      condtmx <- dist(conmat, method = input$distmeth)

      #condtmy <- dist(t(conmat))
      #disdtmy <- dist(t(dismat))
      
      conhclustx <- hclust(condtmx, method = input$clustmeth)

      #conhclusty <- hclust(condtmy)
      #dishclusty <- hclust(disdtmy)
      
      #conorderx <- conhclustx$labels[conhclustx$order]
    })
    
    conhclusty <- reactive({
      subdtcon <- dt[(((TF1_Family %in% input$family1) & (TF2_Family %in% input$family1)) | ((TF1 %in% input$family1) & (TF2 %in% input$family1))), .(TF1, TF2, Concordant_Intersect, Concordant_PValue, Concordant_Harmony, TF1_Family, TF2_Family)]
      
      subdtcon <- subdtcon[!(is.na(Concordant_Harmony))][
        (abs(Concordant_Harmony) > as.numeric(input$familyharmonycutoff))][
          Concordant_Intersect > as.numeric(input$familyintersectcutoff)][
            Concordant_PValue < as.numeric(input$familypvcutoff)][
              !(is.infinite(Concordant_Harmony))
            ]
      
      
      subdtconwdt <- dcast.data.table(subdtcon[, .(TF1, TF2, Concordant_Harmony)], TF1 ~ TF2, value.var = "Concordant_Harmony", fill = 0)
      
      conmat <- as.matrix(subdtconwdt, rownames = "TF1")
      
      
      
      condtmy <- dist(t(conmat), method = input$distmeth)
      
      #condtmy <- dist(t(conmat))
      #disdtmy <- dist(t(dismat))
      
      conhclusty <- hclust(condtmy, method = input$clustmeth)
    })
    
    dishclustx <- reactive({

      subdtdis <- subdt()[, .(TF1, TF2, Discordant_Intersect, Discordant_PValue, Discordant_Harmony)]
      
      
      
      subdtdis <- subdtdis[!(is.na(Discordant_Harmony)) & (abs(Discordant_Harmony) > as.numeric(input$familyharmonycutoff))][
        Discordant_Intersect > as.numeric(input$familyintersectcutoff)][
          Discordant_PValue < as.numeric(input$familypvcutoff)][
            !(is.infinite(Discordant_Harmony))
          ]
      
      subdtdiswdt <- dcast.data.table(subdtdis[, .(TF1, TF2, Discordant_Harmony)], TF1 ~ TF2, value.var = "Discordant_Harmony", fill = 0)
      
      dismat <- as.matrix(subdtdiswdt, rownames = "TF1")
      
      
      
      disdtmx <- dist(dismat, method = input$distmeth)
      
      #condtmy <- dist(t(conmat))
      #disdtmy <- dist(t(dismat))
      
      dishclustx <- hclust(disdtmx, method = input$clustmeth)
      
      #conhclusty <- hclust(condtmy)
      #dishclusty <- hclust(disdtmy)
      
      
    })
    
    dishclusty <- reactive({
      
      subdtdis <- subdt()[, .(TF1, TF2, Discordant_Intersect, Discordant_PValue, Discordant_Harmony)]
      
      
      
      subdtdis <- subdtdis[!(is.na(Discordant_Harmony)) & (abs(Discordant_Harmony) > as.numeric(input$familyharmonycutoff))][
        Discordant_Intersect > as.numeric(input$familyintersectcutoff)][
          Discordant_PValue < as.numeric(input$familypvcutoff)][
            !(is.infinite(Discordant_Harmony))
          ]
      
      subdtdiswdt <- dcast.data.table(subdtdis[, .(TF1, TF2, Discordant_Harmony)], TF1 ~ TF2, value.var = "Discordant_Harmony", fill = 0)
      
      dismat <- as.matrix(subdtdiswdt, rownames = "TF1")
      
      
      
      disdtmy <- dist(t(dismat), method = input$distmeth)
      
      #condtmy <- dist(t(conmat))
      #disdtmy <- dist(t(dismat))
      
      dishclusty <- hclust(disdtmy, method = input$clustmeth)
      
      #conhclusty <- hclust(condtmy)
      #dishclusty <- hclust(disdtmy)
      
      
    })
    
    
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
      #conordery <- conhclustx()$labels[conhclust]
      
      # allcols <- viridis(length(input$family1))
      # names(allcols) <- unique(input$family1)
      # subfam <- tfswithfamilies[Name %in% unique(c(subdtcon$TF1, subdtcon$TF2))]
      # subfam[, Colors := allcols[Family]]
      # 
      #ycols <- tfswithfamilies[Name %in% unique(c(subdtcon$TF1, subdtcon$TF2))]
      #xcols <- axiscolors[unique(subdtcon$TF2)]
      
      #print(length(conorder) == length(unique(c(subdtcon$TF1, subdtcon$TF2))))
      
       subdtcon$TF1 <- factor(subdtcon$TF1, levels = conorderx)
       subdtcon$TF2 <- factor(subdtcon$TF2, levels = conordery)
       setkey(subdtcon, TF1)
       #print(head(subdtcon$TF1))
       
       
       #print(head(subdtdis$TF1))
       
       
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

        scale_fill_viridis_c(option = "a") +
        
        #geom_abline(slope = 1, intercept = 0) + 
        
        ggtitle("Concordant Harmony")
      
      ggplotly(hm) 
      
    })
     
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

        scale_fill_viridis_c(option = "e") +

        #geom_abline(slope = 1, intercept = 0) +
        
        labs(title = "Discordant Harmony")
        
      
     ggplotly(hm2)
      })
    
  output$conmotifTanglegram <- renderPlot(
    {
      motifhclust <- motifhclust()
      motifhclust$labels <- as.character(motifhclust$labels)
      dd <-
        dendlist(as.dendrogram(conhclustx()), as.dendrogram(motifhclust))
      tanglegram(dd, sub = "Concordant x Motifs", k_branches = 5, k_labels = 5, sort = T)
      }
    )
  
  output$phyloconTanglegram <- renderPlot(
    {
      my_alignment_sequence <- msa::msaConvert(htmsa, type="seqinr::alignment")
      distance_alignment <- dist.alignment(my_alignment_sequence)
      Tree <- bionj(distance_alignment)
      ph <- as.dendrogram.phylo(Tree)
      selids <- unique(subdt()[, TF1])
      
      subids <- idoptions[!(idoptions %in% selids)]

      #print(cat(unique(subids),sep = "|"))

      motifhclust <- motifhclust()
      motifhclust$labels <- as.character(motifhclust$labels)
      
      if (length(subids) > 0) {
        ph <- phylogram::prune(ph, pattern = subids)
      } 
      
      dd <- dendlist(as.dendrogram(conhclustx()), ph)
      
      tanglegram(dd, sub = "Concordant x Phylogeny", k_branches = 5, k_labels = 5, sort = T)
    }
  )
  output$phylodisTanglegram <- renderPlot(
    {
      my_alignment_sequence <- msa::msaConvert(htmsa, type="seqinr::alignment")
      distance_alignment <- dist.alignment(my_alignment_sequence)
      Tree <- bionj(distance_alignment)
      ph <- as.dendrogram.phylo(Tree)
      selids <- unique(subdt()[, TF1])
      
      subids <- idoptions[!(idoptions %in% selids)]
      
      #print(cat(unique(subids),sep = "|"))
      
      motifhclust <- motifhclust()
      motifhclust$labels <- as.character(motifhclust$labels)
      
      if (length(subids) > 0) {
        ph <- phylogram::prune(ph, pattern = subids)
      } 
      
      dd <- dendlist(as.dendrogram(dishclustx()), ph)
      
      tanglegram(dd, sub = "Discordant x Phylogeny", k_branches = 5, k_labels = 5, sort = T)
    }
  )
  
  output$dismotifTanglegram <- renderPlot({
    motifhclust <- motifhclust()
    motifhclust$labels <- as.character(motifhclust$labels)
    dd <-
      dendlist(as.dendrogram(dishclustx()), as.dendrogram(motifhclust))
    tanglegram(dd, sub = "Discordant x Motifs", k_branches = 5, k_labels = 5, sort = T)
  })
  
  output$condisTanglegram <- renderPlot({
    dd <-
      dendlist(as.dendrogram(conhclustx()), as.dendrogram(dishclustx()))
    tanglegram(dd, sub = "Concordant x Discordant", k_branches = 5, k_labels = 5, sort = T)
  })
  
  
  output$tanglegram <- renderPlot({
    choices <- c("Concordant", "Discordant", "Phylogeny", "Up Motif similarity", "Down Motif similarity")
    
    
    my_alignment_sequence <- msa::msaConvert(htmsa, type="seqinr::alignment")
    distance_alignment <- dist.alignment(my_alignment_sequence)
    Tree <- bionj(distance_alignment)
    ph <- as.dendrogram.phylo(Tree)
    selids <- unique(subdt()[, TF1])
    
    subids <- idoptions[!(idoptions %in% selids)]

    motifupclust <- motifupclust()
    motifupclust$labels <- sub("_up", "", motifupclust$labels)
    motifupclust$labels <- as.character(motifupclust$labels)
    
    motifdownclust <- motifdownclust()
    motifdownclust$labels <- sub("_down", "", motifdownclust$labels)
    motifdownclust$labels <- as.character(motifdownclust$labels)
    
    if (length(subids) > 0) {
      ph <- phylogram::prune(ph, pattern = subids)
    } 
    
    dd <- dendlist(as.dendrogram(conhclustx()), as.dendrogram(dishclustx()), ph, as.dendrogram(motifupclust), as.dendrogram(motifdownclust))
    
    tanglegram(dd, sub = paste(input$tanglechoice1, "x", input$tanglechoice2), k_branches = input$kbreaks, k_labels = input$kbreaks, sort = T, which = c(which(choices == input$tanglechoice1), which(choices == input$tanglechoice2)))
  })
  
  motifmat <- reactive({
    subids <- unique(subdt()[, TF1])
    # pattern <- paste0("^(", paste0(subids, collapse = "|"), ")_")
    # 
    # subpwms <- pwms[grepl(pattern, names(pwms))]
    
    #subpwms <- pwms[names(pwms) %in% subids]
    # subids is your vector of TF names like "ABF2", "ZAT7", etc.

    matches <- unlist(lapply(subids, function(tf) {
      which(startsWith(names(pwms), paste0(tf, "_")))
    }))
    
    subpwms <- pwms[unique(matches)]

    motifmat <- compare_motifs(subpwms, method = input$motifcompmethod)
  })
  
  motifhclust <- reactive({
    motifdist <- dist(motifmat(), method = input$distmeth)
    motifhclust <- hclust(motifdist, method = input$clustmeth)
  })
  
  motifupclust <- reactive({
    keep <- grep("_up", rownames(motifmat()), value = TRUE)
    
    submat <- motifmat()[keep, keep]
    motifdist <- dist(submat, method = input$distmeth)
    motifupclust <- hclust(motifdist, method = input$clustmeth)
  })
  
  motifdownclust <- reactive({
    keep <- grep("_down", rownames(motifmat()), value = TRUE)
    #message("subnames: ", keep)
    submat <- motifmat()[keep, keep]
    
    motifdist <- dist(submat, method = input$distmeth)
    motifdownclust <- hclust(motifdist, method = input$clustmeth)
  })
  
  output$motifheatmap <- renderPlotly({
    
    
    motiforder <- motifhclust()$labels[motifhclust()$order]
      
    motifdt <- as.data.table(motifmat(), keep.rownames = "TF1")
    motifnar <- melt.data.table(motifdt, id.vars = "TF1", variable.name = "TF2")
    
    motifnar$TF1 <- factor(motifnar$TF1, levels = motiforder)
    motifnar$TF2 <- factor(motifnar$TF2, levels = motiforder)
    #setkey(motifnar, TF1)
    
    ggplot(motifnar, aes(x = TF1, y = TF2, fill = value)) + geom_raster()
  })
  
  output$motifdend <- renderPlotly({
    dd <- as.dendrogram(motifhclust())
    ddplot <- ggdendrogram(dd)
    ggplotly(ddplot)
  })
  
  output$tfregdt <- renderDataTable({
    setDT(targetsubnar())
  })
  
  output$targetregdt <- renderDataTable({
    setDT(subvs())
  })
  
  targetsubnar <- reactive({
    padjcutoff <- as.numeric(input$regpcutoff)
    l2fccutoff <- as.numeric(input$reglfccutoff)
    
    subnar <- narpv[rn %in% input$allgenes][padj < padjcutoff][abs(log2FoldChange) > l2fccutoff]
  })
  
  subvs <- reactive({
    subnar <- targetsubnar()
    unique(vertices[(V1 %in% subnar$TF) | (V1 %in% subnar$rn) | (V1 %in% subnar$TF_ID), .(V1, label = TF, group = Family, value = TedStrength, shape = shape)])
  })
  
  output$regulators <- renderVisNetwork({
    
    subnar <- targetsubnar()
    #subnar[, c("ID", "TF") := tstrsplit(TF, "_")]
    
    # mygraph <- as_tbl_graph(
    #   subnar[,.(TF, rn, padj, value = log2FoldChange, color = Color)], directed = T, 
    #   )
    
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
    
    
    # ggi <- ggraph(mygraph, "tree") +
    # 
    #   geom_edge_diagonal(
    #     aes(
    #       color = factor(sign(log2FoldChange)),
    #       width = log2FoldChange
    #       ),
    #     arrow = arrow(
    #       #angle = ifelse(sign(newdt$TargetFC) == 1, 30, 90),
    #       type = "closed"
    #     ),
    # 
    #     start_cap = circle(5, 'mm'),
    #     end_cap = circle(5, 'mm'),
    #     #alpha = 0.6,
    #     strength = 0.2
    #   ) +
    # 
    #   geom_node_text(
    #     aes(label = name),
    #     size = 10,
    #     repel = T,
    #     #nudge_x = -0.1
    #   ) +
    # 
    #   geom_node_point(
    #     aes(size = Degree)
    #   ) +
    # 
    #   theme_void(base_size = 20)
    
    #newd3 <- igraph_to_networkD3(mygraph)

    # simpleNetwork(
    #   subnar[,.(TF, rn, padj, log2FoldChange)], 
    #   zoom = T, 
    #   linkColour = subnar$Color, 
    #   nodeColour = subnar$Color, 
    #   opacity = 0.8 
    #   )
      
     
  })

}

# Run the application 
shinyApp(ui = ui, server = server)
