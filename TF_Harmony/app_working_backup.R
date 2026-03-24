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

## Initial data read-in and setup ##

## spinner while graphs load
options(spinner.type = 5, spinner.color = "darkblue")

## setting up color scheme that is used in a couple of the network graphs - specifically the one that overlays
## the known activator/repressor data with TF activity
testcolors <- c("red", "gray", "blue", "gray", "darkred", "darkblue")
names(testcolors) <- c("Activator", "Minimally Active", "Repressor", "Unknown", "Upregulated", "Downregulated")

## This msa is used for the protein sequence dendrogram - needs to be redone to reflect phytozome
htmsa <- Biostrings::readAAMultipleAlignment("Data/msa.fa")

## These are the known transcription effector domains, used in combination with the color scheme for network graphs 
allteds <- fread("Data/allteds.tsv")


## Reading in the pre-calculated harmony - commented out is an older version, current version is 'no batch' and from the neb data **
#dt <- fread("harmonytable_justnobatch.tsv")
dt <- fread("Data/harmonytable_nobatch_nebs.tsv")

# Rename columns - only for nobatch_nebs.tsv
setnames(dt, c("TF", "Intersect_Concordant", "Intersect_Discordant", "PValue_Concordant", "PValue_Discordant", "Correlation_Concordant", "Correlation_Discordant", "Harmony_Concordant", "Harmony_Discordant"), c("TF1", "Concordant_Intersect", "Discordant_Intersect", "Concordant_PValue", "Discordant_PValue", "Concordant_Correlation", "Discordant_Correlation", "Concordant_Harmony", "Discordant_Harmony"))

## Harmony is calculated directionally - so each TF is 'TF1' and compared against all others, meaning I need to remove rows

###### NEW
## Trying out a signed Harmony score
## Convert harmony table to signed harmony
setDT(dt)
dt[, Harmony := Concordant_Harmony - Discordant_Harmony]

####
####

## TF1 is the same as TF2
dt <- dt[!(TF1 == TF2)]
## Cleaning up the TF names - these two steps can be done ahead of time to decrease intial load time
dt[, TF1 := sub("-.*", "", TF1)]
dt[, TF2 := sub("-.*", "", TF2)]

## Reading in the data table that associates TF IDs with TF families - which will then be merged later on
## can also be done ahead of time
tfswithfamilies <- fread("Data/tfidswithfamilies.tsv")

## Merging the harmony datatable with the family data - ie. giving TF1 and TF2 family ids
setkey(dt, TF1)
setkey(tfswithfamilies, Name)
dt <- dt[tfswithfamilies[, .(Name, TF1_Family = Family)]]
setkey(dt, TF2)
dt <- dt[tfswithfamilies[, .(Name, TF2_Family = Family)]]


## Reading in PWMs for motif sorting
pwms <- readRDS("Data/pwms.RDS")
#names(pwms) <- sub(".*?_", "", names(pwms))
## Converting pwms to a different format that one of the motif packages needed
pfms <- convert_motifs(pwms, class = "motifStack-pfm")

#narpv <- fread("nar_agg_no-batch_wTED.tsv")

## This should be all the DEGs for each tf in narrow format - no idea if its faster to read something with a lot of rows or a lot of columns, this could be sped up. 
narpv <- fread("Data/newnarpv_nebs_nobatch.tsv")
## adding a column for colorscheme for network graphs to reference
narpv[, Color := ifelse(sign(log2FoldChange) < 0, "Blue", "Red")]
## I swear this TF is named something different in every dataset. 
narpv[TF == "HSF.A4A", TF_ID := "AT4G18880"]

## Experiment 17 was eventually deemed to not pass QC, so removing those tf/targets. 
exp17nogo <- narpv[EXP == "1-17", unique(TF)]
narpv <- narpv[EXP != "1-17"]

## Creating a set of unique gene ids from the DEGs list, used for a selectable list later. 
allgeneids <- unique(narpv$rn)

## Again - pre-setup for network graphs, determining the overall set of vertices
gis <- narpv[, .(unique(rn))]
vertices <- merge(gis, allteds, by.x = "V1", by.y = "Locus", all.x = T)
vertices <- merge(vertices, unique(narpv[, .(TF_ID, TF)]), by.x = "V1", by.y = "TF_ID", all.x = T)
#vertices[!(is.na(TF)), V1 := TF]

## Reading in a dataset of gene expression by root cell type from the Benfey lab - I suspect this will be removed before publication
cte <- fread("Data/cellTypeExpression.tsv", key = "Gene ID")

## Similarly - reading in just-in-time dataset for roots
jitr <- fread("Data/JITGenes_root.csv", skip = 1, select = c("Gene", "FDR adjusted p-value", "First Response (Just-in-time bin)"), col.names = c("GeneID", "pvalue", "JIT"))
setkey(jitr, GeneID)

## Reading in just-in-time dataset for shoots
jits <- fread("Data/JITGenes_shoot.csv", skip = 1, select = c("AtID", "FDR adjusted p-value", "First Response (Just-in-time bin)"), col.names = c("GeneID", "pvalue", "JIT"))
setkey(jits, GeneID)

## Removing experiment 17 from harmony data
dt <- dt[!(TF1 %in% exp17nogo)]
dt <- dt[!(TF2 %in% exp17nogo)]

#idoptions <- dt[, unique(sub("AT.*?_", "", V1))]
#idoptions <- sub("_DES...*.csv", "", list.files("./DEGS", pattern = "^A.*"))


## ok so the orginal form of this app called DEGs individually instead of reading them all in at once - which is still being done in the 
## pairwise section I believe. So I have a folder with individual degs, and a large dataset with a narrow form of degs. This should be consolidated. Below is reading in the list individual TFs based on the files in the DEG folder. 
### modified when nobatch nebs were added because they were TSVs instead
idoptions <- sub("_DES...*.tsv", "", list.files("Data/DEGS", pattern = "^A.*"))
ids <- as.data.table(idoptions)
ids[, c("Ids", "TF") := tstrsplit(idoptions, "_")]
idoptions <- sub("^A.*?_", "", idoptions)
idoptions <- narpv[order(TF), unique(TF)]

## This is a version of the harmony table, but with fewer columns. 
## This is used to show Pairwise Harmony but just between the selected TFs as opposed to the raw data under Global Analyses
newdt <-
  dt[, .(
    TF1 = sub("AT.*?_", "", TF1),
    TF2 = sub("AT.*?_", "", TF2),
    Concordant = Concordant_Harmony,
    Discordant = Discordant_Harmony
  )][order(-Concordant)]

## Reading in a list of all Arabidopsis TF IDs - so that all TFs can be identified by a triangle in the network graphs
fulltfids <- fread("Data/all_ath_tf_ids.tsv")

## and more pre-setup for network graphs
vertices <- merge(vertices, fulltfids, by.x = "V1", by.y = "Gene_ID", all.x = T)
vertices[!(is.na(Family)), shape := "triangle"]
vertices[(is.na(Family)), shape := "circle"]
vertices[!(is.na(TF)), shape := "triangle"]
#vertices[!(is.na(TF)), V1 := TF]
vertices <- unique(vertices)


# Define UI for application that draws a histogram
ui <- navbarPage("Landscape of TF Harmony",

                 ## This is all the UI setup. Generally hiarchical starting with the top three panels: 
                 ## Global Analyses, Pairwise Analyses, and Target Regulation
                 ## Then each of those panels has sub panels and figures. 
                 
          tabPanel("Global Analyses",

                     inputPanel(
                       ## Defining cutoffs for harmony figures - user can set harmony minimums, pvalues, intersects etc
                       ## and the harmony table is subset based on these input
                       textInput(inputId = "familyharmonycutoff", label = "Harmony Minimum:", value = 0),
                       textInput(inputId = "familypvcutoff", label = "Fisher P-Value Maximum:", value = 0.05),
                       textInput(inputId = "familyintersectcutoff", label = "Intersect Minimum:", value = 0),
                       
                     selectInput(
                       ## User can define what type of distance calculation is used prior to clustering
                       inputId = "distmeth", 
                       label = "Distance Method", choices = c("euclidean", "maximum", "manhattan", "canberra"),
                       selected = "maximum", 
                       multiple = F
                     ),
                     selectInput(
                       ## User can define the clustering method used - this is all for the heatmaps
                       inputId = "clustmeth", 
                       label = "Clustering Method", 
                       choices = c("complete", "average", "ward.D", "single", "mcquitty", "median", "centroid", "ward.D2"),
                       selected = "ward.D2", 
                       multiple = F
                     ),
                     ## Similary user can define what type of algorithm is used for motif comparison
                     ## These are just options from the motif package
                     ## again this may be removed at a later point, I'm not sure how many people have strong opinions on 
                     ## motif clustering algorithms. 
                     selectInput(inputId = "motifcompmethod", label = "Motif Comparison Method:", choices = c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT", "HELL", "SEUCL", "MAN", "ALLR_LL", "WEUCL", "WPCC"), selected = "ALLR_LL", multiple = F, selectize = F)
                     
                   
                     
                   ),
                   selectInput(
                     ## the main input on the front page - asking the used which TF or Family to subset by
                     ## Starts pre-populated with all TF Families listed
                     inputId = "family1", 
                     label = "Search by TF or TF Family:", 
                     choices = unique(c(tfswithfamilies$Family, idoptions)),
                     selected = unique(tfswithfamilies$Family),
                     multiple = T,
                     width = "100%",  
                     selectize = T
                   ),
                   
                   tabsetPanel(
                     ## getting into the visuals below the user input options
                     ## Series of tabs under Global Analyses
                     tabPanel("Harmony Heatmap",
                              ## Main visual you see on load - concordant harmony heatmap, interactive plotly
                              tabsetPanel(
                                tabPanel("Concordant Harmony",
                              plotlyOutput(
                                outputId = "conHarmonyHeatmap",
                                
                                height = 2000
                              ) %>% withSpinner(), ),
                              
                              tabPanel("Discordant Harmony",
                              ## Alternative to main heatmap - show discordant harmony instead, interactive plotly
                              plotlyOutput(
                                outputId = "disHarmonyHeatmap",
                                
                                height = 2000
                              )%>% withSpinner(),) 
                              ),
                     ),
                     tabPanel("Family Harmony", 
                              ## Second tab under Global Analyses
                              ## Graphs look at prevalence of family interaction
                              ## Honestly these didn't see all that effective to me
                              ## They'll likely be removed
                              splitLayout(
                                plotlyOutput( outputId= "relationships", 
                                              height = 1000) %>% withSpinner(),
                                visNetworkOutput( outputId= "familyproportion",
                                                  height = 1000) %>% withSpinner()
                                
                              ),
                     ),
                     tabPanel("Tanglegram",
                              ## Third tab under Global Analyses
                              ## allows the user to select their inputs and compare dendrograms - highlighting similarities
                              ## Concordant / Discordant are harmony dendrograms
                              ## Phylogeny is a pre-computed MSA that is trimmed based on used selections
                              ## Motif similarity is calculated on-the-fly I believe based on pre-computed PWMs which are 
                              ## filtered based on user input (ie. only bZIPs)
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
                                  ## user input that allows to you highlight breaks
                                  inputId = "kbreaks",
                                  label = "Number of breaks:",
                                  choices = c(0,1,2,3,4,5,6,7,8,9),
                                  selected = 1, 
                                  multiple = F
                                )
                              ), 
                              ## tanglegram plot output based on above selections
                              plotOutput(outputId = "tanglegram", height = 1000) %>% withSpinner()

                              

                     ),
                     tabPanel("Harmony Table", 
                              ## raw harmony table for user to peruse and interact with
                              ## still under global analyses
                                DT::dataTableOutput(outputId = "fullharmonydt")
                              ),
                     
                     tabPanel("Motif Heatmap",
                            ## heatmap of motif similarity
                            ## changes based on user selection of motif comparison algorithm
                              plotlyOutput(
                                outputId = "motifheatmap", height = 1000
                                )%>% withSpinner()
                              ),
                     tabPanel("Motif Dendrogram",
                              ## shows the *actual* motifs allowing you to actually see how similar they are
                              ## downside is that it takes forever to load if you try to show All the motifs at once
                              ## should probably be limited to like 5 tfs at a time or so
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
                   ## Done with Global analyses, this is the middle panel up at the top
                   ## Defines all the TF x TF comparisons in more detail
                     inputPanel(
                       ## User input for selecting which tfs they want to compare
                       ## Based on which DEG files are available to read in individually
                       selectizeInput(
                         inputId = "TFs",
                         label = "TFs",
                         choices = idoptions,
                         multiple = T,
                         selected = idoptions[1:2],
                         options = list(maxItems = 10)
                       ), 
                       ## The idea here is that you can read in a list of TFs from your own project and see their defined harmony
                       fileInput(inputId = "tfupload", label = "Upload list of TFs"),

                       ## Users can adjust size of pairwise harmony plot
                       sliderInput(
                         inputId = "harmonyplotheight",
                         min = 500,
                         max = 10000,
                         step = 100,
                         value = 1000,
                         label = "Plot height:"

                       ),
                       ## users can decide if they want separate pairwise plots for each combination of TF (faceted) or 
                       ## one plot with all the pairwise combinations shown in different colors
                       checkboxInput(inputId = "harmonyinterfacet", label = "Facet by interaction:", value = T),
                     ),
                   tabsetPanel(id = "pairwise",
                   tabPanel("TF x TF",
                   ## First pairwise plot
                   ## Shows all shared DEGs, direction of regulation, and correlation of Log2FC, computed on the fly
                   ## On-the-fly was picked so that users can upload their own DEG datasets and compare to ours
                   splitLayout(
                        plotlyOutput(outputId = "harmonyplotly", height = "100%") %>% withSpinner(),
                        DT::dataTableOutput("harmonyTable")
                   )

                   
          ),
          ## More motif comparisons, but this time pairwise
          ## Idea being that you can see if a lot of shared regulation is associated with a similar motif
          tabPanel("Motif Comparison", value = "pairwisemotifs",
                   selectInput("motiftreestyle", label = "Tree Style:", choices = c("stack", "tree", "radialPhylog"), multiple = F),
                   #uiOutput(outputId = "tfmotif"),
                   imageOutput(
                     outputId = "motifplot", 
                     width = "80%", 
                     height = 800
                   )%>% withSpinner()
                   
                   
          ),
          ## Are the shared DEGs expressed in a specific root cell type at baseline? 
          ## Overlay of Benfey cell type data subsetted to match shared DEGs
          tabPanel("Cell Type Specificity",
                   plotlyOutput(
                     outputId = "ctplot",  
                     height = 1000
                   )
          ),
          ## Are the shared DEGs differentially expressed at a given time period in the NxTime datasets? 
          ## Overlay of NxTime data - subsetted to match shared DEGs
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
          ## On the fly calculation of GO Terms using the GOST server
          ## Subsetted DEGs are submitted to the server at time of selection - runs even when page isn't selected
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
          
          ## Do the TFs also regulate one another, or do they just share DEGs
          ## I can't remember if the color scheme was correct on this one, might need updating
          ## Users can adjust the style of the network graph shown
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
                   ## Are we subsetting to on shared DEGs? ** TODO - check and see what this subsets? 
                   checkboxInput(
                     inputId = "networkdegs", 
                     label = "DEGs only:", 
                     value = TRUE
                     ), 
                   ),
            ## vis Network plot
            visNetworkOutput(
              outputId = "networkplot", 
              height = 1000
              )%>% withSpinner()
          ), 
          
          ## Separate panel showing overlapping DEGs between TFs - harmony excluded
          ## Idea here was whether we could see if we could identify the primary TF in cascading networks
          ## Degree cutoff is to limit to targets that are downstream of a minimum number of TFs
          ## ie. only keep a target is all 3 of my selected TFs regulate it
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
          ## Final overarching panel - Target Regulation
          ## This is the reverse of the prevous two panels. 
          ## Instead of looking at which targets a TF regulates we start with the targets and look at all TFs that regulate it. 
          ## You can also upload a list of targets (ie. DEGs from your own experiment) and see all the TFs
          ## that regulate those targets. 
          tabPanel("Target Regulation",
                   inputPanel(
                     fileInput(inputId = "targetupload", label = "Upload list of Targets"),
                     textInput(inputId = "regpcutoff", label = "Adjusted P-value Cutoff:", value = "0.05"),
                     textInput(inputId = "reglfccutoff", label = "Log2FoldChange Cutoff:", value = "0"),
                     textInput(inputId = "networkhcutoff", label = "Harmony Cutoff:", value = "0.5"),
                     
                     
                     
                     
                     # selectInput(
                     #   inputId = "targetnetworkstyle",
                     #   label = "Network Style:",
                     #   choices = c(
                     #     "layout_with_sugiyama",
                     #     "layout_on_sphere",
                     #     "layout_with_fr",
                     #     "layout_with_gem",
                     #     "layout_with_kk",
                     #     "layout_with_lgl",
                     #     "layout_with_mds",
                     #     "layout_in_circle",
                     #     "layout_nicely",
                     #     "layout_as_star",
                     #     "layout_as_tree",
                     #     "layout_on_grid",
                     #     "layout_with_dh"
                     #   ),
                     #   selected = "layout_with_sugiyama",
                     #   multiple = F
                     # )
                   ),
                   
                   ## This is where the list of all targets from NarPV is used
                   selectInput(
                     inputId = "allgenes",
                     label = "Targets: ",
                     choices = allgeneids,
                     multiple = T, 
                     width = "100%",
                     selected = allgeneids[1]
                   ),
                   
                   
                            
                            fluidRow(
                              column(12,
                                     # your inputs here
                              )
                            ),
                            
                            fluidRow(
                              column(12,
                                     plotOutput("tfHarmonyScatter", height = "350px")
                              )
                            ),
                   fluidRow(
                     column(2,
                            selectInput(
                              "module_select",
                              "Select module",
                              choices = "All",
                              selected = "All"
                            )
                            ),
                     column(2,
                   sliderInput(inputId = "harmtrans", label = "Harmony Transparency:", min = 0, max = 1, value = 0.25, step = 0.05, ticks = FALSE)),
                   column(4, 
                   sliderInput(inputId = "logtrans", label = "Log2FC Transparency:", min = 0, max = 1, value = 1, step = 0.05, ticks = FALSE)
                   ),
                   column(2,
                          selectInput(
                            "louvain_weight_mode",
                            "Module definition",
                            choices = c(
                              "Absolute Harmony" = "abs_harmony",
                              "Absolute Harmony squared" = "abs_harmony_sq",
                              "Concordant Harmony" = "concordant",
                              "Discordant Harmony" = "discordant"
                            ),
                            selected = "abs_harmony"
                          ))
                     
                   ),
                     
                     
                            fluidRow(
                              column(12,
                                     visNetworkOutput("regulators", height = "850px")
                              )
                            ),
                   fluidRow(
                     plotOutput("targetRegHeatmap", height = "700px")
                   )
                   
                            
                            # fluidRow(
                            #   column(6,
                            #          dataTableOutput("tfregdt")
                            #   ),
                            #   column(6,
                            #          dataTableOutput("targetregdt")
                            #   )
                            # )
                   ,

                   
                   
          ## I had to include this after adding some packages - never figured out which. 
          ## I added the spinner, and motif, and tanglegram packages all around the same time
          ## Then when I got everything up and working all the basic formatting had changed
          ## I tried experimenting and removing individual packages at a time but couldn't determine the culprit
          ## So did a work around where I just redefined the basics 
          includeCSS("~/Desktop/github/HTT_144TFs/TF_Harmony/www/style.css")

          )
          
        )
        
## Defining a handful of functions between the frontend and backend sections. 

## This is the Fisher Test function. You provide the number of shared targets, and then two 
## sets of targets, it sets up the contingency table and runs the fisher exact test.
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


## Function to return a linked list based on a given pattern (pat)
## Subsets NarPV - the DEGs data based on the pattern
## Only grabs rows where the padj is less than 0.05 - grabs the AtID (rn), Log2FC, and padj. 
## Names this list based on the pattern
## For loop cycles through a series of TF ids - so all targets for multiple tfs will be returned in a named list
## if old TF ids are still present (oldll) it doesn't recalculate them
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

## Given multiple differential-expression result tables, compare every pair by merging on gene id, then tag each gene as 
## concordant/discordant in direction of effect, caching previous pair computations if provided.

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


## I played with several different ways of displaying the pairwise motifs
## Can't remember if this is still used... 
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

# Define server logic 
## Backend computation
  server <- function(input, output, session) {
  
  ## DEG datasets for TF1 & TF2 from pairwise analysis
  ## This was the original design when I only allowed two tfs... 
  #Tar1 <- reactive({findfile(input$TF1)})
  #Tar2 <- reactive({findfile(input$TF2)})
  
  ## Above should have been replaced with the ability to call multiple tfs. So now its a selectize input limited to TF ids
  ## findfile was refactored to use the for loop and cycle through for as many tfs as needed. 
  ## Uses NarPV to calculate overlapping DEGs on the fly.   
  #print(input$TFs)
  tarfs <- reactive({
    if (exists("tarfs()")) findfile(input$TFs, tarfs())
    else findfile(input$TFs)
    })
  matched <- reactive({
    if (exists("matched()")) matchtidy(tarfs(), matched())
    else matchtidy(tarfs())
    })
  
  ## Restricts Harmony to only use the input TFs
  ## Uses newdt - the harmony data with fewer columns 
  subgraph <- reactive({
    #newdt[(TF1 == input$TF1 | TF2 == input$TF2 | TF1 == input$TF2 | TF2 == input$TF1) & !(is.na(Concordant) & is.na(Discordant))]
    newdt[(TF1 %in% input$TFs)][TF2 %in% input$TFs]#[!(is.na(Concordant) & is.na(Discordant))]
    })
  
  
  ## The total number of shared DEGs 
  shared <- nrow(matched())
  
  # concordantshared <- reactive({nrow(matched()[Harmony == "Concordant"])})
  # discordantshared <- reactive({nrow(matched()[Harmony == "Discordant"])})
  # 
  # concorrelation <- reactive({cor(matched()[Harmony == "Concordant", log2FC.x], matched()[Harmony == "Concordant", log2FC.y])})
  # discorrelation <- reactive({cor(matched()[Harmony == "Discordant", log2FC.x], matched()[Harmony == "Discordant", log2FC.y])})
  # 
  # conftest <- reactive({ftestfun(concordantshared(), Tar1(), Tar2())})
  # disftest <- reactive({ftestfun(discordantshared(), Tar1(), Tar2())})
  
  
  ## So sometimes the motifs don't show, added this to try to force a 'redraw' 
  ## every time you click the pairwise tab
  trigger_motif_redraw <- reactiveVal(0)
  
  observeEvent(input$tabs, {
    if (input$tabs == "pairwise") {
      trigger_motif_redraw(trigger_motif_redraw() + 1)
    }
  })

  
  ## This is the motif plot from Pairwise Analysis
  ## Subsets the pfms dataset and then grabs the pre-existing PNG file and displays it 
  ## Users can select what type of visual - tree vs stack vs no dendrogram
  ## So it compares motifs based on pfms and then renders images from pre-computed motifs
  
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
  
  ## Another motif browser - but this one is under Global Harmony
  ## Instead of loading pre-computed images it subsets the pfms based on user input
  ## and then generates a motif on the fly
  ## In the latest version the tab panel has been commented out - so this shouldn't be called
  
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
  
  
  ## Showing subset harmony table in an interactive DataTable with scroll bars
  ## This is under *pairwise harmony*
  ## Only shows harmony for subsetted TF list based on user input
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
    
    ## Reactive function to subset harmony data on Global Analyses based on user input
    ## Most future Global Analysis harmony references use subdt**
    ## Subsetted in such a way that you can specify TF or TF Family
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
    
    ## Showing subset harmony table in an interactive DataTable with scroll bars
    ## This is under Global Analyses harmony
    ## Only shows harmony for subsetted TF list based on user input
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
    })
    
    
    ## Rendering plot output of GO terms
    output$goplot <- renderPlotly({
      
      ggplotly(gostplot(goout(), interactive = F), height = input$goplotheight)
      
    })
    
    ## Rendering GO term table from output
    ## Needs fixing - table does not scroll and long tables extend beyond the screen
    ## Likely - gosttable will need to be parsed independently and repackaged into a DataTable
    output$gotable <- renderPlot({
      #print(head(goout()))
      publish_gosttable(goout(), ggplot = T)
    })
    
    ## Telling the app to submit the GO terms in the background to decrease lag time
    outputOptions(output, "goplot", suspendWhenHidden = FALSE)
    
    
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
    
    ## User option to upload a list of targets
    ## Target regulation
    ## Updates selectInput to match upload
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
    
    
    ## I don't remember if this is still used - it was needed because the harmony table would change when the user selected certain rows
    ## I think I removed that feature
    proxy <- dataTableProxy(outputId = "harmonyTable")
    
    ## For the harmony cutoff user input
    ## The slider is updated automatically to use the table minimum and maximum as the extreme values
    observeEvent(input$TFs, {
      meltdt <- melt.data.table(subgraph())
      meltdt <- meltdt[!(is.na(value) | is.infinite(value))]
      proxy %>% selectRows(newdt[TF1 %in% input$TFs & TF2 %in% input$TFs, which = T],)
      updateSliderInput(session, inputId = "HarmonyRange", value = 0, min = 0, max = max(abs(meltdt$value)))
    })
    
    ## Same thing as above but for TF2 now... actually now that I think of it 
    ## this probably isn't in use anymore either when I switched over to the multi-tf input... 
    ## going to try commenting it out to see what happens
    # observeEvent(input$TF2, {
    #   meltdt <- melt.data.table(subgraph())
    #   meltdt <- meltdt[!(is.na(value) | is.infinite(value))]
    #   proxy %>% selectRows(newdt[TF1 == input$TF1 & TF2 == input$TF2, which = T],)
    #   updateSliderInput(session, inputId = "HarmonyRange", value = 0, min = 0, max = max(meltdt$value))
    # })
    
    
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
      #mygraph <- as_tbl_graph(subnar[,.(TF, rn, padj, value = log2FoldChange, group)], directed = T)
      
      #print(head(subnar))
      
      ## Creating the graph data structure from the DEG subset
      mygraph <- graph_from_data_frame(
        d = subnar[, .(TF_ID, rn, label = TF, value = log2FoldChange, padj, group)], 
        vertices = unique(vertices[(V1 %in% subnar$TF) | (V1 %in% subnar$rn) | (V1 %in% subnar$TF_ID), .(V1, label = TF, group = Activity, value = TedStrength, shape = shape)])
        )
      
      
      ## Calculating network degree
      V(mygraph)$Degree <- igraph::degree(mygraph)
      ## Pruning the network to only show vertices above the user cutoff for degree minimum
      mygraph <- delete.vertices(mygraph, V(mygraph)$Degree < degreecutoff)
      
      #V(mygraph)$group <- ifelse(names(V(mygraph)) %in% narpv$TF, narpv[TF %in% names(V(mygraph)),  Activity], "Unknown")
      #V(mygraph)$value <- ifelse(names(V(mygraph)) %in% narpv$TF, narpv[TF %in% names(V(mygraph)),  TedStrength], 0)
      
      ## Setting the colors based on the color scheme established in the setup sectin - before UI
      V(mygraph)$color <- testcolors[V(mygraph)$group]
      E(mygraph)$color <- testcolors[E(mygraph)$group]
      
      #print(head(mygraph))
      
      ## Showing the network graph and defining a handful of options
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
    
    ## This is a different network graph showing harmony - TF Regulation tab under Pairwise Analyses
    ## Generates a network graph from subgraph - which is the subsetted harmony datatable
    ## The graph itself shows concordant vs. discordant harmony by colored arrows assuming there is significant harmony between two tfs
    
    output$networkplot <- renderVisNetwork({
      ## Melting datatable, subsetting, adding color labels
      melted <- melt.data.table(subgraph(), variable.name = "group")
      melted <- melted[!(is.na(value) | is.infinite(value)) & abs(value) > input$HarmonyRange ]
      melted[, color := ifelse(group == "Concordant", "red", "blue")]
      
      
      ## Checks to see if DEGs is checked
      ## Grabs all the TF2s -> TFs that have harmony with TF1 from the melted datatable
      ## Goes through all the TFs in the input and double filters - must be a TF and must be an AtID in tarffs
      ## Then subsets melted to only be the tfs in tf2ids
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
      
      
      ## Creating graph data structure, ensuring it's directional to match harmony 
      mygraph <- as_tbl_graph(as.data.table(melted), directed = T)
      ## calculating network degree
      V(mygraph)$Degree <- igraph::degree(mygraph)
      
      ## Displaying the network graph
      visIgraph(mygraph, layout = input$tfregnetworkstyle) %>% 
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), collapse = T, width = "100%") %>%
        #visEdges(arrows = "to" ) %>%
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
        
        geom_smooth(aes(group = interaction(Inter, Harmony)), method = "lm", se = F, size = 1) +
        
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
      
      
      
    })
    
    ## Reactive object to manually calculate harmony clustering 
    ## This is the concordant clustering
    ## Used a few times under Global Analyses tab
    ## Subsets harmony data
    ## Converts to wide format, then changes to a matrix
    ## Calculates distances - based on user algorithm
    ## Then clusters - based on user algorithm
    ## This is done separately for both X and Y axes as harmony is directional
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
    
    
    ## Reactive object to calculate clustering order for Y axis concordant harmony
    ## Global Analyses
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
    
    
    ## Reactive object to calculate discordant harmony clustering for the x axis
    ## And yes - clustering should be refactored to it's own function to speed things up
    ## Global Analyses
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
    
    
    ## reactive object to calculate discordant harmony clustering for the y axis
    ## Global Analyses
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

        scale_fill_viridis_c(option = "A") +
        
        #geom_abline(slope = 1, intercept = 0) + 
        
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

        scale_fill_viridis_c(option = "e") +

        #geom_abline(slope = 1, intercept = 0) +
        
        labs(title = "Discordant Harmony")
        
      
     ggplotly(hm2)
      })
    
    
  ## Creating the tanglegram for concordant harmony vs. motifs
    ## Uses cluster object, and then creates dendrogram for each
     
  output$conmotifTanglegram <- renderPlot(
    {
      motifhclust <- motifhclust()
      motifhclust$labels <- as.character(motifhclust$labels)
      dd <-
        dendlist(as.dendrogram(conhclustx()), as.dendrogram(motifhclust))
      tanglegram(dd, sub = "Concordant x Motifs", k_branches = 5, k_labels = 5, sort = T)
      }
    )
  
  ## Creating tanglegram for protein sequence alignment vs. concordant harmony
  ## uses precomputed msa 
  ## calculates dendrogram for msa 
  ## prunes dendrogram based on user inputs
  
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
  
  ## Creating tanglegram for protein sequence vs discordant phylogeny
  ## Not sure why motifhclust is called here.... need to look into this
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
  
  ## Creating tanglegram for motifs vs. discordant harmony
  output$dismotifTanglegram <- renderPlot({
    motifhclust <- motifhclust()
    motifhclust$labels <- as.character(motifhclust$labels)
    dd <-
      dendlist(as.dendrogram(dishclustx()), as.dendrogram(motifhclust))
    tanglegram(dd, sub = "Discordant x Motifs", k_branches = 5, k_labels = 5, sort = T)
  })
  
  ## Creating tanglegram for Concordant vs. Discordant Harmony
  ## Just creates two dendrograms based on clustered objects
  output$condisTanglegram <- renderPlot({
    dd <-
      dendlist(as.dendrogram(conhclustx()), as.dendrogram(dishclustx()))
    tanglegram(dd, sub = "Concordant x Discordant", k_branches = 5, k_labels = 5, sort = T)
  })
  
  ## I think this is the updated version of the above code. Instead of specific comparisons pre-computed I switched to a choice input
  ## calculating all the different dendrograms necessary and then changing up the graph based on the input choice
  ## I will try commenting out the previous later. 
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
  
  ## This creates the matrix of motif similarity 
  ## Subsets the pwms object based on user input then calls compare_motifs
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
  
  ## Calculates clustering object for motifs
  motifhclust <- reactive({
    motifdist <- dist(motifmat(), method = input$distmeth)
    motifhclust <- hclust(motifdist, method = input$clustmeth)
  })
  
  ## Specifically calculating the motif similarity but for motifs calculated on only UP Regulated targets
  ## I don't remember where this is used
  motifupclust <- reactive({
    keep <- grep("_up", rownames(motifmat()), value = TRUE)
    
    submat <- motifmat()[keep, keep]
    motifdist <- dist(submat, method = input$distmeth)
    motifupclust <- hclust(motifdist, method = input$clustmeth)
  })
  
  ## Calculating clustered object but for motifs calculated only on DOWN Regulated targets
  motifdownclust <- reactive({
    keep <- grep("_down", rownames(motifmat()), value = TRUE)
    #message("subnames: ", keep)
    submat <- motifmat()[keep, keep]
    
    motifdist <- dist(submat, method = input$distmeth)
    motifdownclust <- hclust(motifdist, method = input$clustmeth)
  })
  
  
  ## Heatmap of motif similarity from precomputed clustered objects
  
  output$motifheatmap <- renderPlotly({
    
    
    motiforder <- motifhclust()$labels[motifhclust()$order]
      
    motifdt <- as.data.table(motifmat(), keep.rownames = "TF1")
    motifnar <- melt.data.table(motifdt, id.vars = "TF1", variable.name = "TF2")
    
    motifnar$TF1 <- factor(motifnar$TF1, levels = motiforder)
    motifnar$TF2 <- factor(motifnar$TF2, levels = motiforder)
    #setkey(motifnar, TF1)
    
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
  output$tfregdt <- renderDataTable({
    setDT(targetsubnar())
  })
  
  ## This is a datatable in Target Regulation - third tab
  ## Takes the target subsetted/filtered data from DEGs data
  ## and then filters the vertices data structure based on DEG target regulation
  ## Returns the datatable showing TF, Family, Transcriptional Enhancer Strength, and shape - is a tf or not
  output$targetregdt <- renderDataTable({
    setDT(subvs())
  })
  
  ## reactive function to create targetsubnar datatable
  ## uses User supplied cutoffs to filter and subset DEG data
  ## Restricts based on a list of Targets
  targetsubnar <- reactive({
    padjcutoff <- as.numeric(input$regpcutoff)
    l2fccutoff <- as.numeric(input$reglfccutoff)
    
    narpv[rn %in% input$allgenes][padj < padjcutoff][abs(log2FoldChange) > l2fccutoff]
  })
  
  ###################
  ### Newly added ###
######################
  
  tf_lookup <- reactive({
    unique(narpv[, .(TF_ID, TF)])
  })
  
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
      .(
        TF1_ID,
        TF2_ID,
        abs_harmony,
        abs_harmony_sq,
        concordant_weight,
        discordant_weight
      )
    ]
    
    mods <- compute_louvain_modules(h, weight_col = weight_col, cutoff = 0.02)
    
    if (is.null(mods)) {
      return(data.table(
        TF_ID = hit,
        module = seq_along(hit)
      ))
    }
    
    mods[]
  })
  
    # Pick a single TF key to join on everywhere (TF_ID preferred; fall back to TF)
  tf_key_col <- reactive({
    if ("TF_ID" %in% names(targetsubnar())) "TF_ID" else "TF"
  })
  
  # A per-TF summary for the uploaded gene set
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
    
    # clean NaN from empty subsets
    clean_cols <- c("mean_up_l2fc", "mean_down_l2fc")
    for (j in clean_cols) {
      out[is.nan(get(j)), (j) := 0]
    }
    
    out[]
  })

  
  # Mean Harmony of each "hit TF" to the other "hit TFs"
  tf_mean_harmony <- reactive({
    reg <- tf_reg_summary()
    hit <- unique(reg$TF_ID)
    req(length(hit) >= 1)
    
    hdt <- target_harmony_dt()
    
    h <- hdt[
      TF1_ID %in% hit & TF2_ID %in% hit & TF1_ID != TF2_ID,
      .(TF1_ID, TF2_ID, Concordant_Harmony, Discordant_Harmony, Harmony)
    ]
    
    # outgoing from TF1 perspective
    out_sum <- h[, .(
      sum_concordant_harmony_out = sum(Concordant_Harmony[!is.na(Concordant_Harmony)], na.rm = TRUE),
      sum_discordant_harmony_out = sum(Discordant_Harmony[!is.na(Discordant_Harmony)], na.rm = TRUE),
      mean_concordant_harmony_out = mean(Concordant_Harmony[!is.na(Concordant_Harmony) & Concordant_Harmony != 0], na.rm = TRUE),
      mean_discordant_harmony_out = mean(Discordant_Harmony[!is.na(Discordant_Harmony) & Discordant_Harmony != 0], na.rm = TRUE),
      n_harmony_connections_out = uniqueN(TF2_ID[!is.na(Harmony) & Harmony != 0])
    ), by = .(TF_ID = TF1_ID)]
    
    # incoming to TF2 perspective
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
  
  tf_scatter_dt <- reactive({
    reg  <- tf_reg_summary()
    mh   <- tf_mean_harmony()
    mods <- tf_modules()
    
    out <- merge(reg, mh, by = "TF_ID", all.x = TRUE)
    out <- merge(out, mods, by = "TF_ID", all.x = TRUE)
    
    fill_zero_cols <- c(
      "sum_concordant_harmony_out",
      "sum_discordant_harmony_out",
      "mean_concordant_harmony_out",
      "mean_discordant_harmony_out",
      "n_harmony_connections_out",
      "sum_concordant_harmony_in",
      "sum_discordant_harmony_in",
      "mean_concordant_harmony_in",
      "mean_discordant_harmony_in",
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
    
  heatmap_dt <- reactive({
    subnar <- targetsubnar()
    reg    <- tf_scatter_dt()
    
    d <- merge(
      subnar,
      reg[, .(TF_ID, label, module)],
      by = "TF_ID",
      all.x = TRUE
    )
    
    d[]
  })
  
  output$targetRegHeatmap <- renderPlot({
    d <- heatmap_dt()
    req(nrow(d) > 0)
    
    # build matrix
    mat_dt <- d[, .(value = mean(log2FoldChange, na.rm = TRUE)),
                by = .(label, rn, module)]
    
    mat_wide <- data.table::dcast(
      mat_dt,
      label + module ~ rn,
      value.var = "value",
      fill = NA
    )
    
    row_meta <- mat_wide[, .(label, module)]
    mat <- as.matrix(mat_wide[, !c("label", "module")])
    rownames(mat) <- mat_wide$label
    
    # clustering matrix (NA → 0 only for clustering)
    mat_clust <- mat
    mat_clust[is.na(mat_clust)] <- 0
    
    # cluster targets using correlation distance
    if (ncol(mat_clust) > 1) {
      col_cor <- cor(mat_clust, use = "pairwise.complete.obs")
      col_dist <- as.dist(1 - col_cor)
      col_hc <- hclust(col_dist, method = "complete")
      gene_order <- colnames(mat_clust)[col_hc$order]
    } else {
      gene_order <- colnames(mat_clust)
    }
    
    # long format for plotting
    plot_dt <- data.table::melt(
      data.table(label = rownames(mat), mat, check.names = FALSE),
      id.vars = "label",
      variable.name = "rn",
      value.name = "log2FoldChange"
    )
    
    plot_dt <- merge(plot_dt, row_meta, by = "label", all.x = TRUE)
    
    # order TFs within module
    plot_dt[, label := factor(label, levels = unique(plot_dt[order(module, label)]$label))]
    plot_dt[, rn := factor(rn, levels = gene_order)]
    plot_dt[, reg_dir := sign(log2FoldChange)]
    
    ggplot(plot_dt, aes(x = rn, y = label, fill = factor(reg_dir))) +
      geom_tile(color = NA) +
      
      # facet by module
      facet_grid(module ~ ., scales = "free_y", space = "free_y") +
      
      scale_fill_manual(
        values = c(
          "-1" = "blue",
          "0" = "white",
          "1" = "red"
        ),
        name = "Regulation"
      )+
      
      labs(
        x = "Target genes (clustered)",
        y = "TFs (grouped by module)",
        fill = "log2FC"
      ) +
      
      theme_bw() +
      theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })
  
  
  #####################
  ## New scatter plot
  #######################
  
  output$tfHarmonyScatter <- renderPlot({
    d <- tf_scatter_dt()
    req(nrow(d) > 0)
    
    d[, harmony_y := sum_concordant_harmony_out - sum_discordant_harmony_out]
    
    ggplot(d, aes(x = n_targets, y = harmony_y, color = factor(module))) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(size = n_targets), alpha = 0.85) +
      ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
      labs(
        x = "# regulated target genes",
        y = "Outgoing Harmony balance",
        color = "Module",
        size = "# targets"
      ) +
      theme_bw()
  })
  
  
  
  ###################
  ###################
  ###################
  
  
  ## Reactive function that subsets and filters vertice datastructure
  ## Restricts vertices to only the TFs regulating a given set of targets
  ## Relabels resulting output columns for readability
  subvs <- reactive({
    subnar <- targetsubnar()
    unique(vertices[(V1 %in% subnar$TF) | (V1 %in% subnar$rn) | (V1 %in% subnar$TF_ID), .(V1, label = TF, group = Family, value = TedStrength, shape = shape)])
  })
  
    
  
  ##########################
  ## Module Detection ######
  ##########################
  
  compute_louvain_modules <- function(h, weight_col = "weight", cutoff = 0.02) {
    
    h <- copy(h)
    h <- h[!is.na(get(weight_col)) & get(weight_col) >= cutoff]
    
    if (nrow(h) == 0) {
      return(NULL)
    }
    
    h[, pair_id := ifelse(
      TF1_ID < TF2_ID,
      paste(TF1_ID, TF2_ID, sep = "__"),
      paste(TF2_ID, TF1_ID, sep = "__")
    )]
    
    h <- h[, .(
      TF1_ID = first(TF1_ID),
      TF2_ID = first(TF2_ID),
      weight = max(get(weight_col), na.rm = TRUE)
    ), by = pair_id]
    
    g_mod <- igraph::graph_from_data_frame(
      d = h[, .(from = TF1_ID, to = TF2_ID, weight)],
      directed = FALSE
    )
    
    cl <- igraph::cluster_louvain(g_mod, weights = igraph::E(g_mod)$weight)
    memb <- igraph::membership(cl)
    
    data.table(
      TF_ID = names(memb),
      module = as.integer(memb)
    )
  }
  

  
  tf_modules <- reactive({
    reg <- tf_reg_summary()
    hit <- unique(reg$TF_ID)
    req(length(hit) >= 1)
    
    hdt <- target_harmony_dt()
    harmony_cutoff <- 0.05
    
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
      .(
        TF1_ID,
        TF2_ID,
        abs_harmony,
        abs_harmony_sq,
        concordant_weight,
        discordant_weight
      )
    ]
    
    # keep only rows with meaningful weight in the chosen mode
    h <- h[!is.na(get(weight_col)) & get(weight_col) >= harmony_cutoff]
    
    if (nrow(h) == 0) {
      return(data.table(
        TF_ID = hit,
        module = seq_along(hit)
      ))
    }
    
    h[, pair_id := ifelse(
      TF1_ID < TF2_ID,
      paste(TF1_ID, TF2_ID, sep = "__"),
      paste(TF2_ID, TF1_ID, sep = "__")
    )]
    
    h <- h[, .(
      TF1_ID = first(TF1_ID),
      TF2_ID = first(TF2_ID),
      weight = max(get(weight_col), na.rm = TRUE)
    ), by = pair_id]
    
    g_mod <- igraph::graph_from_data_frame(
      d = h[, .(from = TF1_ID, to = TF2_ID, weight)],
      directed = FALSE,
      vertices = data.frame(name = hit)
    )
    
    cl <- igraph::cluster_louvain(g_mod, weights = igraph::E(g_mod)$weight)
    memb <- igraph::membership(cl)
    
    data.table(
      TF_ID = names(memb),
      module = as.integer(memb)
    )
  })
  
  observe({
    d <- tf_scatter_dt()
    req(nrow(d) > 0)
    
    mods <- sort(unique(d$module))
    
    updateSelectInput(
      session,
      "module_select",
      choices = c("All", as.character(mods)),
      selected = "All"
    )
  })
  
  
    

  
  #############################
  ## New regulators gene network
  ##############################
  
  output$regulators <- renderVisNetwork({
    
    reg_all <- tf_scatter_dt()

    req(nrow(reg_all) > 0)
    
    # Ensure stable label column exists
    reg_all[, label := fifelse(is.na(TF) | TF == "", TF_ID, TF)]
    
    # Optional module filter
    reg <- copy(reg_all)
    if (!is.null(input$module_select) && input$module_select != "All") {
      reg <- reg[module == as.integer(input$module_select)]
    }
    
    req(nrow(reg) > 0)
    
    hit <- unique(reg$TF_ID)
    
    subnar <- copy(targetsubnar())
    subnar <- subnar[TF_ID %in% hit]
    req(nrow(subnar) > 0)
    
    # -------------------------
    # Nodes
    # -------------------------
    module_colors <- c(
      "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
      "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
      "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"
    )
    
    tf_nodes <- reg[, .(
      id = TF_ID,
      label = label,
      group = "TF",
      shape = "triangle",
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
    
    gene_nodes <- data.table(
      id = unique(subnar$rn),
      label = unique(subnar$rn),
      group = "Gene",
      shape = "box",
      value = 1,
      color = "#d9d9d9"
    )
    
    nodes <- rbind(tf_nodes, gene_nodes, fill = TRUE)
    
    # -------------------------
    # TF -> Gene edges
    # -------------------------
    tf_gene_edges <- subnar[, .(
      from = TF_ID,
      to = rn,
      arrows = "to",
      value = abs(log2FoldChange),
      width = pmax(1, 3 * abs(log2FoldChange)),
      title = paste0(
        "<b>", fifelse(is.na(TF) | TF == "", TF_ID, TF), " → ", rn, "</b><br>",
        "Direction: ", ifelse(log2FoldChange > 0, "Positive", "Negative"), "<br>",
        "log2FC: ", round(log2FoldChange, 3), "<br>",
        "padj: ", signif(padj, 3)
      ),
      edge_col = ifelse(
        log2FoldChange > 0,
        paste0("rgba(139,0,0,", input$logtrans, ")"),   # darkred
        paste0("rgba(0,0,139,", input$logtrans, ")")    # darkblue
      )
    )]
    
    tf_gene_edges[, color := lapply(
      edge_col,
      function(cc) list(color = cc, highlight = cc, hover = cc, inherit = FALSE)
    )]
    tf_gene_edges[, edge_col := NULL]
    

    # -------------------------
    # TF -> TF Harmony edges (DIRECTED)
    # -------------------------
    harmony_cutoff <- input$networkhcutoff
    k <- 5L
    hdt <- target_harmony_dt()
    
    h <- hdt[
      TF1_ID %in% hit & TF2_ID %in% hit & TF1_ID != TF2_ID &
        !is.na(Harmony) & abs(Harmony) >= harmony_cutoff,
      .(TF1_ID, TF2_ID, Harmony, Concordant_Harmony, Discordant_Harmony)
    ]
    
    if (nrow(h) > 0) {
      
      # keep top-k outgoing harmony edges per TF1
      h <- h[order(TF1_ID, -abs(Harmony))]
      h <- h[, head(.SD, k), by = TF1_ID]
      
      tf_tf_edges <- h[, .(
        from = TF1_ID,
        to = TF2_ID,
        arrows = "to",
        value = abs(Harmony),
        width = pmax(0.1, 0.8 * abs(Harmony)),
        dashes = TRUE,
        smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.15),
        title = paste0(
          "<b>", TF1_ID, " → ", TF2_ID, "</b><br>",
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
      
      tf_tf_edges[, color := lapply(
        edge_col,
        function(cc) list(color = cc, highlight = cc, hover = cc, inherit = FALSE)
      )]
      tf_tf_edges[, edge_col := NULL]
      
    } else {
      tf_tf_edges <- data.table(
        from = character(),
        to = character(),
        arrows = character(),
        value = numeric(),
        width = numeric(),
        dashes = logical(),
        smooth = list(),
        title = character(),
        color = list()
      )
    }
    
    # -------------------------
    # Combine edges
    # -------------------------
    edges <- rbind(tf_gene_edges, tf_tf_edges, fill = TRUE)
    
    # -------------------------
    # Render network
    # -------------------------
    visNetwork(nodes, edges, width = "100%", height = "900px") %>%
      visGroups(groupname = "TF", shape = "triangle") %>%
      visGroups(groupname = "Gene", shape = "box") %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = TRUE,
        collapse = TRUE
      ) %>%
      visInteraction(
        dragNodes = TRUE,
        dragView = TRUE,
        zoomView = TRUE,
        navigationButtons = TRUE
      ) %>%
      #visHierarchicalLayout(enabled = TRUE, direction = "UD", sortMethod = "directed") %>%
      visPhysics(
        solver = "forceAtlas2Based",
        stabilization = list(enabled = TRUE, iterations = 800, fit = TRUE)
      ) %>%
      visEvents(
        stabilizationIterationsDone = "function () { this.setOptions({physics: false}); }"
      )
    
    
  
  })
    
  
  
      
     
 

}

# Run the application 
shinyApp(ui = ui, server = server)
