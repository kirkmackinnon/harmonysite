# ui.R
# UI definition for TF Harmony Shiny app.
# Requires: tfswithfamilies, idoptions, allgeneids from data_loading.R
#
## UI is hierarchical starting with the top three panels:
## Global Analyses, Pairwise Analyses, and Target Regulation
## Then each of those panels has sub panels and figures.

## Network layout choices reused across multiple inputs
network_layouts <- c(
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
)

ui <- navbarPage("Landscape of TF Harmony",

  tabPanel("Global Analyses",

    inputPanel(
      ## Defining cutoffs for harmony figures - user can set harmony minimums, pvalues, intersects etc
      ## and the harmony table is subset based on these inputs
      numericInput(inputId = "familyharmonycutoff", label = "Harmony Minimum:", value = 0, min = 0, step = 0.1),
      numericInput(inputId = "familypvcutoff", label = "Fisher P-Value Maximum:", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput(inputId = "familyintersectcutoff", label = "Intersect Minimum:", value = 0, min = 0, step = 1),

      ## User can define what type of distance calculation is used prior to clustering
      selectInput(
        inputId = "distmeth",
        label = "Distance Method",
        choices = c("euclidean", "maximum", "manhattan", "canberra"),
        selected = "maximum",
        multiple = F
      ),
      ## User can define the clustering method used - this is all for the heatmaps
      selectInput(
        inputId = "clustmeth",
        label = "Clustering Method",
        choices = c("complete", "average", "ward.D", "single", "mcquitty", "median", "centroid", "ward.D2"),
        selected = "ward.D2",
        multiple = F
      ),
      ## Similarly user can define what type of algorithm is used for motif comparison
      ## These are just options from the motif package
      ## again this may be removed at a later point, I'm not sure how many people have strong opinions on
      ## motif clustering algorithms.
      selectInput(
        inputId = "motifcompmethod",
        label = "Motif Comparison Method:",
        choices = c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT", "HELL", "SEUCL", "MAN", "ALLR_LL", "WEUCL", "WPCC"),
        selected = "ALLR_LL",
        multiple = F,
        selectize = F
      )
    ),

    ## The main input on the front page - asking the user which TF or Family to subset by
    ## Starts pre-populated with all TF Families listed
    selectInput(
      inputId = "family1",
      label = "Search by TF or TF Family:",
      choices = unique(c(tfswithfamilies$Family, idoptions)),
      selected = unique(tfswithfamilies$Family),
      multiple = T,
      width = "100%",
      selectize = T
    ),

    ## Series of tabs under Global Analyses
    tabsetPanel(
      tabPanel("Harmony Heatmap",
        ## Main visual you see on load - concordant harmony heatmap, interactive plotly
        tabsetPanel(
          tabPanel("Concordant Harmony",
            plotlyOutput(outputId = "conHarmonyHeatmap", height = 2000) %>% withSpinner()
          ),
          ## Alternative to main heatmap - show discordant harmony instead, interactive plotly
          tabPanel("Discordant Harmony",
            plotlyOutput(outputId = "disHarmonyHeatmap", height = 2000) %>% withSpinner()
          )
        )
      ),

      ## Second tab under Global Analyses
      ## Graphs look at prevalence of family interaction
      ## Honestly these didn't seem all that effective to me
      ## They'll likely be removed
      tabPanel("Family Harmony",
        splitLayout(
          plotlyOutput(outputId = "relationships", height = 1000) %>% withSpinner(),
          visNetworkOutput(outputId = "familyproportion", height = 1000) %>% withSpinner()
        )
      ),

      ## Third tab under Global Analyses
      ## Allows the user to select their inputs and compare dendrograms - highlighting similarities
      ## Concordant / Discordant are harmony dendrograms
      ## Phylogeny is a pre-computed MSA that is trimmed based on user selections
      ## Motif similarity is calculated on-the-fly based on pre-computed PWMs which are
      ## filtered based on user input (ie. only bZIPs)
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
          ## User input that allows you to highlight breaks
          selectInput(
            inputId = "kbreaks",
            label = "Number of breaks:",
            choices = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
            selected = 1,
            multiple = F
          )
        ),
        ## Tanglegram plot output based on above selections
        plotOutput(outputId = "tanglegram", height = 1000) %>% withSpinner()
      ),

      ## Raw harmony table for user to peruse and interact with
      tabPanel("Harmony Table",
        DT::dataTableOutput(outputId = "fullharmonydt") %>% withSpinner()
      ),

      ## Heatmap of motif similarity
      ## Changes based on user selection of motif comparison algorithm
      tabPanel("Motif Heatmap",
        plotlyOutput(outputId = "motifheatmap", height = 1000) %>% withSpinner()
      ),

      ## Shows the *actual* motifs allowing you to actually see how similar they are
      ## Downside is that it takes forever to load if you try to show all the motifs at once
      ## Should probably be limited to like 5 tfs at a time or so
      tabPanel("Motif Dendrogram",
        plotlyOutput(outputId = "motifdend", height = 1000) %>% withSpinner()
      )
    )
  ),

  ## Pairwise Analyses - the middle panel up at the top
  ## Defines all the TF x TF comparisons in more detail
  tabPanel("Pairwise Analyses",
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
        min = 500, max = 10000, step = 100, value = 1000,
        label = "Plot height:"
      ),
      ## Users can decide if they want separate pairwise plots for each combination of TF (faceted) or
      ## one plot with all the pairwise combinations shown in different colors
      checkboxInput(inputId = "harmonyinterfacet", label = "Facet by interaction:", value = T)
    ),

    tabsetPanel(id = "pairwise",
      ## First pairwise plot
      ## Shows all shared DEGs, direction of regulation, and correlation of Log2FC, computed on the fly
      ## On-the-fly was picked so that users can upload their own DEG datasets and compare to ours
      tabPanel("TF x TF",
        splitLayout(
          plotlyOutput(outputId = "harmonyplotly", height = "100%") %>% withSpinner(),
          DT::dataTableOutput("harmonyTable")
        )
      ),

      ## More motif comparisons, but this time pairwise
      ## Idea being that you can see if a lot of shared regulation is associated with a similar motif
      tabPanel("Motif Comparison", value = "pairwisemotifs",
        selectInput("motiftreestyle", label = "Tree Style:", choices = c("stack", "tree", "radialPhylog"), multiple = F),
        imageOutput(outputId = "motifplot", width = "80%", height = 800) %>% withSpinner()
      ),

      ## Are the shared DEGs expressed in a specific root cell type at baseline?
      ## Overlay of Benfey cell type data subsetted to match shared DEGs
      tabPanel("Cell Type Specificity",
        plotlyOutput(outputId = "ctplot", height = 1000) %>% withSpinner()
      ),

      ## Are the shared DEGs differentially expressed at a given time period in the NxTime datasets?
      ## Overlay of NxTime data - subsetted to match shared DEGs
      tabPanel("NxTime Overlay",
        splitLayout(
          cellWidths = c("50%", "50%"),
          plotlyOutput(outputId = "nxtplot1", height = 1000),
          plotlyOutput(outputId = "nxtplot2", height = 1000)
        )
      ),

      ## On the fly calculation of GO Terms using the GOST server
      ## Subsetted DEGs are submitted to the server at time of selection - runs even when page isn't selected
      tabPanel("GO Terms",
        sliderInput(inputId = "goplotheight", min = 500, max = 10000, step = 100, value = 1000, label = "Plot height:"),
        fluidRow(
          column(6, plotlyOutput(outputId = "goplot", width = "100%", height = 1000) %>% withSpinner()),
          column(6, div(style = "overflow-x: auto;", DT::DTOutput(outputId = "gotable", width = "100%") %>% withSpinner()))
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
            min = 0, max = 215000, value = 0,
            round = T, animate = T, width = 800
          ),
          selectInput(
            inputId = "tfregnetworkstyle",
            label = "Network Style:",
            choices = network_layouts,
            selected = "layout_with_sugiyama",
            multiple = F
          ),
          ## Are we subsetting to only shared DEGs? ** TODO - check and see what this subsets?
          checkboxInput(inputId = "networkdegs", label = "DEGs only:", value = TRUE)
        ),
        visNetworkOutput(outputId = "networkplot", height = 1000) %>% withSpinner()
      ),

      ## Separate panel showing overlapping DEGs between TFs - harmony excluded
      ## Idea here was whether we could see if we could identify the primary TF in cascading networks
      ## Degree cutoff is to limit to targets that are downstream of a minimum number of TFs
      ## ie. only keep a target if all 3 of my selected TFs regulate it
      tabPanel("DEG Networks",
        inputPanel(
          numericInput(inputId = "tfregpcutoff", label = "Adjusted P-value Cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(inputId = "tfreglfccutoff", label = "Log2FoldChange Cutoff:", value = 0, min = 0, step = 0.1),
          numericInput(inputId = "tfdegreecutoff", label = "Degree Cutoff:", value = 0, min = 0, step = 1),
          selectInput(
            inputId = "degnetworkstyle",
            label = "Network Style:",
            choices = network_layouts,
            selected = "layout_with_sugiyama",
            multiple = F
          )
        ),
        visNetworkOutput(outputId = "tfnetworkplot", height = 1000) %>% withSpinner()
      )
    )
  ),

  ## Final overarching panel - Target Regulation
  ## This is the reverse of the previous two panels.
  ## Instead of looking at which targets a TF regulates we start with the targets and look at all TFs that regulate it.
  ## You can also upload a list of targets (ie. DEGs from your own experiment) and see all the TFs
  ## that regulate those targets.
  tabPanel("Target Regulation",
    inputPanel(
      fileInput(inputId = "targetupload", label = "Upload list of Targets"),
      numericInput(inputId = "regpcutoff", label = "Adjusted P-value Cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput(inputId = "reglfccutoff", label = "Log2FoldChange Cutoff:", value = 0, min = 0, step = 0.1),
      selectInput(
        inputId = "targetnetworkstyle",
        label = "Network Style:",
        choices = network_layouts,
        selected = "layout_with_sugiyama",
        multiple = F
      )
    ),

    ## This is where the list of all targets from NarPV is used
    ## Server-side selectize for performance with large gene list
    selectizeInput(
      inputId = "allgenes",
      label = "Targets: ",
      choices = NULL,
      multiple = T,
      width = "100%"
    ),

    splitLayout(
      visNetworkOutput(outputId = "regulators", width = "50%", height = 1000) %>% withSpinner(),
      verticalLayout(
        DT::DTOutput(outputId = "tfregdt", height = 500, width = "50%"),
        DT::DTOutput(outputId = "targetregdt", height = 500, width = "50%")
      )
    ),

    includeCSS("www/style.css")
  )
)
