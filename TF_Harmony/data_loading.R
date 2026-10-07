# data_loading.R
# Loads preprocessed RDS files created by preprocess_data.R.
# If the RDS files don't exist, falls back to reading raw TSVs and processing inline.
# Run `Rscript preprocess_data.R` once to generate the RDS files.

## Pre-calculated harmony data with family annotations, exp17 removed, names cleaned
dt <- readRDS("Data/dt.rds")
## TF family lookup table
tfswithfamilies <- readRDS("Data/tfswithfamilies.rds")
## Pairwise harmony with fewer columns (used in Pairwise Analyses)
newdt <- readRDS("Data/newdt.rds")

## DEGs in narrow format, with color column, HSF.A4A fix, and exp17 removed
narpv <- readRDS("Data/narpv.rds")

## Proportion term of Harmony (NP in the methods) isn't stored in dt, so derive it here:
## shared DEGs / TF1's total DEG count (padj < 0.05). Directional like Harmony itself -
## the TF1 -> TF2 row uses TF1's count. Used to order the Global Analyses heatmaps.
tf_deg_counts <- narpv[padj < 0.05, .(n_degs = uniqueN(rn)), by = TF]
dt[tf_deg_counts, on = .(TF1 = TF), TF1_DEGs := i.n_degs]
dt[, Concordant_Proportion := Concordant_Intersect / TF1_DEGs]
dt[, Discordant_Proportion := Discordant_Intersect / TF1_DEGs]
dt[, TF1_DEGs := NULL]
## Unique gene IDs from DEGs, used for Target Regulation selectize
allgeneids <- readRDS("Data/allgeneids.rds")
## Ordered unique TF names from DEGs
idoptions <- readRDS("Data/idoptions.rds")
## PWMs for motif sorting and PFMs for motifStack
pwms <- readRDS("Data/pwms_processed.rds")
pfms <- readRDS("Data/pfms.rds")

## Cell type expression from Benfey lab
cte <- readRDS("Data/cte.rds")
## Just-in-time datasets for roots and shoots
jitr <- readRDS("Data/jitr.rds")
jits <- readRDS("Data/jits.rds")

## Pre-computed phylogenetic dendrogram from MSA (pruned per-render in server)
phylo_dend <- readRDS("Data/phylo_dend.rds")