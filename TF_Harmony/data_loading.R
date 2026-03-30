# data_loading.R
# Loads preprocessed RDS files created by preprocess_data.R.
# If the RDS files don't exist, falls back to reading raw TSVs and processing inline.
# Run `Rscript preprocess_data.R` once to generate the RDS files.

## Color scheme for network graphs overlaying activator/repressor data with TF activity
testcolors <- c("red", "gray", "blue", "gray", "darkred", "darkblue")
names(testcolors) <- c("Activator", "Minimally Active", "Repressor", "Unknown", "Upregulated", "Downregulated")

## Pre-calculated harmony data with family annotations, exp17 removed, names cleaned
dt <- readRDS("Data/dt.rds")
## TF family lookup table
tfswithfamilies <- readRDS("Data/tfswithfamilies.rds")
## Pairwise harmony with fewer columns (used in Pairwise Analyses)
newdt <- readRDS("Data/newdt.rds")

## DEGs in narrow format, with color column, HSF.A4A fix, and exp17 removed
narpv <- readRDS("Data/narpv.rds")
## Unique gene IDs from DEGs, used for Target Regulation selectize
allgeneids <- readRDS("Data/allgeneids.rds")
## Ordered unique TF names from DEGs
idoptions <- readRDS("Data/idoptions.rds")
## TF ID to TF name mapping, derived from narpv
ids <- readRDS("Data/ids.rds")

## Known transcription effector domains
allteds <- readRDS("Data/allteds.rds")
## Pre-built vertices for network graphs with shapes assigned
vertices <- readRDS("Data/vertices.rds")

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