# preprocess_data.R
# Run this script ONCE to convert raw TSV/CSV data files into preprocessed RDS files.
# This performs all the column renaming, merging, filtering, and cleanup that
# data_loading.R used to do on every app startup.
#
# Usage: Rscript preprocess_data.R   (from the TF_Harmony directory)
#
# After running, the app will load the RDS files instead of the raw TSVs,
# cutting startup time significantly (especially for the 669MB narpv and 95MB cte files).

library(data.table)
library(universalmotif)

cat("Preprocessing data files into RDS format...\n")

## --- Harmony table ---
cat("  Processing harmony table...\n")
dt <- fread("Data/harmonytable_nobatch_nebs.tsv")

setnames(dt,
  c("TF", "Intersect_Concordant", "Intersect_Discordant",
    "PValue_Concordant", "PValue_Discordant",
    "Correlation_Concordant", "Correlation_Discordant",
    "Harmony_Concordant", "Harmony_Discordant"),
  c("TF1", "Concordant_Intersect", "Discordant_Intersect",
    "Concordant_PValue", "Discordant_PValue",
    "Concordant_Correlation", "Discordant_Correlation",
    "Concordant_Harmony", "Discordant_Harmony"))

## Signed Harmony score
setDT(dt)
dt[, Harmony := Concordant_Harmony - Discordant_Harmony]

dt <- dt[!(TF1 == TF2)]
dt[, TF1 := sub("-.*", "", TF1)]
dt[, TF2 := sub("-.*", "", TF2)]

## --- TF families ---
tfswithfamilies <- fread("Data/tfidswithfamilies.tsv")

setkey(dt, TF1)
setkey(tfswithfamilies, Name)
dt <- dt[tfswithfamilies[, .(Name, TF1_Family = Family)]]
setkey(dt, TF2)
dt <- dt[tfswithfamilies[, .(Name, TF2_Family = Family)]]

## --- DEGs (narpv) ---
cat("  Processing DEGs (669MB, this may take a moment)...\n")
narpv <- fread("Data/newnarpv_nebs_nobatch.tsv")
narpv[, Color := ifelse(sign(log2FoldChange) < 0, "Blue", "Red")]
narpv[TF == "HSF.A4A", TF_ID := "AT4G18880"]

## Remove experiment 17 from both narpv and dt
exp17nogo <- narpv[EXP == "1-17", unique(TF)]
narpv <- narpv[EXP != "1-17"]
dt <- dt[!(TF1 %in% exp17nogo)]
dt <- dt[!(TF2 %in% exp17nogo)]

## --- Derived objects from narpv ---
allgeneids <- unique(narpv$rn)
idoptions <- narpv[order(TF), unique(TF)]

## ids table (TF ID to TF name mapping, derived from narpv instead of DEGS folder)
ids <- unique(narpv[, .(Ids = TF_ID, TF)])

## --- Vertices for network graphs ---
cat("  Processing vertices...\n")
allteds <- fread("Data/allteds.tsv")
gis <- narpv[, .(unique(rn))]
vertices <- merge(gis, allteds, by.x = "V1", by.y = "Locus", all.x = T)
vertices <- merge(vertices, unique(narpv[, .(TF_ID, TF)]), by.x = "V1", by.y = "TF_ID", all.x = T)

fulltfids <- fread("Data/all_ath_tf_ids.tsv")
vertices <- merge(vertices, fulltfids, by.x = "V1", by.y = "Gene_ID", all.x = T)
vertices[!(is.na(Family)), shape := "triangle"]
vertices[(is.na(Family)), shape := "circle"]
vertices[!(is.na(TF)), shape := "triangle"]
vertices <- unique(vertices)

## --- Pairwise harmony (newdt) ---
newdt <- dt[, .(
  TF1 = sub("AT.*?_", "", TF1),
  TF2 = sub("AT.*?_", "", TF2),
  Concordant = Concordant_Harmony,
  Discordant = Discordant_Harmony
)][order(-Concordant)]

## --- Cell type expression ---
cat("  Processing cell type expression (95MB)...\n")
cte <- fread("Data/cellTypeExpression.tsv", key = "Gene ID")

## --- JIT datasets ---
cat("  Processing JIT datasets...\n")
jitr <- fread("Data/JITGenes_root.csv", skip = 1,
  select = c("Gene", "FDR adjusted p-value", "First Response (Just-in-time bin)"),
  col.names = c("GeneID", "pvalue", "JIT"))
setkey(jitr, GeneID)

jits <- fread("Data/JITGenes_shoot.csv", skip = 1,
  select = c("AtID", "FDR adjusted p-value", "First Response (Just-in-time bin)"),
  col.names = c("GeneID", "pvalue", "JIT"))
setkey(jits, GeneID)

## --- PWMs / PFMs ---
cat("  Processing motif data...\n")
pwms <- readRDS("Data/pwms.RDS")
pfms <- convert_motifs(pwms, class = "motifStack-pfm")

## --- Phylogenetic tree from MSA ---
cat("  Computing phylogenetic tree from MSA...\n")
htmsa <- Biostrings::readAAMultipleAlignment("Data/msa.fa")
phylo_alignment <- msa::msaConvert(htmsa, type = "seqinr::alignment")
phylo_distance <- seqinr::dist.alignment(phylo_alignment)
phylo_tree <- ape::bionj(phylo_distance)
phylo_dend <- as.dendrogram.phylo(phylo_tree)

## --- Save all preprocessed objects ---
cat("  Saving RDS files to Data/...\n")

saveRDS(dt, "Data/dt.rds")
saveRDS(narpv, "Data/narpv.rds")
saveRDS(cte, "Data/cte.rds")
saveRDS(jitr, "Data/jitr.rds")
saveRDS(jits, "Data/jits.rds")
saveRDS(newdt, "Data/newdt.rds")
saveRDS(vertices, "Data/vertices.rds")
saveRDS(ids, "Data/ids.rds")
saveRDS(allgeneids, "Data/allgeneids.rds")
saveRDS(idoptions, "Data/idoptions.rds")
saveRDS(tfswithfamilies, "Data/tfswithfamilies.rds")
saveRDS(allteds, "Data/allteds.rds")
saveRDS(pwms, "Data/pwms_processed.rds")
saveRDS(pfms, "Data/pfms.rds")
saveRDS(phylo_dend, "Data/phylo_dend.rds")

cat("Done! RDS files saved to Data/\n")
cat("You can now restart the app — it will use the preprocessed files.\n")
