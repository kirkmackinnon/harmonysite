# helpers.R
# Utility functions used by the Shiny server.

## This is the Fisher Test function. You provide the number of shared targets, and then two
## sets of targets, it sets up the contingency table and runs the fisher exact test.
ftestfun <- function(shared, Tar1, Tar2) {
  totalgenes <- 32031
  inNewOnly <- nrow(Tar1) - shared
  inOldOnly <- nrow(Tar2) - shared
  notInEither <- totalgenes - (shared + inNewOnly + inOldOnly)

  contTable <- data.table(
    notA = c(notInEither, inNewOnly),
    inA = c(inOldOnly, shared))

  fisher.test(contTable, alternative = "greater")
}

## Function to return a linked list based on a given pattern (pat)
## Subsets NarPV - the DEGs data based on the pattern
## Only grabs rows where the padj is less than 0.05 - grabs the AtID (rn), Log2FC, and padj.
## Names this list based on the pattern
## For loop cycles through a series of TF ids - so all targets for multiple tfs will be returned in a named list
## if old TF ids are still present (oldll) it doesn't recalculate them
findfile <- function(pat, oldll = list()) {
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

## Harmony clustering helper function.
## Subsets harmony data, filters by user cutoffs, casts to wide matrix, computes distance, and clusters.
## This is done separately for both X and Y axes as harmony is directional.
## harmony_type: "Concordant" or "Discordant" - which harmony columns to use
## transpose: if TRUE, transposes the matrix before distance calculation (Y axis clustering)
## subdt: the pre-subsetted harmony datatable
## harmony_cutoff, intersect_cutoff, pv_cutoff: user-supplied filter values
## dist_method, clust_method: algorithm choices for dist() and hclust()
harmony_hclust <- function(subdt, harmony_type, transpose = FALSE,
                           harmony_cutoff, intersect_cutoff, pv_cutoff,
                           dist_method, clust_method) {
  harmony_col <- paste0(harmony_type, "_Harmony")
  intersect_col <- paste0(harmony_type, "_Intersect")
  pvalue_col <- paste0(harmony_type, "_PValue")

  sub <- subdt[, .(TF1, TF2,
    Harmony = get(harmony_col),
    Intersect = get(intersect_col),
    PValue = get(pvalue_col)
  )]

  sub <- sub[!is.na(Harmony)][
    abs(Harmony) > as.numeric(harmony_cutoff)][
      Intersect > as.numeric(intersect_cutoff)][
        PValue < as.numeric(pv_cutoff)][
          !is.infinite(Harmony)]

  wdt <- dcast.data.table(sub[, .(TF1, TF2, Harmony)], TF1 ~ TF2, value.var = "Harmony", fill = 0)
  mat <- as.matrix(wdt, rownames = "TF1")

  if (transpose) mat <- t(mat)

  d <- dist(mat, method = dist_method)
  hclust(d, method = clust_method)
}

## Leiden community detection on a harmony edge table (via igraph's native implementation,
## no Python/leidenAlg dependency). Uses modularity objective for stable, interpretable
## results on small graphs.
## h: data.table with TF1_ID, TF2_ID, and one or more weight columns
## weight_col: which column to use as edge weight
## cutoff: minimum weight to keep an edge
## method: kept as an argument for forward-compatibility (e.g. swapping to CPM objective
##         or another algorithm later); currently only "leiden" is implemented.
## Returns a data.table with TF_ID and module columns, or NULL if no edges pass the cutoff.
compute_modules <- function(h, weight_col = "weight", cutoff = 0.02,
                            method = "leiden") {
  h <- copy(h)
  h <- h[!is.na(get(weight_col)) & get(weight_col) >= cutoff]

  if (nrow(h) == 0) return(NULL)

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

  cl <- igraph::cluster_leiden(
    g_mod,
    weights = igraph::E(g_mod)$weight,
    objective_function = "modularity",
    n_iterations = 2
  )
  memb <- igraph::membership(cl)

  data.table(
    TF_ID = names(memb),
    module = as.integer(memb)
  )
}

## Recalculate harmony for a user-supplied target subset (Kirk's preferred behavior).
## Implements H = N x NP x P x R from the methods, computed per TF1 x TF2 pair
## using only the genes in `target_genes` and only TFs that regulate >=1 such gene.
##
## Per Kirk's preference, NP uses the SUBSET proportion: NP(A) = N / |TF_A's DEGs INSIDE the upload|
## (not the global TF_A DEG count). We symmetrize the asymmetric NP via geometric mean.
##
## Concordant / Discordant variants split shared genes by whether log2FC signs agree.
##
## fisher_universe controls the Fisher test background. Two reasonable defaults:
##   - 32031: Arabidopsis genome size (matches ftestfun's global convention)
##   - length(target_genes): the upload-restricted universe (paired with subset NP)
## We default to genome size; toggle via the argument.
##
## Returns a data.table with the same columns as the precomputed `dt`:
##   TF1, TF2, Concordant_Intersect, Discordant_Intersect, Concordant_PValue, Discordant_PValue,
##   Concordant_Correlation, Discordant_Correlation, Concordant_Harmony, Discordant_Harmony,
##   TF1_Family, TF2_Family
recalculate_subset_harmony <- function(narpv, target_genes, padj_cutoff = 0.05,
                                       fisher_universe = 32031) {
  if (length(target_genes) == 0) return(NULL)

  ## DEGs of each TF that fall inside the upload
  sub <- narpv[rn %in% target_genes & padj < padj_cutoff,
               .(TF_ID, TF, Family, rn, log2FoldChange)]
  if (nrow(sub) == 0) return(NULL)

  ## Per-TF total DEG count within the upload (the subset NP denominator)
  tf_totals <- sub[, .(n_in_upload = uniqueN(rn),
                        TF = first(TF), Family = first(Family)), by = TF_ID]
  tf_totals <- tf_totals[n_in_upload >= 1]
  if (nrow(tf_totals) < 2) return(NULL)

  ## Self-merge on rn to enumerate every (TF1, TF2) pair that shares >=1 gene
  pairs <- merge(sub, sub, by = "rn", allow.cartesian = TRUE,
                  suffixes = c(".A", ".B"))
  pairs <- pairs[TF_ID.A < TF_ID.B]
  if (nrow(pairs) == 0) return(NULL)

  ## Tag concordant vs discordant per shared gene
  pairs[, agree := sign(log2FoldChange.A) == sign(log2FoldChange.B)]

  ## Aggregate per pair: counts + per-variant correlations
  pair_stats <- pairs[, {
    n_con <- sum(agree)
    n_dis <- sum(!agree)
    r_con <- if (n_con >= 3) suppressWarnings(cor(log2FoldChange.A[agree], log2FoldChange.B[agree])) else NA_real_
    r_dis <- if (n_dis >= 3) suppressWarnings(cor(log2FoldChange.A[!agree], log2FoldChange.B[!agree])) else NA_real_
    .(TF1 = first(TF.A), TF2 = first(TF.B),
      TF1_Family = first(Family.A), TF2_Family = first(Family.B),
      Concordant_Intersect = n_con, Discordant_Intersect = n_dis,
      Concordant_Correlation = r_con, Discordant_Correlation = r_dis)
  }, by = .(TF1_ID = TF_ID.A, TF2_ID = TF_ID.B)]

  ## Bring in each TF's denominator (subset-proportion)
  pair_stats <- merge(pair_stats,
    tf_totals[, .(TF1_ID = TF_ID, n_TF1 = n_in_upload)], by = "TF1_ID")
  pair_stats <- merge(pair_stats,
    tf_totals[, .(TF2_ID = TF_ID, n_TF2 = n_in_upload)], by = "TF2_ID")

  ## Vectorized one-sided Fisher (hypergeometric) p-value per variant
  ## phyper(q-1, m, n-m, k, lower.tail=FALSE) where q=shared, m=TF1, k=TF2, n=universe
  fisher_p <- function(shared, m_a, m_b, universe) {
    p <- phyper(shared - 1, m_a, universe - m_a, m_b, lower.tail = FALSE)
    pmin(pmax(p, .Machine$double.xmin), 1)
  }

  pair_stats[, Concordant_PValue := fisher_p(Concordant_Intersect, n_TF1, n_TF2, fisher_universe)]
  pair_stats[, Discordant_PValue := fisher_p(Discordant_Intersect, n_TF1, n_TF2, fisher_universe)]

  ## H = N x NP x P x R
  ## NP symmetrized via geometric mean of NP(A) and NP(B). |R| because Harmony is positive.
  pair_stats[, NP_con := sqrt((Concordant_Intersect / n_TF1) * (Concordant_Intersect / n_TF2))]
  pair_stats[, NP_dis := sqrt((Discordant_Intersect / n_TF1) * (Discordant_Intersect / n_TF2))]

  ## Cap -log10(p) at a sane value to avoid Inf when p underflows
  pair_stats[, P_con := pmin(-log10(Concordant_PValue), 300)]
  pair_stats[, P_dis := pmin(-log10(Discordant_PValue), 300)]

  pair_stats[, Concordant_Harmony := fifelse(
    Concordant_Intersect > 0 & !is.na(Concordant_Correlation),
    Concordant_Intersect * NP_con * P_con * abs(Concordant_Correlation),
    0
  )]
  pair_stats[, Discordant_Harmony := fifelse(
    Discordant_Intersect > 0 & !is.na(Discordant_Correlation),
    Discordant_Intersect * NP_dis * P_dis * abs(Discordant_Correlation),
    0
  )]

  ## Match the column layout of the precomputed dt
  pair_stats[, c("n_TF1", "n_TF2", "NP_con", "NP_dis", "P_con", "P_dis") := NULL]

  setcolorder(pair_stats, c(
    "TF1", "TF2",
    "Concordant_Intersect", "Discordant_Intersect",
    "Concordant_PValue", "Discordant_PValue",
    "Concordant_Correlation", "Discordant_Correlation",
    "Concordant_Harmony", "Discordant_Harmony",
    "TF1_Family", "TF2_Family",
    "TF1_ID", "TF2_ID"
  ))
  pair_stats[]
}

## Given multiple differential-expression result tables, compare every pair by merging on gene id, then tag each gene as
## concordant/discordant in direction of effect, caching previous pair computations if provided.
matchtidy <- function(Tar1, oldmatch = data.table()) {
  if (length(Tar1) < 2) return(data.table())
  combos <- c()
  listlength <- choose(length(Tar1), 2)
  tempagg <- vector("list", listlength)
  count <- 1

  for (i in seq_along(Tar1)) {
    for (j in i:length(Tar1)) {
      if (i == j) next
      inter <- paste0(names(Tar1)[i], "_", names(Tar1)[j])
      combos <- c(combos, inter)
      if (inter %in% oldmatch$Inter) next
      temp <- merge.data.table(Tar1[[i]], Tar1[[j]], by = "rn")
      temp[, Inter := inter]
      tempagg[[count]] <- temp
      count <- count + 1
    }
  }

  matched <- rbindlist(tempagg)
  if (nrow(matched) > 0) {
    matched[, Harmony := ifelse(sign(`log2FC.x`) == sign(`log2FC.y`), "Concordant", "Discordant")]
  }
  matched <- rbind(matched, oldmatch)
  matched <- matched[Inter %in% combos]
  return(matched)
}
