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

## Louvain community detection on a harmony edge table.
## h: data.table with TF1_ID, TF2_ID, and one or more weight columns
## weight_col: which column to use as edge weight
## cutoff: minimum weight to keep an edge
## Returns a data.table with TF_ID and module columns, or NULL if no edges pass the cutoff.
compute_louvain_modules <- function(h, weight_col = "weight", cutoff = 0.02) {
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

  cl <- igraph::cluster_louvain(g_mod, weights = igraph::E(g_mod)$weight)
  memb <- igraph::membership(cl)

  data.table(
    TF_ID = names(memb),
    module = as.integer(memb)
  )
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
