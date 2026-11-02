# ============================================================
# 11_ari_from_precomputed.R
# Recomputes ARI (Adjusted Rand Index) from precomputed CSVs,
# without requiring the 28 H5 files.
#
# Inputs:
#   - results/snrnaseq/10_expr_mat.csv       (42 genes x 13,737 cells)
#   - results/snrnaseq/10_labels.csv         (orig.ident, barcode_clean, zheng)
#   - results/snrnaseq/11_nCount.csv        (orig.ident, barcode_clean, nCount_RNA)
#
# Output:
#   - results/snrnaseq/11_ari_from_precomputed.csv
#   - logs/snrnaseq/11_ari_from_precomputed_log.txt
#
# Verifies: ARI (k=2) vs Zheng = 0.402 (manuscript: 0.40)
# Normalization: log1p(counts / nCount_RNA * 1e4), PCA scale.=FALSE
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

SEP <- "============================================================"
dir.create("logs/snrnaseq", recursive = TRUE, showWarnings = FALSE)
sink("logs/snrnaseq/11_ari_from_precomputed_log.txt", split = TRUE); on.exit(sink())

cat(SEP, "\n10b: ARI from precomputed CSVs (no H5 needed)\n")
cat("Date:", format(Sys.time()), "\n", SEP, "\n\n")

expr <- fread("results/snrnaseq/10_expr_mat.csv", data.table = FALSE,
              check.names = FALSE, header = TRUE)
rownames(expr) <- expr[[1]]; expr <- expr[, -1]
cat("Expression matrix:", nrow(expr), "x", ncol(expr), "\n")

labels <- fread("results/snrnaseq/10_labels.csv", data.table = FALSE, check.names = FALSE)
labels$key <- paste(labels$orig.ident, labels$barcode_clean, sep = "_")
cat("Labels:", nrow(labels), "cells\n")

ncount <- fread("results/snrnaseq/11_nCount.csv", data.table = FALSE,
                header = FALSE, col.names = c("orig.ident", "barcode_clean", "nCount_RNA"))
ncount$key <- paste(ncount$orig.ident, ncount$barcode_clean, sep = "_")
ncount$nCount_RNA <- suppressWarnings(as.numeric(ncount$nCount_RNA))
cat("nCount rows:", nrow(ncount), "\n\n")

common <- intersect(colnames(expr), labels$key)
common <- intersect(common, ncount$key)
common <- common[!is.na(common)]
cat("Cells in expr + labels + nCount:", length(common), "\n")

expr      <- expr[, common, drop = FALSE]
labels_m  <- labels[match(common, labels$key), ]
ncount_v  <- ncount$nCount_RNA[match(common, ncount$key)]

ok <- !is.na(ncount_v) & ncount_v > 0 & labels_m$zheng %in% c("MOL2_like", "MOL56_like")
expr     <- expr[, ok, drop = FALSE]
labels_m <- labels_m[ok, ]
ncount_v <- ncount_v[ok]
cat("After filter (MOL2/MOL56 + valid nCount):", ncol(expr), "cells\n")
cat("  MOL2_like:", sum(labels_m$zheng == "MOL2_like"),
    " | MOL56_like:", sum(labels_m$zheng == "MOL56_like"), "\n\n")

sd_vec <- apply(expr, 1, sd)
expr <- expr[sd_vec > 0, , drop = FALSE]
cat("Genes with variance:", nrow(expr), "\n\n")

cat("--- Normalization: log1p(counts / nCount_RNA * 1e4) ---\n")
expr_norm <- log1p(t(t(as.matrix(expr)) / ncount_v * 1e4))

cat("--- PCA (10 PCs, scale.=FALSE) ---\n")
set.seed(42)
pca <- prcomp(t(expr_norm), center = TRUE, scale. = FALSE)
var_exp <- 100 * sum(pca$sdev[1:10]^2) / sum(pca$sdev^2)
cat(sprintf("  PC1-10 var explained: %.1f%%\n\n", var_exp))

ari_fn <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sc <- function(x) sum(choose(x, 2))
  idx <- sc(tab); ei <- sc(rowSums(tab)); ej <- sc(colSums(tab))
  (idx - ei * ej / choose(n, 2)) / ((ei + ej) / 2 - ei * ej / choose(n, 2))
}
zheng_bin <- ifelse(labels_m$zheng == "MOL2_like", 1, 2)

cat("--- k-means k=2..6 on PC1-10 (default params: nstart=1, iter.max=10) ---\n")
aris <- numeric(5)
for (k in 2:6) {
  set.seed(42)
  km <- kmeans(pca$x[, 1:10], centers = k)
  aris[k - 1] <- ari_fn(km$cluster, zheng_bin)
  cat(sprintf("  k=%d: ARI=%.3f\n", k, aris[k - 1]))
}
cat("\n")

out <- data.frame(
  metric = c("ari_k2", "ari_k3", "ari_k4", "ari_k5", "ari_k6",
             "n_cells", "n_mol2", "n_mol56", "n_genes", "pc_var_explained_pct"),
  value  = c(aris, ncol(expr), sum(labels_m$zheng == "MOL2_like"),
             sum(labels_m$zheng == "MOL56_like"), nrow(expr), var_exp)
)
write.csv(out, "results/snrnaseq/11_ari_from_precomputed.csv", row.names = FALSE)
cat("[OK] results/snrnaseq/11_ari_from_precomputed.csv\n\n")

cat(SEP, "\nVERIFICATION\n", SEP, "\n")
expected <- 0.402
delta <- abs(aris[1] - expected)
cat(sprintf("  ARI (k=2) actual    = %.4f\n", aris[1]))
cat(sprintf("  ARI (k=2) expected  = %.3f (manuscript)\n", expected))
cat(sprintf("  delta               = %.4f %s\n", delta,
            if (delta < 0.005) "OK" else "MISMATCH"))
cat("\n", SEP, "\nFIN\n", SEP, "\n")
