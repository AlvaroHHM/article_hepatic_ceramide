suppressPackageStartupMessages({
  library(dplyr)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

setwd("~/HARDVARD/recount3_analysis")

f <- "data/GSE312673/GSE312673_counts_matrix.tsv.gz"
if (!file.exists(f)) stop("Missing: ", f)

counts <- read.delim(gzfile(f), check.names = FALSE)
cat("[OK] counts:", nrow(counts), "x", ncol(counts), "\n")
cat("[OK] columnas:\n"); print(colnames(counts))

# Extraer gene_id y strip versión
gene_ids <- sub("\\..*$", "", as.character(counts[[1]]))
mat <- as.matrix(counts[, -1])
rownames(mat) <- gene_ids

# Mapear solo los genes que necesitamos
targets_ens <- mapIds(org.Mm.eg.db,
                      keys = c("Smpd3", "Sgms1"),
                      column = "ENSEMBL",
                      keytype = "SYMBOL",
                      multiVals = "first")

cat("\n[OK] ENSEMBL para Smpd3:", targets_ens["Smpd3"])
cat("\n[OK] ENSEMBL para Sgms1:", targets_ens["Sgms1"], "\n")

# Columnas de grupos
col_nc_pre  <- grep("pre-onset.*NC",  colnames(mat), value = TRUE)
col_eae_pre <- grep("pre-onset.*EAE", colnames(mat), value = TRUE)
col_nc_pk   <- grep("peak.*NC",       colnames(mat), value = TRUE)
col_eae_pk  <- grep("peak.*EAE",      colnames(mat), value = TRUE)

cat("\n[OK] Pre-onset NC:",  paste(col_nc_pre,  collapse = ", "), "\n")
cat("[OK] Pre-onset EAE:", paste(col_eae_pre, collapse = ", "), "\n")
cat("[OK] Peak NC:",       paste(col_nc_pk,   collapse = ", "), "\n")
cat("[OK] Peak EAE:",      paste(col_eae_pk,  collapse = ", "), "\n")

# CPM
library_size <- colSums(mat)
cpm <- t(t(mat) / library_size * 1e6)

# LOO function
loo_log2fc <- function(values_nc, values_eae) {
  base <- log2(mean(values_eae) + 1) - log2(mean(values_nc) + 1)
  loo_results <- c()
  for (i in seq_along(values_nc)) {
    loo_results <- c(loo_results,
      log2(mean(values_eae) + 1) - log2(mean(values_nc[-i]) + 1))
  }
  for (i in seq_along(values_eae)) {
    loo_results <- c(loo_results,
      log2(mean(values_eae[-i]) + 1) - log2(mean(values_nc) + 1))
  }
  list(base = base, loo = loo_results, 
       range = range(loo_results),
       flips = any(sign(loo_results) != sign(base)))
}

out <- data.frame()
for (gene in c("Smpd3", "Sgms1")) {
  ens <- targets_ens[gene]
  if (is.na(ens) || !(ens %in% rownames(cpm))) {
    cat("[WARN]", gene, "no encontrado en counts\n")
    next
  }
  
  for (tp in c("pre-onset", "peak")) {
    nc_cols  <- if (tp == "pre-onset") col_nc_pre  else col_nc_pk
    eae_cols <- if (tp == "pre-onset") col_eae_pre else col_eae_pk
    
    nc_vals  <- cpm[ens, nc_cols]
    eae_vals <- cpm[ens, eae_cols]
    
    res <- loo_log2fc(nc_vals, eae_vals)
    
    cat(sprintf("\n=== %s / %s ===\n", gene, tp))
    cat("  NC  (n=", length(nc_vals), "): ", 
        paste(round(nc_vals, 1), collapse = ", "), "\n")
    cat("  EAE (n=", length(eae_vals), "): ", 
        paste(round(eae_vals, 1), collapse = ", "), "\n")
    cat(sprintf("  log2FC base = %+.3f\n", res$base))
    cat(sprintf("  LOO range   = [%+.3f, %+.3f]\n", res$range[1], res$range[2]))
    cat(sprintf("  Sign flip in any LOO: %s\n", res$flips))
    
    out <- rbind(out, data.frame(
      gene = gene, timepoint = tp,
      cpm_nc = paste(round(nc_vals, 1), collapse = ";"),
      cpm_eae = paste(round(eae_vals, 1), collapse = ";"),
      log2FC_base = round(res$base, 3),
      log2FC_LOO_min = round(res$range[1], 3),
      log2FC_LOO_max = round(res$range[2], 3),
      sign_flip = res$flips,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(out, "results/bulk/35_LOO_EAE.csv", row.names = FALSE)
cat("\n[OK] results/bulk/35_LOO_EAE.csv\n")
