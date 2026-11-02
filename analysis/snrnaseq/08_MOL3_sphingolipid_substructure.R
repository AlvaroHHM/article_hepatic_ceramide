# ============================================================
# 44b_MOL3_sphingolipid_fixed.R
# Version corregida: usa sparseMatrix() en lugar de
# acceso directo a indptr (bug CSC vs CSR resuelto)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Matrix)
  library(hdf5r)
})

sink("logs/44b_MOL3_sphingolipid_log.txt", split = TRUE)
on.exit(sink())

SEP <- "============================================================"
cat(SEP, "\n")
cat("PERFIL DE ESFINGOLIPIDOS POR SUBPOBLACION DE MOL3_EAE (v2)\n")
cat("Fecha:", format(Sys.time()), "\n")
cat(SEP, "\n\n")

# ---- 1. Cargar anotacion ----
cat("--- 1. CARGANDO ANOTACION ---\n")
annot <- fread("data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz",
               data.table = FALSE, check.names = FALSE)
barcodes <- annot[[1]]
rownames(annot) <- barcodes
annot <- annot[, -1]

mol3_barcodes <- barcodes[grepl("MOL3_EAE", annot$cellTypeFalcao2018)]
cat(sprintf("  MOL3_EAE: %d celulas\n\n", length(mol3_barcodes)))

mol3_meta <- annot[mol3_barcodes, c("cellType", "orig.ident",
                                     "DamageScore", "InterferonScore")]
mol3_meta$subpop <- case_when(
  mol3_meta$cellType == "MOL2"  ~ "MOL2_like",
  mol3_meta$cellType == "MOL56" ~ "MOL56_like",
  TRUE                          ~ "Other"
)
mol3_meta$stage <- case_when(
  grepl("^Early", mol3_meta$orig.ident) ~ "Early",
  grepl("^Peak",  mol3_meta$orig.ident) ~ "Peak",
  grepl("^Late",  mol3_meta$orig.ident) ~ "Late",
  TRUE ~ "Other"
)
mol3_meta$barcode <- rownames(mol3_meta)

# ---- 2. Cargar mapeo GSM -> orig.ident ----
cat("--- 2. CARGANDO MAPEO GSM ---\n")
gsm_map <- read.csv("results/34_GSM_mapping_v2.csv", stringsAsFactors = FALSE)
gsm_col <- grep("GSM|gsm", colnames(gsm_map), value = TRUE)[1]
orig_col <- "orig.ident"
cat(sprintf("  Mapeo: %d GSMs\n\n", nrow(gsm_map)))

# ---- 3. Genes de interes ----
genes_target <- c("B4galt6", "Smpd3", "Sgms1", "Ugcg", "S1pr1",
                  "Sgpl1", "Asah1", "B4galt5")
cat("Genes:", paste(genes_target, collapse = ", "), "\n\n")

# ---- 4. Funcion corregida de extraccion ----
extract_genes_from_h5 <- function(h5_file, genes_want) {
  h5 <- H5File$new(h5_file, mode = "r")
  on.exit(tryCatch(h5$close_all(), error = function(e) NULL))

  # Leer estructura con h5read (mas robusto)
  gene_names <- as.character(h5[["matrix/features/name"]][])
  cell_bcs   <- as.character(h5[["matrix/barcodes"]][])
  shape      <- as.integer(h5[["matrix/shape"]][])
  n_genes    <- shape[1]
  n_cells    <- shape[2]

  # Filtrar genes
  gene_idx <- which(gene_names %in% genes_want)
  if (length(gene_idx) == 0) return(NULL)

  # Leer componentes CSC completos
  data_vec <- as.numeric(h5[["matrix/data"]][])
  indices  <- as.integer(h5[["matrix/indices"]][])   # 0-indexed ROW indices
  indptr   <- as.integer(h5[["matrix/indptr"]][])    # length = n_cells + 1

  # Construir sparse matrix CSC con sparseMatrix
  # 'i' = row indices (0-based -> 1-based)
  # 'p' = column pointers (length n_cells + 1)
  # 'x' = values
  # dims = c(n_genes, n_cells)
  mat <- sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = data_vec,
    dims = c(n_genes, n_cells),
    index1 = TRUE,
    giveCsparse = TRUE
  )

  # Extraer solo los genes de interes (por filas)
  sub_mat <- mat[gene_idx, , drop = FALSE]
  rownames(sub_mat) <- gene_names[gene_idx]
  colnames(sub_mat) <- cell_bcs

  # Liberar memoria
  rm(mat); gc(verbose = FALSE)

  return(list(counts = sub_mat, barcodes = cell_bcs))
}

# ---- 5. Loop sobre H5 ----
cat("--- 3. PROCESANDO 28 H5 ---\n")
h5_files <- list.files("data/GSE250589/h5", pattern = "\\.h5$",
                       full.names = TRUE)
cat(sprintf("  Archivos: %d\n\n", length(h5_files)))

all_counts <- list()
all_meta   <- list()

for (i in seq_along(h5_files)) {
  f <- h5_files[i]
  gsm <- sub("_.*$", "", basename(f))
  orig <- gsm_map[[orig_col]][gsm_map[[gsm_col]] == gsm]
  if (length(orig) == 0) {
    cat(sprintf("[%2d/%d] %s [sin mapeo]\n", i, length(h5_files), gsm))
    next
  }
  orig <- orig[1]

  # Filtrar MOL3 de este orig.ident
  mol3_this <- mol3_meta[mol3_meta$orig.ident == orig, ]
  if (nrow(mol3_this) == 0) {
    cat(sprintf("[%2d/%d] %s [0 MOL3]\n", i, length(h5_files), gsm))
    next
  }

  # Extraer
  res <- tryCatch(extract_genes_from_h5(f, genes_target),
                  error = function(e) {
                    cat(sprintf("[%2d/%d] %s [error: %s]\n",
                                i, length(h5_files), gsm, e$message))
                    NULL
                  })
  if (is.null(res)) next

  # Matching barcodes (con/sin sufijo)
  bc_h5 <- res$barcodes
  bc_annot <- mol3_this$barcode
  match_idx <- match(bc_annot, bc_h5)
  if (all(is.na(match_idx))) {
    # Intentar strip de sufijo
    bc_h5_strip <- sub("-.*$", "", bc_h5)
    bc_annot_strip <- sub("-.*$", "", bc_annot)
    match_idx <- match(bc_annot_strip, bc_h5_strip)
  }
  matched <- !is.na(match_idx)
  n_matched <- sum(matched)
  cat(sprintf("[%2d/%d] %s: %d/%d barcodes matched\n",
              i, length(h5_files), gsm, n_matched, length(bc_annot)))

  if (n_matched > 0) {
    counts_matched <- res$counts[, match_idx[matched], drop = FALSE]
    colnames(counts_matched) <- bc_annot[matched]
    all_counts[[gsm]] <- counts_matched
    all_meta[[gsm]] <- mol3_this[matched, ]
  }
}

# ---- 6. Combinar ----
cat("\n--- 4. COMBINANDO ---\n")
if (length(all_counts) == 0) {
  stop("No se pudo extraer. Ver barcodes.")
}

counts_all <- do.call(cbind, all_counts)
meta_all   <- do.call(rbind, all_meta)
cat(sprintf("  Total celulas MOL3 con counts: %d\n", ncol(counts_all)))
cat(sprintf("  Genes: %s\n", paste(rownames(counts_all), collapse = ", ")))

# ---- 7. Positividad por subpop x stage ----
cat("\n--- 5. POSITIVIDAD POR SUBPOP x STAGE ---\n")

results <- list()
for (g in rownames(counts_all)) {
  df <- data.frame(
    barcode = colnames(counts_all),
    gene = g,
    count = as.numeric(counts_all[g, ]),
    positive = as.numeric(counts_all[g, ]) > 0,
    subpop = meta_all$subpop,
    stage = meta_all$stage
  )
  results[[g]] <- df %>%
    group_by(gene, subpop, stage) %>%
    summarise(
      n_cells = n(),
      n_positive = sum(positive),
      pct_positive = 100 * mean(positive),
      .groups = "drop"
    )
}

results_df <- do.call(rbind, results)
write.csv(results_df, "results/44b_MOL3_sphingolipid_substructure.csv",
          row.names = FALSE)

cat("\n--- B4GALT6 por subpop x stage ---\n")
print(results_df[results_df$gene == "B4galt6", ])

# ---- 8. Veredicto ----
cat("\n", SEP, "\n")
cat("VEREDICTO\n")
cat(SEP, "\n\n")

b4 <- results_df[results_df$gene == "B4galt6" & results_df$stage == "Late", ]
if (nrow(b4) >= 2) {
  pct_mol2  <- b4$pct_positive[b4$subpop == "MOL2_like"]
  pct_mol56 <- b4$pct_positive[b4$subpop == "MOL56_like"]

  if (length(pct_mol2) > 0 && length(pct_mol56) > 0) {
    cat(sprintf("B4galt6+ en Late:\n"))
    cat(sprintf("  MOL3 ∩ MOL2_like:  %.1f%% (n=%d)\n",
                pct_mol2, b4$n_cells[b4$subpop == "MOL2_like"]))
    cat(sprintf("  MOL3 ∩ MOL56_like: %.1f%% (n=%d)\n",
                pct_mol56, b4$n_cells[b4$subpop == "MOL56_like"]))

    delta <- pct_mol56 - pct_mol2
    cat(sprintf("  Diferencia: %+.1f puntos\n\n", delta))

    if (delta > 10) {
      cat(">> MOL56_like tiene MAS B4galt6+\n")
      cat(">> CONSISTENTE CON HIPOTESIS: MOL56_like produce LacCer\n")
      cat(">> y tiene DamageScore bajo -> dano paracrino plausible\n")
    } else if (delta < -10) {
      cat(">> MOL2_like tiene MAS B4galt6+\n")
      cat(">> La fraccion danada tambien es la productora\n")
    } else {
      cat(">> B4galt6 similar en ambas\n")
      cat(">> B4galt6 es 'capa metabolica' que se superpone\n")
    }
  }
}

cat("\n[OK] results/44b_MOL3_sphingolipid_substructure.csv\n")
cat("\n", SEP, "\n")
