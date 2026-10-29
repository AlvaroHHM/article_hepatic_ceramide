# 35b v2 - procesa cada GSM y guarda checkpoint (reanudable)
suppressPackageStartupMessages({
  library(hdf5r); library(Matrix); library(dplyr); library(tidyr); library(ggplot2)
})
dir.create("logs", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/35b_checkpoint", showWarnings = FALSE, recursive = TRUE)
sink("logs/35b_acquisition_log.txt", split = TRUE)
on.exit(sink())

message("=== 35b v2 (checkpoints) ===")
message("Date: ", Sys.time())

genes_int <- c("B4galt6", "S1pr1", "Irf7", "B4galt5", "Ugcg",
               "Smpd3", "Sgms1", "Asah1", "Sgpl1")

map <- read.csv("results/34_GSM_mapping_v2.csv", stringsAsFactors = FALSE)
map_use <- map %>% filter(!is.na(orig.ident), !grepl("Naive", timepoint))
message("GSMs: ", nrow(map_use))

ann <- read.csv(gzfile("data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"),
                colClasses = "character", check.names = FALSE)
if (colnames(ann)[1] == "") colnames(ann)[1] <- "barcode"
ann_mol3 <- ann %>% filter(cellTypeFalcao2018 == "MOL3_EAE") %>%
  select(barcode, orig.ident) %>% mutate(bc_clean = sub("-\\d+$", "", barcode))
rm(ann); gc(verbose = FALSE)
message("MOL3_EAE totales: ", nrow(ann_mol3))

# --- FASE 1: procesar cada GSM, guardar checkpoint ---
message("\n=== FASE 1: procesar GSMs ===")
for (i in seq_len(nrow(map_use))) {
  gsm <- map_use$GSM[i]; oid <- map_use$orig.ident[i]
  tp <- map_use$timepoint[i]; cond <- map_use$condition[i]
  ckpt <- file.path("results/35b_checkpoint", paste0(gsm, ".rds"))

  if (file.exists(ckpt)) {
    message("[", i, "/", nrow(map_use), "] ", gsm, ": SKIP (ya existe)")
    next
  }

  bc_mol3 <- ann_mol3$bc_clean[ann_mol3$orig.ident == oid]
  if (length(bc_mol3) == 0) {
    saveRDS(data.frame(), ckpt)
    message("[", i, "/", nrow(map_use), "] ", gsm, ": 0 MOL3 (empty marker)")
    next
  }

  ok <- tryCatch({
    h5_files <- list.files("data/GSE250589/h5",
                           pattern = paste0("^", gsm, "_"), full.names = TRUE)
    h5_file <- h5_files[grepl("\\.h5$", h5_files)][1]
    h5 <- H5File$new(h5_file, mode = "r")
    features <- h5[["matrix/features/name"]]$read()
    barcodes <- h5[["matrix/barcodes"]]$read()
    d <- h5[["matrix/data"]]$read()
    idx <- h5[["matrix/indices"]]$read()
    iptr <- h5[["matrix/indptr"]]$read()
    shp <- h5[["matrix/shape"]]$read()
    h5$close_all()

    bc_clean <- sub("-\\d+$", "", barcodes)
    matched <- intersect(bc_mol3, bc_clean)
    if (length(matched) == 0) {
      saveRDS(data.frame(), ckpt)
      rm(d, idx, iptr); gc(verbose = FALSE)
      message("[", i, "/", nrow(map_use), "] ", gsm, ": no match")
      FALSE
    } else {
      col_idx <- which(bc_clean %in% matched)
      j_all <- rep(seq_len(shp[2]), diff(iptr))
      keep <- j_all %in% col_idx
      mat <- sparseMatrix(i = idx[keep] + 1,
                          j = match(j_all[keep], col_idx),
                          x = d[keep], dims = c(shp[1], length(col_idx)))
      rownames(mat) <- features
      rm(d, idx, iptr, j_all, keep); gc(verbose = FALSE)

      total_umi <- colSums(mat)
      g_avail <- intersect(genes_int, features)
      mat_g <- as.matrix(mat[g_avail, , drop = FALSE])
      rm(mat); gc(verbose = FALSE)

      cpm_g <- sweep(mat_g, 2, ifelse(total_umi > 0, total_umi, 1), "/") * 1e6
      df_meta <- data.frame(GSM = gsm, orig.ident = oid, timepoint = tp,
                            condition = cond, total_umi = total_umi,
                            stringsAsFactors = FALSE)
      df_genes <- as.data.frame(t(cpm_g), check.names = FALSE)
      colnames(df_genes) <- g_avail
      out <- cbind(df_meta, df_genes)
      saveRDS(out, ckpt)
      rm(df_genes, cpm_g, mat_g); gc(verbose = FALSE)
      message("[", i, "/", nrow(map_use), "] ", gsm, " (", oid, "): ",
              length(matched), " cells OK")
      TRUE
    }
  }, error = function(e) {
    message("[", i, "/", nrow(map_use), "] ", gsm, ": ERROR ", conditionMessage(e))
    FALSE
  })
  Sys.sleep(0.2)
}

# --- FASE 2: agregar y analizar ---
message("\n=== FASE 2: analizar ===")
ckpts <- list.files("results/35b_checkpoint", pattern = "\\.rds$", full.names = TRUE)
cells_list <- lapply(ckpts, function(c) {
  x <- tryCatch(readRDS(c), error = function(e) NULL)
  if (is.data.frame(x) && nrow(x) > 0) x else NULL
})
cells_list <- cells_list[!sapply(cells_list, is.null)]
if (length(cells_list) == 0) stop("No valid checkpoints")
cells_df <- bind_rows(cells_list)
message("Total cells: ", nrow(cells_df))

long <- cells_df %>% filter(condition == "EAE") %>%
  pivot_longer(cols = all_of(intersect(genes_int, colnames(cells_df))),
               names_to = "gene", values_to = "CPM")

by_tp <- long %>% group_by(timepoint, gene) %>%
  summarise(mean_CPM = mean(CPM), median_CPM = median(CPM),
            n_cells = n(), .groups = "drop")
write.csv(by_tp, "results/35b_MOL3_CPM_by_timepoint.csv", row.names = FALSE)

early_peak <- long %>% filter(timepoint %in% c("Early", "Peak")) %>%
  group_by(gene, timepoint) %>%
  summarise(mean_CPM = mean(CPM), n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = timepoint,
              values_from = c(mean_CPM, n_cells)) %>%
  mutate(log2FC_mean = log2(mean_CPM_Peak + 1) - log2(mean_CPM_Early + 1)) %>%
  arrange(desc(log2FC_mean))
write.csv(early_peak, "results/35b_MOL3_early_vs_peak.csv", row.names = FALSE)

stats_list <- list()
for (g in unique(long$gene)) {
  e <- long$CPM[long$gene == g & long$timepoint == "Early"]
  p <- long$CPM[long$gene == g & long$timepoint == "Peak"]
  if (length(e) >= 3 && length(p) >= 3) {
    wt <- tryCatch(wilcox.test(p, e, exact = FALSE), error = function(er) NULL)
    if (!is.null(wt)) stats_list[[g]] <- data.frame(
      gene = g, n_Early = length(e), n_Peak = length(p),
      median_Early = median(e), median_Peak = median(p),
      log2FC_median = log2(median(p) + 1) - log2(median(e) + 1),
      p_value = wt$p.value)
  }
}
stats_df <- bind_rows(stats_list) %>% mutate(padj = p.adjust(p_value, "BH"))
write.csv(stats_df, "results/35b_MOL3_wilcoxon.csv", row.names = FALSE)

print(early_peak)
message("\n=== FIN ===")
