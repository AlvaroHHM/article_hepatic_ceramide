# 03_deseq2_kegg_pipeline.R (v11 - contraste manual + fix factor levels)
suppressPackageStartupMessages({
  library(recount3)
  library(DESeq2)
  library(SummarizedExperiment)
  library(org.Mm.eg.db)
  library(fgsea)
  library(KEGGREST)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Uso: Rscript 03_deseq2_kegg_pipeline.R <SRP_ID> [TEST_GROUP] [REF_GROUP]")
project_id <- args[1]
test_override <- if (length(args) >= 2) args[2] else NULL
ref_override  <- if (length(args) >= 3) args[3] else NULL

dir.create("rds", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

sink(file.path("logs", paste0(project_id, "_log.txt")), split = TRUE)
on.exit(sink())

extract_attr <- function(attrs, key) {
  sapply(attrs, function(attr) {
    if (is.na(attr) || attr == "") return(NA_character_)
    pattern <- paste0(gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", key), ";;([^|]+)")
    m <- regmatches(attr, regexpr(pattern, attr))
    if (length(m) == 0) return(NA_character_)
    trimws(sub(paste0(key, ";;"), "", m))
  })
}

tryCatch({
  message("=== PROCESANDO: ", project_id, " ===")
  if (!is.null(test_override)) message("Test override: ", test_override)
  if (!is.null(ref_override))  message("Ref override:  ", ref_override)

  rse_file <- file.path("rds", paste0(project_id, "_rse.rds"))
  if (file.exists(rse_file)) {
    rse <- readRDS(rse_file)
    message("[1/5] RSE cargado desde caché")
  } else {
    mp <- available_projects(organism = "mouse")
    info <- subset(mp, project == project_id & project_type == "data_sources")
    rse <- create_rse(info)
    saveRDS(rse, rse_file)
  }

  count_assay <- intersect(c("counts", "raw_counts"), assayNames(rse))[1]
  message("[2/5] Assay: ", count_assay)

  # --- Parseo flexible ---
  message("[3/5] Parseando condiciones")
  cd <- as.data.frame(colData(rse))
  cd$condition <- NA_character_

  if ("sra.sample_attributes" %in% colnames(cd)) {
    attrs <- cd$sra.sample_attributes
    keys_to_try <- c("treated with", "treatment", "condition",
                     "genotype/variation", "genotype", "group", "disease")
    best_key <- NULL; best_n <- 0
    for (key in keys_to_try) {
      vals <- extract_attr(attrs, key)
      n_vals <- sum(!is.na(vals) & vals != "")
      uniq <- length(unique(vals[!is.na(vals) & vals != ""]))
      if (n_vals >= 4 && uniq >= 2 && uniq > best_n) {
        best_key <- key; best_n <- uniq
      }
    }
    if (!is.null(best_key)) {
      message("Usando key: '", best_key, "' (", best_n, " grupos)")
      raw_vals <- extract_attr(attrs, best_key)
      cd$condition <- sapply(raw_vals, function(x) {
        if (is.na(x) || x == "") return(NA_character_)
        if (grepl("wild[ _]type", x, ignore.case = TRUE)) return("WT")
        if (grepl("knockout|knock-out|-/-| KO", x, ignore.case = TRUE)) return("KO")
        if (grepl("^WT$", x, ignore.case = TRUE)) return("WT")
        if (grepl("^HO$", x, ignore.case = TRUE)) return("HO")
        if (grepl("^KO$", x, ignore.case = TRUE)) return("KO")
        if (grepl("APAP", x, ignore.case = TRUE)) return("APAP")
        if (grepl("Saline", x, ignore.case = TRUE)) return("Saline")
        if (grepl("inflamed", x, ignore.case = TRUE)) return("inflamed")
        if (grepl("nodule", x, ignore.case = TRUE)) return("nodule")
        if (grepl("^day[0-9]+$", x, ignore.case = TRUE)) return(tolower(x))
        return(x)
      })
    }
  }

  if (all(is.na(cd$condition)) && "sra.sample_title" %in% colnames(cd)) {
    titles <- cd$sra.sample_title
    message("Fallback: parseando sra.sample_title")
    cd$condition <- dplyr::case_when(
      grepl("_nodule_", titles) ~ "nodule",
      grepl("_inflamed_", titles) ~ "inflamed",
      grepl("_WT_", titles) ~ "WT",
      grepl("^K[0-9]", titles) ~ "KO",
      grepl("^W[0-9]", titles) ~ "WT",
      grepl("DSSd[0-9]+", titles) ~ sub("DSSd([0-9]+).*", "day\\1", titles),
      TRUE ~ NA_character_
    )
  }

  cd$condition <- factor(cd$condition)
  message("\nCondición detectada:")
  print(table(cd$condition, useNA = "ifany"))

  # --- Elegir contraste ---
  group_counts <- table(cd$condition[!is.na(cd$condition)])
  message("\nGrupos disponibles:")
  print(group_counts)

  if (!is.null(test_override) && !is.null(ref_override)) {
    ref <- ref_override; test <- test_override
    if (!(ref %in% names(group_counts))) stop("Ref '", ref, "' no existe")
    if (!(test %in% names(group_counts))) stop("Test '", test, "' no existe")
    message("Usando contraste manual: ", test, " vs ", ref)
  } else {
    levels_cond <- names(sort(group_counts, decreasing = TRUE))
    if ("Saline" %in% levels_cond && "APAP" %in% levels_cond) {
      ref <- "Saline"; test <- "APAP"
    } else if ("WT" %in% levels_cond && "KO" %in% levels_cond) {
      ref <- "WT"; test <- "KO"
    } else if ("WT" %in% levels_cond && "HO" %in% levels_cond) {
      ref <- "WT"; test <- "HO"
    } else if ("WT" %in% levels_cond && "inflamed" %in% levels_cond) {
      ref <- "WT"; test <- "inflamed"
    } else {
      ref <- levels_cond[length(levels_cond)]; test <- levels_cond[1]
    }
    message("Contraste auto: ", test, " vs ", ref)
  }

  # --- CRÍTICO: filtrar SOLO los 2 grupos del contraste ANTES de redefinir niveles ---
  keep <- cd$condition %in% c(ref, test)
  cd_filt <- cd[keep, ]
  rse_filt <- rse[, keep]

  message("\nMuestras en el contraste: ", ncol(rse_filt),
          " (", ref, "=", sum(cd_filt$condition == ref),
          ", ", test, "=", sum(cd_filt$condition == test), ")")

  # --- DESeq2 ---
  counts_mat <- as.matrix(assay(rse_filt, count_assay))
  max_val <- max(counts_mat, na.rm = TRUE)
  if (max_val > 2^31 - 1) counts_mat <- counts_mat / ceiling(max_val / (2^31 - 1))
  counts_mat <- round(counts_mat)
  storage.mode(counts_mat) <- "integer"

  cd_filt$condition <- factor(as.character(cd_filt$condition), levels = c(ref, test))
  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                 colData = cd_filt,
                                 design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", test, ref))
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  res_df$ensembl_clean <- sub("\\..*$", "", res_df$gene_id)
  write.csv(res_df,
            file.path("results", paste0(project_id, "_DESeq2_", test, "_vs_", ref, ".csv")),
            row.names = FALSE)
  message("[4/5] DESeq2 completado: ",
          sum(res_df$padj < 0.05, na.rm = TRUE), " genes padj<0.05")

  # --- GSEA ---
  message("\n[5/5] GSEA KEGG")
  kegg_links <- keggLink("pathway", "mmu")
  kegg_pathways <- split(names(kegg_links), kegg_links)
  names(kegg_pathways) <- sub("^path:", "", names(kegg_pathways))
  kegg_pathways <- lapply(kegg_pathways, function(x) sub("^mmu:", "", x))
  kegg_names <- keggList("pathway", "mmu")
  names(kegg_names) <- sub("^path:", "", names(kegg_names))

  all_kegg_entrez <- unique(unlist(kegg_pathways))
  entrez_to_ens <- suppressMessages(
    mapIds(org.Mm.eg.db, keys = all_kegg_entrez,
           column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first")
  )
  entrez_to_ens <- entrez_to_ens[!is.na(entrez_to_ens)]

  kegg_pathways_ens <- lapply(kegg_pathways, function(genes) {
    ens <- entrez_to_ens[genes]; unique(ens[!is.na(ens)])
  })
  kegg_pathways_ens <- kegg_pathways_ens[sapply(kegg_pathways_ens, length) >= 15]

  ranks <- res_df$stat
  names(ranks) <- res_df$ensembl_clean
  ranks <- ranks[!is.na(ranks)]
  set.seed(42)
  ranks <- ranks + rnorm(length(ranks), 0, 1e-10)
  ranks <- sort(ranks, decreasing = TRUE)

  message("Corriendo fgsea con ", length(ranks), " genes...")
  fgsea_res <- suppressWarnings(
    fgsea(pathways = kegg_pathways_ens, stats = ranks,
          minSize = 15, maxSize = 500, nPermSimple = 10000)
  )

  fgsea_df <- as.data.frame(fgsea_res, stringsAsFactors = FALSE)
  fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, function(x) paste(x, collapse = ";"))
  fgsea_df$pathway_name <- kegg_names[fgsea_df$pathway]
  fgsea_df <- fgsea_df[order(fgsea_df$padj), c("pathway", "pathway_name",
                                               "size", "ES", "NES", "pval", "padj", "leadingEdge")]
  write.csv(fgsea_df,
            file.path("results", paste0(project_id, "_GSEA_", test, "_vs_", ref, ".csv")),
            row.names = FALSE)
  message("Pathways significativos (padj<0.05): ",
          sum(fgsea_df$padj < 0.05, na.rm = TRUE))

  hipotesis <- fgsea_df[fgsea_df$padj < 0.05 &
                          grepl("bile|sphingolipid|ceramide|PPAR|fatty acid|lipid|TNF|NF-kappa|chemokine|inflamm",
                                fgsea_df$pathway_name, ignore.case = TRUE), ]

  if (nrow(hipotesis) > 0) {
    message("\n★ RUTAS RELEVANTES (", test, " vs ", ref, "):")
    print(head(hipotesis[, c("pathway_name", "NES", "pval", "padj")], 20))
  } else {
    message("(Ninguna ruta de tu hipótesis con padj<0.05)")
    sig_all <- fgsea_df[fgsea_df$padj < 0.05, ]
    message("Top 10 rutas significativas:")
    print(head(sig_all[, c("pathway_name", "NES", "padj")], 10))
  }

  summary_df <- data.frame(
    project_id = project_id,
    comparison = paste0(test, "_vs_", ref),
    n_ref = sum(cd_filt$condition == ref),
    n_test = sum(cd_filt$condition == test),
    n_sig_genes = sum(res_df$padj < 0.05, na.rm = TRUE),
    n_sig_strict = sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE),
    n_pathways_total = nrow(fgsea_df),
    n_pathways_sig = sum(fgsea_df$padj < 0.05, na.rm = TRUE)
  )
  write.csv(summary_df,
            file.path("results", paste0(project_id, "_summary_", test, "_vs_", ref, ".csv")),
            row.names = FALSE)
  message("\n=== RESUMEN FINAL ===")
  print(summary_df)

}, error = function(e) {
  message("\n❌ ERROR: ", conditionMessage(e))
})
