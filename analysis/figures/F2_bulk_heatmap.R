# ============================================================
# F2_bulk_heatmap.R
# Figure 2: SMPD3 induced, SMS1 flat across liver and CNS
# Style: theme_cowplot + AAAS palette, no title
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})

# --- repo root ---
if (dir.exists("article_hepatic_ceramide")) {
  REPO <- "article_hepatic_ceramide"
} else if (basename(getwd()) == "article_hepatic_ceramide") {
  REPO <- "."
} else stop("Ejecuta desde ~/HARDVARD/recount3_analysis/ o desde el repo")

source(file.path(REPO, "analysis/figures/00_figure_style.R"))
message("[repo] ", normalizePath(REPO))

OUT    <- file.path(REPO, "figures/main/F2_bulk_heatmap.pdf")
OUTCSV <- file.path(REPO, "figures/main/F2_bulk_heatmap_data.csv")

find_col <- function(df, cand) {
  hit <- intersect(cand, colnames(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}
gene_cands <- c("gene","Gene","symbol","Symbol","ensembl_clean","gene_id","X")

read_de <- function(path, label) {
  full <- file.path(REPO, path)
  if (!file.exists(full)) {
    message("  [SKIP] ", path, " (no existe)"); return(NULL)
  }
  df <- read.csv(full, check.names = FALSE, stringsAsFactors = FALSE)
  g <- find_col(df, gene_cands)
  l <- find_col(df, c("log2FoldChange","log2FC","logFC"))
  p <- find_col(df, c("padj","adj.P.Val"))
  if (any(is.na(c(g,l,p)))) {
    message("  [SKIP] ", basename(path), " — cols: ",
            paste(colnames(df), collapse=", "))
    return(NULL)
  }
  gene_vals <- as.character(df[[g]])

  # Auto-detect ENSEMBL IDs and map to symbols
  if (length(gene_vals) > 0 && grepl("^ENSMUSG", gene_vals[1], ignore.case = TRUE)) {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
        requireNamespace("AnnotationDbi", quietly = TRUE)) {
      ens <- sub("[.].*$", "", gene_vals)
      sym <- suppressMessages(AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        keys    = unique(ens),
        column  = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      ))
      gene_vals <- unname(sym[ens])
      n_mapped <- sum(!is.na(gene_vals))
      message("  [map]  ", basename(path), ": ENSEMBL -> SYMBOL (",
              n_mapped, "/", length(gene_vals), " mapped)")
    } else {
      warning("org.Mm.eg.db or AnnotationDbi missing for ", basename(path))
    }
  }

  message("  [OK]   ", basename(path), " (", nrow(df), " genes)")
  data.frame(gene    = gene_vals,
             log2FC  = suppressWarnings(as.numeric(df[[l]])),
             padj    = suppressWarnings(as.numeric(df[[p]])),
             condition = label, stringsAsFactors = FALSE)
}

dfs <- list(
  read_de("results/bulk/SRP074136_DESeq2_inflamed_vs_WT.csv", "Liver (Mdr2-KO)"),
  read_de("results/bulk/SRP071304_DESeq2_pEAE_vs_WT.csv",     "CNS bulk (pEAE)"),
  read_de("results/bulk/SRP197582_DESeq2_day06_vs_day00.csv", "Colon (DSS d06)")
)
dfs <- Filter(Negate(is.null), dfs)

# --- EAE liver from 34 ---
f_pre  <- file.path(REPO, "results/bulk/34_DESeq2_liver_EAE_Pre_EAE_vs_Ctrl.csv")
f_peak <- file.path(REPO, "results/bulk/34_DESeq2_liver_EAE_Peak_EAE_vs_Ctrl.csv")
if (file.exists(f_pre) && file.exists(f_peak)) {
  dpre  <- read.csv(f_pre,  check.names = FALSE)
  dpeak <- read.csv(f_peak, check.names = FALSE)
  gp <- find_col(dpre, gene_cands);  gk <- find_col(dpeak, gene_cands)
  lp <- find_col(dpre, c("log2FoldChange","log2FC","logFC"))
  lk <- find_col(dpeak,c("log2FoldChange","log2FC","logFC"))
  pp <- find_col(dpre, c("padj","adj.P.Val"))
  pk <- find_col(dpeak,c("padj","adj.P.Val"))
  if (!any(is.na(c(gp,gk,lp,lk,pp,pk)))) {
    m <- merge(dpre[,  c(gp,lp,pp)],
               dpeak[, c(gk,lk,pk)],
               by.x = gp, by.y = gk)
    colnames(m) <- c("gene","pre_log2FC","pre_padj","peak_log2FC","peak_padj")
    dfs[[length(dfs)+1]] <- data.frame(gene=m$gene, log2FC=m$pre_log2FC,
                                        padj=m$pre_padj,
                                        condition="Liver EAE (pre-onset)",
                                        stringsAsFactors=FALSE)
    dfs[[length(dfs)+1]] <- data.frame(gene=m$gene, log2FC=m$peak_log2FC,
                                        padj=m$peak_padj,
                                        condition="Liver EAE (peak)",
                                        stringsAsFactors=FALSE)
    message("  [OK]   34 merged (", nrow(m), " genes)")
  }
}

data_long <- bind_rows(dfs)
message("[data] ", nrow(data_long), " filas")

genes_sel  <- c("Smpd3","Sgms1","Sgms2","B4galt6","B4galt5","Cers2")
main_conds <- c("Liver (Mdr2-KO)","Liver EAE (pre-onset)",
                "Liver EAE (peak)","CNS bulk (pEAE)")

panel <- data_long %>%
  filter(gene %in% genes_sel, condition %in% main_conds) %>%
  mutate(
    condition = factor(condition, levels = main_conds),
    gene      = factor(gene,      levels = genes_sel),
    stars = case_when(
      is.na(padj)  ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    ),
    label = sprintf("%.2f%s", log2FC, stars)
  )

# --- Figure: NO TITLE, publication format ---
p <- ggplot(panel, aes(x = condition, y = gene, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.6, fontface = "bold") +
  scale_fill_aaas_gradient(midpoint = 0, name = expression(log[2]*FC)) +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(x = NULL, y = NULL)

ggsave(OUT, p, width = 7, height = 3.2, dpi = 300)
write.csv(panel, OUTCSV, row.names = FALSE)
message("[OK] ", OUT)
