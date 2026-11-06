suppressPackageStartupMessages({ library(dplyr) })
setwd("~/HARDVARD/recount3_analysis")

# Paths verificados
f_mdr2 <- "article_hepatic_ceramide/results/bulk/04_liver_check.csv"
f_eae  <- "article_hepatic_ceramide/results/bulk/06_DESeq2_liver_EAE_panel.csv"
f_mash <- "results/bulk/13_GSE243976_MASH_vs_Chow.csv"

for (f in c(f_mdr2, f_eae, f_mash)) {
  if (!file.exists(f)) stop("Missing: ", f)
}

mdr2 <- read.csv(f_mdr2)
eae  <- read.csv(f_eae)
mash <- read.csv(f_mash)

# Extraer Smpd3/Sgms1 con nombres de columna exactos confirmados
get_mdr2 <- function(gene) {
  r <- mdr2[mdr2$gene == gene, ]
  c(lfc = r$log2FoldChange[1], padj = r$padj[1])
}

get_eae <- function(gene, which = "pre") {
  r <- eae[eae$gene == gene, ]
  if (which == "pre") {
    c(lfc = r$log2FC_Pre_EAE_vs_Ctrl[1], padj = r$padj_Pre_EAE_vs_Ctrl[1])
  } else {
    c(lfc = r$log2FC_Peak_EAE_vs_Ctrl[1], padj = r$padj_Peak_EAE_vs_Ctrl[1])
  }
}

get_mash <- function(gene) {
  r <- mash[mash$gene == gene, ]
  c(lfc = r$log2FC_MASH_vs_Chow[1], padj = r$padj[1])
}

results <- data.frame()

for (g in c("Smpd3", "Sgms1")) {
  v_mdr2 <- get_mdr2(g)
  v_pre  <- get_eae(g, "pre")
  v_peak <- get_eae(g, "peak")
  v_mash <- get_mash(g)
  
  padj <- c(v_mdr2["padj"], v_pre["padj"], v_peak["padj"], v_mash["padj"])
  lfc  <- c(v_mdr2["lfc"],  v_pre["lfc"],  v_peak["lfc"],  v_mash["lfc"])
  models <- c("Mdr2-KO", "EAE_pre", "EAE_peak", "MASH")
  
  cat("\n=== ", g, " ===\n")
  for (i in seq_along(models)) {
    cat(sprintf("  %-10s: log2FC = %+7.3f, padj = %.3e\n", 
                models[i], lfc[i], padj[i]))
    results <- rbind(results, data.frame(
      gene = g, model = models[i],
      log2FC = round(lfc[i], 4),
      padj = signif(padj[i], 4),
      stringsAsFactors = FALSE
    ))
  }
  
  # Fisher sobre padj (proxy conservador)
  valid <- !is.na(padj) & padj > 0 & padj < 1
  X2 <- -2 * sum(log(padj[valid]))
  df <- 2 * sum(valid)
  p_comb <- pchisq(X2, df = df, lower.tail = FALSE)
  
  dir_consistent <- all(sign(lfc[valid]) == sign(lfc[valid][1]))
  
  cat(sprintf("  Fisher combined: X² = %.2f, df = %d, p = %.3e\n", X2, df, p_comb))
  cat(sprintf("  Direction consistent: %s\n", dir_consistent))
  
  results <- rbind(results, data.frame(
    gene = g, model = "Fisher_combined",
    log2FC = NA,
    padj = signif(p_comb, 4),
    stringsAsFactors = FALSE
  ))
}

write.csv(results, "results/bulk/36_fisher_3models.csv", row.names = FALSE)
cat("\n[OK] results/bulk/36_fisher_3models.csv\n")
