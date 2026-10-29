# ============================================================
# F5_convergence_panel.R
# Figure 5: metabolic relay liver -> CNS
# 4 panels (A-D) labeled via plot_grid; no titles inside panels
# Style: theme_cowplot + AAAS
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})

if (dir.exists("article_hepatic_ceramide")) {
  REPO <- "article_hepatic_ceramide"
} else if (basename(getwd()) == "article_hepatic_ceramide") {
  REPO <- "."
} else stop("Ejecuta desde ~/HARDVARD/recount3_analysis/ o desde el repo")

source(file.path(REPO, "analysis/figures/00_figure_style.R"))
message("[repo] ", normalizePath(REPO))

OUT <- file.path(REPO, "figures/main/F6_convergence_panel.pdf")

find_col <- function(df, cand) {
  hit <- intersect(cand, colnames(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}
gene_cands <- c("gene","gene_id","Gene","symbol","Symbol","X")

# ---- A) Liver B4GALT6 pre vs peak ----
A_data <- data.frame(
  phase  = factor(c("Pre-onset","Peak"), levels = c("Pre-onset","Peak")),
  log2FC = c(1.49, -0.11), sig = c("***","ns"),
  stringsAsFactors = FALSE)
f33 <- file.path(REPO, "results/bulk/33_liver_EAE_log2FC.csv")
if (file.exists(f33)) {
  d <- read.csv(f33, check.names = FALSE)
  g  <- find_col(d, gene_cands)
  pc <- find_col(d, c("Pre_EAE_vs_Ctrl","log2FC_pre","log2FC_Pre","pre_log2FC"))
  kc <- find_col(d, c("Peak_EAE_vs_Ctrl","log2FC_peak","log2FC_Peak","peak_log2FC"))
  if (!is.na(g) && !is.na(pc) && !is.na(kc)) {
    row <- d[d[[g]] == "B4galt6", , drop = FALSE]
    if (nrow(row) == 1) {
      A_data$log2FC <- c(as.numeric(row[[pc]]), as.numeric(row[[kc]]))
      fpre  <- file.path(REPO, "results/bulk/34_DESeq2_liver_EAE_Pre_EAE_vs_Ctrl.csv")
      fpeak <- file.path(REPO, "results/bulk/34_DESeq2_liver_EAE_Peak_EAE_vs_Ctrl.csv")
      if (file.exists(fpre) && file.exists(fpeak)) {
        dp <- read.csv(fpre,  check.names=FALSE)
        dk <- read.csv(fpeak, check.names=FALSE)
        gp <- find_col(dp, gene_cands); gk <- find_col(dk, gene_cands)
        lp <- find_col(dp, c("log2FoldChange","log2FC","logFC"))
        lk <- find_col(dk, c("log2FoldChange","log2FC","logFC"))
        pp <- find_col(dp, c("padj","adj.P.Val"))
        pk <- find_col(dk, c("padj","adj.P.Val"))
        if (!any(is.na(c(gp,gk,lp,lk,pp,pk)))) {
          rp <- dp[dp[[gp]] == "B4galt6", , drop=FALSE]
          rk <- dk[dk[[gk]] == "B4galt6", , drop=FALSE]
          if (nrow(rp)==1 && nrow(rk)==1) {
            A_data$log2FC <- c(as.numeric(rp[[lp]]), as.numeric(rk[[lk]]))
            A_data$sig <- c(
              if (as.numeric(rp[[pp]]) < 0.05) "***" else "ns",
              if (as.numeric(rk[[pk]]) < 0.05) "*"   else "ns")
          }
        }
      }
    }
  }
}
A_data$vjust <- ifelse(A_data$log2FC > 0, -0.4, 1.4)
A_data$fill_group <- factor(ifelse(A_data$log2FC > 0, "up","down"),
                             levels = c("up","down"))

pA <- ggplot(A_data, aes(x = phase, y = log2FC)) +
  geom_col(aes(fill = fill_group), width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = paste0(sprintf("%+.2f", log2FC), sig),
                vjust = vjust), size = 3) +
  scale_fill_manual(values = c(up = AAAS$red, down = AAAS$blue)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  coord_cartesian(ylim = c(min(A_data$log2FC) - 0.6,
                            max(A_data$log2FC) + 0.8)) +
  theme_pub(base_size = 9) +
  labs(x = NULL, y = expression(log[2]*FC))

# ---- B) MOL3_EAE % B4galt6+ Early vs Peak ----
B_data <- data.frame(
  time = factor(c("Early","Peak"), levels = c("Early","Peak")),
  pct  = c(7.06, 17.9), stringsAsFactors = FALSE)
f36 <- file.path(REPO, "results/snrnaseq/36_positivity.csv")
if (file.exists(f36)) {
  d <- read.csv(f36, check.names = FALSE)
  if (all(c("gene","timepoint","pct_positive") %in% colnames(d))) {
    db4 <- d[d$gene == "B4galt6", , drop = FALSE]
    ev <- db4$pct_positive[db4$timepoint == "Early"]
    kv <- db4$pct_positive[db4$timepoint == "Peak"]
    if (length(ev)==1 && length(kv)==1) B_data$pct <- c(ev, kv)
  }
}
B_data$vjust <- -0.4
pB <- ggplot(B_data, aes(x = time, y = pct)) +
  geom_col(fill = AAAS$red, alpha = 0.85, width = 0.55) +
  geom_text(aes(label = sprintf("%.1f%%", pct), vjust = vjust), size = 3) +
  coord_cartesian(ylim = c(0, max(B_data$pct) * 1.3)) +
  theme_pub(base_size = 9) +
  labs(x = NULL, y = expression("% B4galt6"^"+"))

# ---- C) Damage vs Interferon ----
C_data <- data.frame(
  time            = factor(c("Early","Peak","Late"), levels=c("Early","Peak","Late")),
  Damage          = c(-0.002, 0.072, 0.038),
  Interferon      = c(0.343, 0.562, 0.112),
  stringsAsFactors = FALSE)
f36d <- file.path(REPO, "results/snrnaseq/36_damage_mol3_by_timepoint.csv")
if (file.exists(f36d)) {
  d <- read.csv(f36d, check.names = FALSE)
  t  <- find_col(d, c("timepoint","time","Timepoint","Time"))
  dc <- find_col(d, c("mean_Damage","DamageScore","damage","Damage"))
  ic <- find_col(d, c("mean_Interferon","InterferonScore","interferon","Interferon"))
  if (!any(is.na(c(t,dc,ic)))) {
    C_data <- data.frame(
      time       = factor(as.character(d[[t]]),
                          levels = c("Early","Peak","Late")),
      Damage     = as.numeric(d[[dc]]),
      Interferon = as.numeric(d[[ic]]),
      stringsAsFactors = FALSE)
  }
}
C_long <- C_data %>%
  pivot_longer(c(Damage, Interferon),
               names_to = "score", values_to = "value")

pC <- ggplot(C_long, aes(x = time, y = value,
                         color = score, group = score, shape = score)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(Damage = AAAS$red,
                                 Interferon = AAAS$blue),
                     name = NULL) +
  scale_shape_manual(values = c(Damage = 16, Interferon = 17),
                     name = NULL) +
  theme_pub(base_size = 9) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = NULL, y = "Score")

# (panel D removed: text-box schematic was redundant with A and B)

# ---- Combine with panel labels A) B) C) D) ----
final <- plot_grid(pA, pB, pC,
                   ncol           = 3,
                   labels         = c("A)", "B)", "C)"),
                   label_size     = 11,
                   label_fontface = "bold",
                   align          = "h",
                   rel_widths     = c(1, 1, 1))

ggsave(OUT, final, width = 12, height = 3.8, dpi = 300)
message("[OK] ", OUT)
