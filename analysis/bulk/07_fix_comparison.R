# Fix del error de select() y generacion del CSV comparativo
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr)
})

panel <- read.csv("results/34_DESeq2_liver_EAE_panel.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

mdr2_reference <- tribble(
  ~gene,    ~log2FC_Mdr2KO, ~padj_Mdr2KO,
  "Smpd3",  4.80,           2e-58,
  "Sgms1",  0.01,           0.98,
  "Sgms2",  -0.23,          0.15,
  "B4galt6",3.77,           1e-43,
  "B4galt5",1.89,           1e-20,
  "Cers2",  -0.15,          0.43
)

comparison <- mdr2_reference %>%
  left_join(panel %>%
              dplyr::select(gene,
                            log2FC_Peak_EAE_vs_Ctrl, padj_Peak_EAE_vs_Ctrl,
                            sig_Peak_EAE_vs_Ctrl,
                            log2FC_Pre_EAE_vs_Ctrl,  padj_Pre_EAE_vs_Ctrl,
                            sig_Pre_EAE_vs_Ctrl),
            by = "gene")

cat("\n=== MDR2-KO vs EAE ===\n")
print(as.data.frame(comparison))

write.csv(comparison, "results/34_comparison_Mdr2KO_vs_EAE.csv",
          row.names = FALSE)
cat("\n  -> results/34_comparison_Mdr2KO_vs_EAE.csv\n")
