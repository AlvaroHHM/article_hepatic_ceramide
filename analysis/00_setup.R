# ============================================================
# 00_setup.R
# Instala y verifica las dependencias del repositorio.
# ============================================================
pkgs_bioc <- c("recount3","DESeq2","fgsea","KEGGREST","GEOquery",
               "limma","Seurat","tximeta")
pkgs_cran <- c("dplyr","tidyr","ggplot2","cowplot","tibble",
               "readr","pdftools","data.table","ggsci")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("[install] ", p)
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("[install] ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

message("\n=== Versiones ===")
message("R:           ", getRversion())
message("Bioconductor:", as.character(BiocManager::version()))
for (p in c(pkgs_bioc, pkgs_cran)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "missing")
  message(sprintf("  %-18s %s", p, v))
}
message("\n[OK] Setup completo")
