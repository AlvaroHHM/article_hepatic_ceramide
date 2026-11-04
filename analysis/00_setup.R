# ============================================================
# 00_setup.R
# Install and verify all dependencies required by the repository.
# Usage: Rscript analysis/00_setup.R
# ============================================================

# --- Repositories: CRAN + Bioconductor 3.18 ---
options(repos = c(
  CRAN     = "https://cloud.r-project.org",
  BioCsoft = "https://bioconductor.org/packages/3.18/bioc",
  BioCann  = "https://bioconductor.org/packages/3.18/data/annotation",
  BioCexp  = "https://bioconductor.org/packages/3.18/data/experiment"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# --- Bioconductor packages ---
pkgs_bioc <- c(
  # Core pipeline
  "recount3", "SummarizedExperiment", "DESeq2",
  "fgsea", "KEGGREST",
  # Bulk / cellular datasets
  "GEOquery", "limma",
  # snRNA-seq
  "Seurat", "hdf5r",
  # Annotation
  "org.Mm.eg.db", "AnnotationDbi"
)

# --- CRAN packages ---
pkgs_cran <- c(
  # Data manipulation
  "dplyr", "tidyr", "tibble", "data.table", "readr",
  # Visualization / reporting
  "ggplot2", "cowplot", "ggsci",
  # Statistics
  "lme4", "boot",
  # Excel reading (GSE139531 supplementary file)
  "readxl"
)

# --- Installer with per-package verification ---
install_if_missing <- function(pkg, source) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("[install] ", pkg)
    tryCatch({
      if (source == "bioc") {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg)
      }
      if (!requireNamespace(pkg, quietly = TRUE))
        warning("Package ", pkg, " did NOT install correctly")
    }, error = function(e) warning("Error installing ", pkg, ": ", conditionMessage(e)))
  }
}

for (p in pkgs_bioc) install_if_missing(p, "bioc")
for (p in pkgs_cran) install_if_missing(p, "cran")

# --- Session info for reproducibility ---
dir.create("logs", showWarnings = FALSE)
sink("logs/session_info.txt")
tryCatch({
  cat("=== SESSION INFO ===\n")
  cat("Date:", format(Sys.time()), "\n")
  cat("R:", R.version.string, "\n")
  cat("Bioconductor:", as.character(BiocManager::version()), "\n\n")
  sessionInfo()
}, finally = { sink() })

message("\n=== Key package versions ===")
for (p in c(pkgs_bioc, pkgs_cran)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
  message(sprintf("  %-22s %s", p, v))
}

# --- Final verification ---
missing <- c(pkgs_bioc, pkgs_cran)[
  !sapply(c(pkgs_bioc, pkgs_cran), requireNamespace, quietly = TRUE)
]
if (length(missing) == 0) {
  message("\n[OK] All packages installed successfully")
  message("[OK] session_info saved to logs/session_info.txt")
} else {
  message("\n[MISSING] Packages not installed: ", paste(missing, collapse = ", "))
  message("[WARNING] Scripts using these packages will fail")
}
