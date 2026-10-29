# Consolidar TODOS los GSEA de resultados
library(dplyr)

files <- list.files("results", pattern = "_GSEA_.*\\.csv$", full.names = TRUE)
message("Archivos GSEA encontrados: ", length(files))

rutas_hipotesis <- c(
  "Bile secretion",
  "Primary bile acid biosynthesis",
  "PPAR signaling pathway",
  "Fatty acid degradation",
  "Fatty acid metabolism",
  "Sphingolipid metabolism",
  "Sphingolipid signaling pathway",
  "NF-kappa B signaling pathway",
  "TNF signaling pathway",
  "Chemokine signaling pathway",
  "Inflammatory bowel disease"
)

resultados <- list()
for (f in files) {
  # Nombre del archivo: SRPXXXXXX_GSEA_TEST_vs_REF.csv
  filename <- basename(f)
  proj_test_ref <- sub("_GSEA_", "|", sub("\\.csv$", "", filename))
  parts <- strsplit(proj_test_ref, "\\|")[[1]]
  proj <- parts[1]
  comparison <- parts[2]
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  for (ruta in rutas_hipotesis) {
    match_row <- df[grepl(ruta, df$pathway_name, ignore.case = TRUE), ]
    if (nrow(match_row) > 0) {
      resultados[[length(resultados) + 1]] <- data.frame(
        project = proj,
        comparison = comparison,
        pathway = ruta,
        NES = round(match_row$NES[1], 2),
        padj = match_row$padj[1],
        sig = ifelse(match_row$padj[1] < 0.05, "✓", "ns")
      )
    }
  }
}

consolidated <- bind_rows(resultados)
write.csv(consolidated, "results/00_CONSOLIDATED_all_hypothesis_pathways.csv",
          row.names = FALSE)

# Mostrar tabla pivote
pivot <- consolidated %>%
  tidyr::pivot_wider(
    id_cols = pathway,
    names_from = c(project, comparison),
    values_from = NES
  )
print(pivot, width = Inf)
