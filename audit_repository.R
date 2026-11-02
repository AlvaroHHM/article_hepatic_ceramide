# ============================================================
# audit_repository.R
# Auditoría en 4 fases del compendio article_hepatic_ceramide/
# Uso: Rscript audit_repository.R
# ============================================================

sep <- paste(rep("=", 60), collapse = "")
cat(sep, "\n")
cat("AUDITORÍA DEL REPOSITORIO article_hepatic_ceramide/\n")
cat("Fecha:", format(Sys.time()), "\n")
cat(sep, "\n\n")

# ---- Setup ----
repo <- getwd()
scripts <- list.files("analysis", pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
csvs <- list.files("results", pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
figs <- list.files("figures", pattern = "[.]pdf$", recursive = TRUE, full.names = TRUE)

cat("Scripts encontrados: ", length(scripts), "\n", sep = "")
cat("CSVs encontrados:    ", length(csvs), "\n", sep = "")
cat("Figuras encontradas: ", length(figs), "\n\n", sep = "")

# ============================================================
# FASE 1: SYNTAX CHECK
# ============================================================
cat("--- FASE 1: SYNTAX CHECK ---\n")
syntax_fail <- character(0)
for (s in scripts) {
  ok <- tryCatch({ parse(s); TRUE }, error = function(e) FALSE)
  if (!ok) {
    syntax_fail <- c(syntax_fail, s)
    cat(sprintf("  ❌ %s\n", s))
  }
}
if (length(syntax_fail) == 0) {
  cat(sprintf("  ✅ Todos los %d scripts parsean OK\n", length(scripts)))
} else {
  cat(sprintf("  ❌ %d/%d scripts con error de sintaxis\n",
              length(syntax_fail), length(scripts)))
}
cat("\n")

# ============================================================
# FASE 2: DEPENDENCIES
# ============================================================
cat("--- FASE 2: DEPENDENCIES ---\n")
# Extraer todos los library()/require() de todos los scripts
extract_libs <- function(path) {
  txt <- readLines(path, warn = FALSE)
  pat <- "(?:library|require)\\s*\\(\\s*['\"]?([A-Za-z0-9._]+)['\"]?\\s*\\)"
  hits <- unlist(regmatches(txt, gregexpr(pat, txt, perl = TRUE)))
  hits <- gsub(pat, "\\1", hits, perl = TRUE)
  unique(hits)
}
all_libs <- unique(unlist(lapply(scripts, extract_libs)))
all_libs <- sort(all_libs[all_libs != "NULL" & !is.na(all_libs) & nchar(all_libs) > 0])
cat(sprintf("  Paquetes únicos importados: %d\n", length(all_libs)))
missing <- character(0)
for (pkg in all_libs) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (ok) {
    cat(sprintf("  ✅ %-20s\n", pkg))
  } else {
    cat(sprintf("  ❌ %-20s MISSING\n", pkg))
    missing <- c(missing, pkg)
  }
}
cat("\n")
if (length(missing) > 0) {
  cat(sprintf("  ⚠️  %d paquete(s) faltante(s): %s\n",
              length(missing), paste(missing, collapse = ", ")))
  cat(sprintf("  Fix: Rscript -e 'install.packages(c(%s))'\n",
              paste(sprintf("\"%s\"", missing), collapse = ", ")))
} else {
  cat("  ✅ Todos los paquetes instalados\n")
}
cat("\n")

# ============================================================
# FASE 3: INPUTS / OUTPUTS
# ============================================================
cat("--- FASE 3: INPUTS / OUTPUTS ---\n")

# Función para extraer rutas de read*/write* de un script
extract_paths <- function(path, op) {
  txt <- readLines(path, warn = FALSE)
  pat <- paste0(op, "\\s*\\(\\s*['\"]([^'\"]+)['\"]")
  hits <- unlist(regmatches(txt, gregexpr(pat, txt, perl = TRUE)))
  if (length(hits) == 0) return(character(0))
  unique(gsub(pat, "\\1", hits, perl = TRUE))
}

# Inputs
all_inputs <- unique(unlist(lapply(scripts, function(s) extract_paths(s, "read\\.(csv|RDS|delim|table|tsv)"))))
all_inputs <- all_inputs[!grepl("^(https?|ftp)://", all_inputs)]  # excluir URLs
cat(sprintf("  Inputs referenciados: %d\n", length(all_inputs)))
missing_inputs <- all_inputs[!file.exists(all_inputs)]
if (length(missing_inputs) > 0) {
  cat(sprintf("  ❌ %d inputs faltantes:\n", length(missing_inputs)))
  for (m in head(missing_inputs, 10)) cat(sprintf("     - %s\n", m))
} else {
  cat("  ✅ Todos los inputs existen\n")
}

# Outputs
all_outputs <- unique(unlist(lapply(scripts, function(s) extract_paths(s, "write\\.(csv|table|delim)"))))
cat(sprintf("\n  Outputs referenciados: %d\n", length(all_outputs)))
missing_outputs <- all_outputs[!file.exists(all_outputs)]
empty_outputs <- character(0)
for (o in all_outputs[file.exists(all_outputs)]) {
  if (file.info(o)$size == 0) empty_outputs <- c(empty_outputs, o)
}
if (length(missing_outputs) > 0) {
  cat(sprintf("  ⚠️  %d outputs no existen todavía (pueden ser opcionales):\n", length(missing_outputs)))
  for (m in head(missing_outputs, 5)) cat(sprintf("     - %s\n", m))
} else {
  cat("  ✅ Todos los outputs existen\n")
}
if (length(empty_outputs) > 0) {
  cat(sprintf("  ❌ %d outputs vacíos (0 bytes):\n", length(empty_outputs)))
  for (e in head(empty_outputs, 5)) cat(sprintf("     - %s\n", e))
} else {
  cat("  ✅ Ningún output vacío\n")
}
cat("\n")

# ============================================================
# FASE 4: MANUSCRIPT NUMBER CROSS-CHECK
# ============================================================
cat("--- FASE 4: MANUSCRIPT NUMBER CROSS-CHECK ---\n")
cat("  Buscando 15 números clave del manuscrito en los CSVs...\n\n")

# Leer todos los CSVs una vez
cat("  Cargando", length(csvs), "CSVs en memoria...\n")
csv_content <- list()
for (f in csvs) {
  csv_content[[f]] <- tryCatch(
    readLines(f, warn = FALSE, n = 1000),
    error = function(e) character(0)
  )
}

# Números clave: valor → descripción
targets <- list(
  "4.80"  = "SMPD3 log2FC (Mdr2-KO)",
  "3.77"  = "B4GALT6 log2FC (Mdr2-KO)",
  "2.88"  = "Smpd3 log2FC (EAE pre-onset)",
  "1.77"  = "Smpd3 log2FC (EAE peak)",
  "1.49"  = "B4galt6 log2FC (EAE pre-onset)",
  "20.3"  = "S1PR1/S1PR3 ratio (astrocitos)",
  "21.1"  = "S1PR1/S1PR3 ratio (microglía)",
  "2.7"   = "S1PR1/S1PR3 ratio (neuronas)",
  "7.06"  = "B4galt6+ % (Early)",
  "17.9"  = "B4galt6+ % (Peak)",
  "59.6"  = "MOL3 % del compartimento EAE-OLG",
  "0.061" = "DamageScore MOL3∩MOL2",
  "0.006" = "DamageScore MOL3∩MOL56",
  "4.3e-06" = "LacCer module p-value",
  "0.40"  = "ARI MOL3 subdivisión"
)

found_count <- 0
for (val_str in names(targets)) {
  desc <- targets[[val_str]]
  # Buscar el valor como string (con tolerancia de formato)
  pattern <- paste0("(^|[^0-9.])", gsub("[.]", "[.]", val_str), "([^0-9]|$)")
  hits <- character(0)
  for (f in names(csv_content)) {
    if (any(grepl(pattern, csv_content[[f]]))) {
      hits <- c(hits, f)
    }
  }
  if (length(hits) > 0) {
    cat(sprintf("  ✅ %-10s (%s) → %d CSV: %s\n",
                val_str, desc, length(hits),
                paste(head(basename(hits), 2), collapse = ", ")))
    found_count <- found_count + 1
  } else {
    cat(sprintf("  ⚠️  %-10s (%s) → NO encontrado\n", val_str, desc))
  }
}

cat(sprintf("\n  Resumen: %d/%d números encontrados en CSVs\n",
            found_count, length(targets)))
cat("\n")

# ============================================================
# RESUMEN FINAL
# ============================================================
cat(sep, "\n")
cat("RESUMEN\n")
cat(sep, "\n")
cat(sprintf("  Scripts:              %d/%d syntax OK\n",
            length(scripts) - length(syntax_fail), length(scripts)))
cat(sprintf("  Paquetes faltantes:   %d\n", length(missing)))
cat(sprintf("  Inputs faltantes:     %d\n", length(missing_inputs)))
cat(sprintf("  Outputs vacíos:       %d\n", length(empty_outputs)))
cat(sprintf("  Números verificados:  %d/%d\n", found_count, length(targets)))
cat("\n")

# Veredicto
issues <- length(syntax_fail) + length(missing) + length(empty_outputs)
if (issues == 0 && found_count >= length(targets) * 0.8) {
  cat("✅ VEREDICTO: Repositorio sano y verificable\n")
} else if (issues <= 3) {
  cat("⚠️  VEREDICTO: Repositorio con issues menores\n")
} else {
  cat("❌ VEREDICTO: Repositorio requiere atención\n")
}
cat(sep, "\n")
