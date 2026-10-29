# Data access

No raw data is stored in this repository (see .gitignore). All
datasets are public and can be downloaded with the commands below.

## recount3 (bulk)

    library(recount3)
    for (srp in c("SRP074136","SRP148279","SRP197582",
                  "SRP071304","SRP161558")) {
      rse <- create_rse_manual("mouse", srp)
      saveRDS(rse, sprintf("data/raw/%s.rds", srp))
    }

## GEO (bulk + cellular + snRNA-seq)

    library(GEOquery)
    for (gse in c("GSE312673","GSE224127","GSE139531",
                  "GSE253318","GSE279707","GSE250589")) {
      getGEO(gse, destdir = "data/raw", GSEMatrix = TRUE)
    }

## Lipidome

PR000633 - human liver lipidome (NAFLD). Available at
https://lipidomer.shinyapps.io/

## Layout

Once downloaded, place files under data/raw/<accession>/ (gitignored):

    data/raw/
      SRP074136/  SRP148279/  SRP197582/  SRP071304/  SRP161558/
      GSE312673/  GSE224127/  GSE139531/
      GSE253318/  GSE279707/  GSE250589/
