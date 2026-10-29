# Analysis scripts — indexed by evidence layer

Numbering restarts at 01 in each folder. Order within a folder is
chronological (initial pipeline -> verification).

## L1 - Bulk (bulk/)

| # | Script | Purpose |
|---|---|---|
| 01 | deseq2_kegg_pipeline.R | DESeq2 + fgsea for recount3 datasets |
| 02 | consolidate_gsea.R | Merge pathway NES |
| 03 | verify_v8_numbers.R | Cross-check v8 numbers |
| 04 | verify_liver_numbers.R | Recompute liver stats |
| 05 | verify_liver_EAE_GSE312673.R | Liver EAE panel |
| 06 | DESeq2_liver_EAE_GSE312673.R | Full DESeq2 liver EAE |
| 07 | fix_comparison.R | Mdr2-KO vs EAE comparison |

## L2 - Cellular (cellular/)

| # | Script | Purpose |
|---|---|---|
| 01 | GSE139531_astrocytes.R | Astrocytes EAE (Chao 2019) |
| 02 | GSE139531_RPKM_by_group.R | Astrocyte RPKM by group |
| 03 | GSE253318_microglia.R | Microglia EAE |
| 04 | GSE279707_neurons.R | NeuN-sorted nuclei EAE |
| 05 | cellular_comparison.R | Cross-cell ratios (SMPD3/SGMS1, S1PR1/S1PR3) |
| 06 | cellular_heatmap_publication.R | Produces Figure F3 |
| 07 | check_composition_GSE279707.R | QC: composition confound |

## L3 - snRNA-seq (snrnaseq/)

| # | Script | Purpose |
|---|---|---|
| 01 | oligodendrocytes_GSE250589_annotation.R | Download + annotate |
| 02 | GSE250589_OLG_scores.R | Pseudobulk + Falcao table |
| 03 | falcao_quantification.R | PI, VI, entropy |
| 04 | falcao_doble_jab_table.R | Double-hit classification |
| 05 | diagnose_h5_and_barcodes.R | h5 diagnostics |
| 06 | test_acquisition_clean.R | MOL3 acquisition test |
| 07 | positivity_analysis.R | Percent positive cells + damage scores |

## L4 - Convergence (convergence/)

| # | Script | Purpose |
|---|---|---|
| 01 | spla2_hepatic_panel.R | sPLA2 panel in liver EAE |

## Figures (`figures/`)

| Source script | Output PDF |
|---|---|
| `00_figure_style.R` | Shared `theme_pub()` + AAAS palette |
| `F1_causal_chain.sh` + `F1_causal_chain.tikz` | `figures/main/F1_causal_chain.pdf` |
| `F2_bulk_heatmap.R` | `figures/main/F2_bulk_heatmap.pdf` |
| (via `cellular/06_cellular_heatmap_publication.R`) | `figures/main/F3_cellular_heatmap.pdf` |
| (via `snrnaseq/03_falcao_quantification.R`) | `figures/main/F4_snrnaseq_topology.pdf` |
| (via `snrnaseq/07_positivity_analysis.R`) | `figures/main/F5_mol3_conversion.pdf` |
| `F5_convergence_panel.R` | `figures/main/F6_convergence_panel.pdf` |
| `F4_snrnaseq_panels.sh` + `.tex` | (deprecated — panels split into F4 + F5) |

## Usage

    Rscript analysis/00_setup.R
    Rscript analysis/figures/F2_bulk_heatmap.R
    Rscript analysis/figures/F5_convergence_panel.R
    bash    analysis/figures/F1_causal_chain.sh
