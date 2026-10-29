# Extended methods - the "why" behind each step

This document complements the Methods section of the manuscript. It
explains the rationale behind each analytical layer, what each step can
and cannot prove, and how the four layers converge into one argument.

## The four layers

| Layer | Folder | Question it answers |
|---|---|---|
| L1 - Bulk | analysis/bulk/ | Is there a signal in the tissue? |
| L2 - Cellular | analysis/cellular/ | Which cell type contributes? |
| L3 - snRNA-seq | analysis/snrnaseq/ | Which individual cells, with what phenotype? |
| L4 - Convergence | analysis/convergence/ | Do the layers agree? |

## Why three layers (and not one)

Bulk transcriptomics integrates averaged signals from millions of cells
and cannot distinguish whether a change comes from hepatocytes, Kupffer
cells, astrocytes, microglia, neurons, oligodendrocytes or infiltrating
lymphocytes. The cellular layer resolves that ambiguity. snRNA-seq
resolves it further: which individual cells, with what phenotype.

The hierarchy is deliberate: bulk generates the hypothesis; the cellular
layer tests it; snRNA-seq resolves it to single-cell granularity.

## What convergence proves (and does not prove)

Proves: the pattern is consistent across models (Mdr2-KO, EAE, DSS,
LPS), temporally coherent (liver before CNS), and cell-type consistent
(astrocytes SMS-dominant, MOL3_EAE B4galt6+).

Does not prove: causality. The step liver -> CNS ceramide transport is
not measured in this study - it is the central prediction of Aim 1.

## Falsifiability

The five predictions (P1-P5) are listed in Table 6 of the manuscript.
Each one is testable with a defined experiment. If P1-P2 fail, the
model falls.

## Reproducibility

- All scripts run under R 4.3.3 / Bioconductor 3.18 (see analysis/00_setup.R).
- Random seeds fixed to 42 for fgsea and any bootstrap.
- sessionInfo() at the end of each log.
