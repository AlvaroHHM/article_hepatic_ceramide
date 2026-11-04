#
# Requires:
#   data/GSE250589/h5/*.h5 (28 files, ~2.3 GB)
#   Download: bash analysis/snrnaseq/00_download_h5.sh
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(hdf5r); library(ggplot2)
})
sink("logs/53c_log.txt", split=TRUE); on.exit(sink())
cat("=== 53c: Re-extraer matriz + PCA + ARI ===\n")

# 1. MOL3
annot <- fread("data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz",
               data.table=FALSE, check.names=FALSE)
rownames(annot) <- annot[[1]]; annot <- annot[,-1]
mol3 <- annot[annot$cellTypeFalcao2018 == "MOL3_EAE",]
mol3$barcode_clean <- sub("-.*$","",rownames(mol3))
mol3$zheng <- ifelse(mol3$cellType=="MOL2", "MOL2_like",
              ifelse(mol3$cellType=="MOL56","MOL56_like","Other"))

# 2. Lista de genes (mouse case)
res <- read.csv("results/52_DESeq2_MOL2_vs_MOL56.csv", check.names=FALSE)
mouse_case <- function(x) paste0(toupper(substr(x,1,1)), tolower(substr(x,2,nchar(x))))
top_up <- head(res$gene[!is.na(res$log2FoldChange) & res$log2FoldChange>0.5 &
                        !is.na(res$padj) & res$padj<0.05], 20)
top_dn <- head(res$gene[!is.na(res$log2FoldChange) & res$log2FoldChange< -0.5 &
                        !is.na(res$padj) & res$padj<0.05], 20)
top_genes <- unique(c("Klk6","Ptgds","B4galt6",
                      mouse_case(top_up), mouse_case(top_dn)))

# 3. Re-extraer
h5_files <- list.files("data/GSE250589/h5", pattern="\\.h5$", full.names=TRUE)
mapping <- fread("results/34_GSM_mapping_v2.csv", data.table=FALSE)
expr_list <- list()
for (h5 in h5_files) {
  gsm <- sub("_.*","",basename(h5))
  samp <- mapping$orig.ident[mapping$GSM==gsm]
  if (length(samp)!=1) next
  bcs_samp <- mol3$barcode_clean[mol3$orig.ident==samp]
  if (length(bcs_samp)==0) next
  f <- H5File$new(h5, mode="r")
  shape <- as.integer(f[["matrix/shape"]][]); idx <- as.integer(f[["matrix/indices"]][])
  ptr <- as.integer(f[["matrix/indptr"]][]); dat <- as.numeric(f[["matrix/data"]][])
  feat <- sub("\\..*$","",f[["matrix/features/name"]][])
  bcs <- sub("-.*$","",f[["matrix/barcodes"]][])
  f$close_all()
  if (max(idx) >= shape[1]) shape <- rev(shape)
  mat <- sparseMatrix(i=idx+1L, p=ptr, x=dat, dims=shape, index1=TRUE)
  rownames(mat) <- feat; colnames(mat) <- bcs
  cells <- intersect(bcs_samp, colnames(mat))
  if (length(cells) < 5) next
  genes_here <- intersect(top_genes, rownames(mat))
  if (length(genes_here)==0) next
  sub_mat <- mat[genes_here, cells, drop=FALSE]
  for (g in genes_here) {
    expr_list[[g]] <- c(expr_list[[g]],
                         setNames(sub_mat[g,], paste(samp, cells, sep="_")))
  }
  cat(sprintf("  [OK] %s\n", samp))
}

common_cells <- Reduce(intersect, lapply(expr_list, names))
expr_mat <- do.call(rbind, lapply(expr_list, function(x) x[common_cells]))
cat(sprintf("\nMatriz: %d genes × %d células\n", nrow(expr_mat), ncol(expr_mat)))

mol3_key <- paste(mol3$orig.ident, mol3$barcode_clean, sep="_")
mol3_sub <- mol3[match(colnames(expr_mat), mol3_key),]

# Guardar
write.csv(expr_mat, "results/53c_expr_mat.csv")
write.csv(mol3_sub[,c("orig.ident","barcode_clean","zheng")],
          "results/53c_labels.csv", row.names=FALSE)
cat("[OK] 53c_expr_mat.csv + 53c_labels.csv\n\n")

# 4. PCA
keep <- apply(expr_mat, 1, sd) > 0
expr_filt <- expr_mat[keep,,drop=FALSE]
cat(sprintf("Genes con varianza: %d\n", nrow(expr_filt)))

nCount <- mol3_sub$nCount_RNA
expr_norm <- log1p(t(t(expr_filt) / nCount * 1e4))

set.seed(42)
pca <- prcomp(t(expr_norm), center=TRUE, scale.=FALSE)
var_exp <- summary(pca)$importance[2,1:10]
cat(sprintf("PC1-10 var explicada: %.1f%%\n\n", 100*sum(var_exp)))

# 5. ARI
ari_fn <- function(a,b) {
  tab <- table(a,b); n <- sum(tab)
  sc <- function(x) sum(choose(x,2))
  idx <- sc(tab); ei <- sc(rowSums(tab)); ej <- sc(colSums(tab))
  if (choose(n,2)==0) return(0)
  (idx - ei*ej/choose(n,2)) / ((ei+ej)/2 - ei*ej/choose(n,2))
}
valid <- mol3_sub$zheng != "Other"
zheng_bin <- ifelse(mol3_sub$zheng=="MOL2_like", 1, 2)

cat("--- ARI vs Zheng para k-means k=2..6 sobre PC1-10 ---\n")
for (k in 2:6) {
  set.seed(42)
  km <- kmeans(pca$x[,1:10], centers=k)
  ari <- ari_fn(km$cluster[valid], zheng_bin[valid])
  cat(sprintf("  k=%d: ARI=%.3f\n", k, ari))
}

# k=2 detalle
set.seed(42)
km2 <- kmeans(pca$x[,1:10], centers=2)
ari_k2 <- ari_fn(km2$cluster[valid], zheng_bin[valid])
cat(sprintf("\n=== k=2 detalle ===\n"))
cat(sprintf("ARI PCA(10PC) vs Zheng: %.3f\n", ari_k2))
tab <- table(km2$cluster, mol3_sub$zheng)
print(tab)
cat(sprintf("\nCluster 1: %d cells | Cluster 2: %d cells\n",
            sum(km2$cluster==1), sum(km2$cluster==2)))
cat(sprintf("Concordancia si cluster 1=MOL2: %.1f%%\n",
            100*mean(mol3_sub$zheng[km2$cluster==1]=="MOL2_like", na.rm=TRUE)))
cat(sprintf("Concordancia si cluster 2=MOL56: %.1f%%\n",
            100*mean(mol3_sub$zheng[km2$cluster==2]=="MOL56_like", na.rm=TRUE)))

# 6. Figura
df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],
                 zheng=mol3_sub$zheng, cluster=factor(km2$cluster))
p <- ggplot(df[df$zheng!="Other",], aes(PC1, PC2, color=zheng)) +
  geom_point(alpha=0.2, size=0.5) +
  scale_color_manual(values=c(MOL2_like="#B2182B", MOL56_like="#2166AC")) +
  theme_minimal() +
  labs(title=sprintf("PCA 42 genes MOL3 | ARI vs Zheng = %.3f", ari_k2),
       subtitle=sprintf("%d células MOL2_like + %d MOL56_like",
                        sum(mol3_sub$zheng=="MOL2_like"),
                        sum(mol3_sub$zheng=="MOL56_like")))
ggsave("results/53c_PCA_MOL3.pdf", p, width=8, height=6)
cat("[OK] 53c_PCA_MOL3.pdf\n")

cat("\n=== FIN 53c ===\n")
