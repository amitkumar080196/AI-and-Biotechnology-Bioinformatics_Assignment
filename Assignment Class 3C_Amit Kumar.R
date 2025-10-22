#-----------------------------------------------
# Differential Expression Analysis: GSE15852
# Breast Cancer vs Normal Tissue
#-----------------------------------------------

# 1️⃣ Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs_needed <- c("GEOquery", "Biobase", "limma", "dplyr", "tibble",
                 "AnnotationDbi", "hgu133plus2.db", "ggplot2", "pheatmap")
for (p in pkgs_needed) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# 2️⃣ Download GSE15852 from GEO
gset_list <- getGEO("GSE15852", GSEMatrix = TRUE)
gset <- gset_list[[1]]  # Typically only one ExpressionSet

# 3️⃣ Extract expression matrix and sample info
expr_matrix <- exprs(gset)
sample_info <- pData(gset)
feature_info <- fData(gset)

# 4️⃣ Define probe IDs
probe_ids <- rownames(expr_matrix)

# 5️⃣ Annotate probes using hgu133plus2.db
annotated_df <- AnnotationDbi::select(hgu133plus2.db,
                                      keys = probe_ids,
                                      columns = c("SYMBOL", "GENENAME"),
                                      keytype = "PROBEID")

# 6️⃣ Merge annotation with expression data
expr_annotated <- data.frame(PROBEID = probe_ids,
                             expr_matrix,
                             check.names = FALSE)
expr_annotated <- merge(annotated_df, expr_annotated, by = "PROBEID", all.x = FALSE)

# 7️⃣ Handle duplicate probes: average expression per gene
expr_annotated <- expr_annotated[!is.na(expr_annotated$SYMBOL) & expr_annotated$SYMBOL != "", ]
expr_collapsed <- expr_annotated %>%
  dplyr::select(-PROBEID, -GENENAME) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("SYMBOL")

# Check dimensions
dim(expr_collapsed)

# 8️⃣ Define group factor (cancer vs normal)
table(sample_info$source_name_ch1)  # Check exact names
group <- factor(ifelse(grepl("normal", sample_info$source_name_ch1, ignore.case = TRUE),
                       "Normal", "Cancer"))

# 9️⃣ Design matrix for limma
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 10️⃣ Fit linear model
fit <- lmFit(expr_collapsed, design)

# 11️⃣ Define contrast: Cancer vs Normal
contrast_matrix <- makeContrasts(Cancer_vs_Normal = Cancer - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 12️⃣ Extract DEG results
deg_results <- topTable(fit2, number = Inf, adjust.method = "fdr")
head(deg_results)

# 13️⃣ Volcano plot
logFC_cutoff <- 1
adjP_cutoff <- 0.05

deg_results$regulation <- "Not Significant"
deg_results$regulation[deg_results$logFC >= logFC_cutoff & deg_results$adj.P.Val < adjP_cutoff] <- "Up"
deg_results$regulation[deg_results$logFC <= -logFC_cutoff & deg_results$adj.P.Val < adjP_cutoff] <- "Down"

volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Cancer vs Normal", x = "log2 Fold Change", y = "-log10 Adjusted P-value")

print(volcano_plot)

# Save volcano plot
if(!dir.exists("Results")) dir.create("Results")
ggsave("Results/volcano_plot.png", volcano_plot, width = 6, height = 5)

# 14️⃣ Save DEG results as CSV
write.csv(deg_results, "Results/deg_all.csv", row.names = TRUE)
write.csv(deg_results[deg_results$regulation == "Up", ], "Results/deg_upregulated.csv", row.names = TRUE)
write.csv(deg_results[deg_results$regulation == "Down", ], "Results/deg_downregulated.csv", row.names = TRUE)

# Use expr_collapsed directly
deg_results_valid <- deg_results[deg_results$SYMBOL %in% rownames(expr_collapsed) & !is.na(deg_results$SYMBOL), ]
top25_genes <- head(deg_results_valid$SYMBOL[order(deg_results_valid$adj.P.Val)], 25)
expr_top25 <- expr_collapsed[top25_genes, ]

library(pheatmap)

# Use top 25 genes from the collapsed matrix
top_genes <- head(rownames(expr_collapsed), 25)
expr_top25 <- expr_collapsed[top_genes, ]

# Check that the matrix is valid
if(nrow(expr_top25) > 0 && ncol(expr_top25) > 0 && any(apply(expr_top25, 1, sd) != 0)){
  pheatmap(expr_top25,
           scale = "row",
           annotation_col = data.frame(Group = group),
           main = "Heatmap of Top 25 DEGs",
           filename = "Results/top25_DEG_heatmap.png")
} else {
  message("Cannot generate heatmap: matrix has zero rows/columns or all rows are constant.")
}

