
# packages
library(Seurat)
library(ggplot2)
library(dplyr)

dat<-readRDS("03_merge_mannual_annotation_sub.RDS")

genes_of_interest <- c("Ccne2","Cdc6","Cdca7","Cdc20","Cdca3","Cdca8","Cdk1","Cenpa","Cenpe","Cenpf","Ckap2","Ckap2l","Ckap5","E2f8","Mcm2","Mcm4",	"Mcm5","Mcm6","Msh2","Mki67","Pcna","Pola1","Pold3","Top2a")

# selected nIPC 
nIPC_cells <- subset(dat, idents = "nIPC")
expr_data <- GetAssayData(nIPC_cells, slot = "data")
metadata <- nIPC_cells@meta.data

results <- data.frame()

# Average expression
for (gene in genes_of_interest) {
  gene_expr <- expr_data[gene, ]
  gene_metadata <- data.frame(cell = names(gene_expr), expression = as.numeric(gene_expr), orig.ident = metadata$orig.ident)
 
  summary_stats <- gene_metadata %>%
    group_by(orig.ident) %>%
    summarize(
      avg_expression = mean(expression, na.rm = TRUE),
      total_cells = n(),
      expressed_cells = sum(expression > 0, na.rm = TRUE),
      .groups = 'drop'
    )
  
  summary_stats$expression_ratio <- (summary_stats$expressed_cells / summary_stats$total_cells) * 100
  
  # data.frame
  for (i in 1:nrow(summary_stats)) {
    results <- rbind(results, data.frame(Gene = gene, Sample = summary_stats$orig.ident[i], 
                                          AvgExpression = summary_stats$avg_expression[i], 
                                          ExpressionRatio = summary_stats$expression_ratio[i]))
  }
}

# sample comparison
samples_to_plot <- list(c("p7_l7_wt", "p7_l7_nestin"), 
                        c("p7_l7_wt", "p7_l7_emxt"), 
                        c("p7_l7_wt", "p7_l7_nestin", "p7_l7_emxt"))


# save PNG
for (samples in samples_to_plot) {
  plot_data <- results %>% filter(Sample %in% samples)
 
  plot_data$Gene <- factor(plot_data$Gene, levels = genes_of_interest)

  plot_data$Sample <- factor(plot_data$Sample, levels = c("p7_l7_wt", setdiff(samples, "p7_l7_wt")))  
  p <- ggplot(plot_data, aes(x = Sample, y = Gene)) +
    geom_point(aes(size =ExpressionRatio , color = AvgExpression), alpha = 0.6) +
    scale_size_continuous(range = c(2, 10), name = "Expression Ratio (%)") +  
    scale_color_gradient(low = "#330066", high = "#FFCC33", name = "Average Expression") + 
    labs(title = paste("Gene Expression in nIPC Cells:", paste(samples, collapse = " vs ")), 
         x = "Sample", y = "Gene") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(filename = paste0("nIPC_bubble_plot_selected_genes_", paste(samples, collapse = "_vs_"), ".pdf"), plot = p, width = 7, height = 6)
}