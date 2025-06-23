library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)
# read data
merge <- readRDS("merge.rds")
meta_data <- merge@meta.data

merge$celltype_merged <- merge$celltype
merge$celltype_merged[merge$celltype %in% c("Neuroblast1", "Neuroblast2")] <- "Neuroblast"

# selected_celltypes
selected_celltypes <- c("Neuroblast", "nIPC", "RGL")
subset_obj <- subset(merge, celltype_merged %in% selected_celltypes)

features <- c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score")
#vilon plot
violin_plots <- lapply(features, function(feature) {
  VlnPlot(subset_obj, features = feature, group.by = "celltype_merged",split.by = "orig.ident",pt.size = 0) +
    ggtitle(paste("Violin plot of", feature)) 
})

# combine
combined_plot <- plot_grid(plotlist = violin_plots, ncol = 1)

# save PDF
ggsave("P7_Neuroblast_violin_plots_combined.pdf", plot = combined_plot, width = 10, height = 18)







