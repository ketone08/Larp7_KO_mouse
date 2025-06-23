# 加载必要的包
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)

seurat_obj <- readRDS("GSE104323_scRNA.rds")
target_ages <- c("P5", "P0", "E16.5")
cells_to_keep <- sc_data@meta.data$characteristics..age %in% target_ages
sc_data <- subset(sc_data, cells = rownames(sc_data@meta.data)[cells_to_keep])
dat_nIPC <- subset(sc_data, subset = characteristics..cell.cluster == "nIPC")

plot_data <- FetchData(
  dat_nIPC,
  vars = c("characteristics..age", "Hexim1", "Hexim2")
)

plot_data_long <- plot_data %>%
  pivot_longer(
    cols = c(Hexim1, Hexim2),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data_long$characteristics..age <- factor(
  plot_data_long$characteristics..age,
  levels = c( "E16.5", "P0","P5")
)

plot_data_long$combined_group <- paste(
  plot_data_long$Gene,
  plot_data_long$characteristics..age,
  sep = "_"
)
group_levels <- c(
  "Hexim1_E16.5", "Hexim1_P0","Hexim1_P5",
  "blank",  # 空白位置
  "Hexim2_E16.5", "Hexim2_P0","Hexim2_P5"
)

blank_data <- data.frame(
  orig.ident = factor(NA, levels = levels(plot_data_long$characteristics..age)),
  Gene = "blank",
  Expression = NA,
  combined_group = "blank"
)

final_plot_data <- bind_rows(plot_data_long, blank_data)
final_plot_data$combined_group <- factor(
  final_plot_data$combined_group,
  levels = group_levels
)

sample_colors <- c(
  "Hexim2" = rgb(123, 175, 222, maxColorValue = 255),
  "Hexim1" = rgb(251, 128, 114, maxColorValue = 255)
)

p <- ggplot(final_plot_data, aes(x = combined_group, y = Expression)) +
  geom_violin(
    aes(fill = Gene),
    scale = "width",
    trim = TRUE,
    na.rm = TRUE
  ) +
  scale_x_discrete(
    labels = c("E16.5", "P0","P5", "", "E16.5", "P0","P5")
  ) +
  labs(
    x = "",
    y = "Expression Level",
    title = "Hexim1 and Hexim2 Expression in nIPC Cells"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  annotate(
    "text",
    x = 2, y = max(final_plot_data$Expression, na.rm = TRUE) * 1.1,
    label = "Hexim1", size = 5
  ) +
  annotate(
    "text",
    x = 6, y = max(final_plot_data$Expression, na.rm = TRUE) * 1.1,
    label = "Hexim2", size = 5
  )

ggsave(
  filename = "nIPC_Hexim_expression of GSE104323_3.pdf",
  plot = p,
  width = 8,
  height = 6,
  device = "pdf"
)