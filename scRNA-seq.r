#mouse ID8_P53KO_NACT
#创建seurat对象

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(homologene)


after<-Read10X("./ID8_P53KO_NACT_after/filtered_feature_bc_matrix/")
after <- CreateSeuratObject(counts =after , project = "ID8_P53KO_NACT_after",min.cells = 3, min.features = 200)

#质量控制
after[["percent.mt"]] <- PercentageFeatureSet(after,pattern = "^mt-")
preQC_after <- VlnPlot(after, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        ncol = 3, 
                        group.by = "orig.ident", 
                        pt.size = 0)
ggsave("./after_preQC_1.pdf", preQC_after, width=7 ,height=6)

#过滤低质量数据
after <- subset(after, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)	
preQC_after <- VlnPlot(after, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        ncol = 3, 
                        group.by = "orig.ident", 
                        pt.size = 0)
						
after[["percent.rb"]] <- PercentageFeatureSet(after, pattern = "^RP[SL]")



#去除线粒体基因
mt.index <- grep(pattern = "^mt-", x = rownames(x = after@assays$RNA@counts), value = FALSE)
after <- after[-mt.index,]


#scale数据
after <- NormalizeData(after,normalization.method = "LogNormalize", scale.factor = 10000)

#特征选择
after <- FindVariableFeatures(after, nfeatures = 4000, selected.method="vst")

# 可视化高变基因
var_plot <- VariableFeaturePlot(after)
p3<-LabelPoints(plot = var_plot, points = head(VariableFeatures(after), repel = TRUE))
ggsave("after_Variable_Genes.pdf", p3,width = 10)

#缩放数据
after <- ScaleData(after,verbose=FALSE)


#细胞周期评分

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
mouse.s.genes<-homologene(s.genes, inTax = 9606, outTax = 10090)
mouse.g2m.genes<-homologene(g2m.genes, inTax = 9606, outTax = 10090)

mouse.s.gene <- as.character(mouse.s.genes[[2]])
mouse.g2m.gene <- as.character(mouse.g2m.genes[[2]])
after <- CellCycleScoring(after, 
                           s.features = mouse.s.gene, 
                           g2m.features = mouse.g2m.gene, 
                           set.ident = TRUE)

# PCA降维
after <- RunPCA(
  after,
  features = VariableFeatures(object = after),
  npcs = 50  # 计算前50个主成分
)

# 可视化PCA结果
p1<-DimPlot(after, reduction = "pca")
ggsave("./after_pca_1.pdf",p1, width=7 ,height=6)
p2<-ElbowPlot(after, ndims = 50)  # 选择主成分数量（拐点位置）
ggsave("./after_Eblowplot_1.pdf", p2, width=7 ,height=6)

#聚类分析
# 确定主成分数（通常10-30）
pca_dims <- 1:30  # 根据ElbowPlot调整

# 构建KNN图
after <- FindNeighbors(after, dims = pca_dims)

# 聚类（分辨率越高，簇越多）
after <- FindClusters(
  after,
  resolution = 0.7  # 通常0.4-1.2（小样本用低值）
)

# 非线性降维 (UMAP/t-SNE)
after <- RunUMAP(after, dims = pca_dims)
after <- RunTSNE(after, dims = pca_dims)

# 可视化聚类结果
umap_plot <- DimPlot(after, reduction = "umap", label = TRUE)
tsne_plot <- DimPlot(after, reduction = "tsne", label = TRUE)
combined_plot <- umap_plot + tsne_plot
ggsave("after_Clustering_UMAP_TSNE.pdf", combined_plot, width = 14)

# ----------------------
# 步骤7: 差异表达基因（标记基因）
# ----------------------
# 寻找每个簇的标记基因
markers <- FindAllMarkers(
  after,
  only.pos = TRUE,     # 只保留阳性标记
  min.pct = 0.25,      # 基因在至少25%的簇细胞中表达
  logfc.threshold = 0.25
)
write.table(markers,file="./after_significant_marker_genes.txt",quote=F,sep="\t",col.names=FALSE)
# 提取每个簇的前5标记基因
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# 可视化标记基因表达
heatmap_plot <- DoHeatmap(
  after, 
  features = top_markers$gene
)
ggsave("after_Marker_Heatmap.pdf", heatmap_plot, width = 12)

saveRDS(after, file = "./01_ID8_p53KO_NACT_after_del_mt.RDS")


#merge
sampleList <- list(before,after)


merge <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
merge <- IntegrateData(anchorset = merge, dims = 1:50)
saveRDS(merge, file = "./merge/01_merge_mt_del.RDS")

DefaultAssay(merge) <- "integrated"


merge <- SCTransform(merge, variable.features.n=4000)
merge <- RunPCA(merge, npcs = 50, verbose = FALSE)
merge <- RunHarmony(merge,group.by.vars='orig.ident',dims=1:50,assay='SCT')
merge <- FindNeighbors(merge, reduction = "harmony", dims = 1:30)
merge <- FindClusters(merge, resolution = 0.5,algorithm=1)
merge <- RunUMAP(merge, reduction = "harmony", dims = 1:30,assay='SCT')

p1 <- DimPlot(merge, reduction = "umap", label = TRUE, repel = TRUE)
p2 <- DimPlot(merge, reduction = "umap", group.by = "orig.ident",label=TRUE)
write.table(merge$seurat_clusters,file="./all_merge_umapCluster.txt",quote=F,sep="\t",col.names=FALSE)
ggsave("./01_merge_umap.pdf", p1, width=12 ,height=10)
ggsave("./01_merge_umap_group.pdf", p2, width=12 ,height=10)
saveRDS(merge, file = "./01_merge_mt_del_processed.RDS")

#Cell Annotation
marker_genes <- list(
  Astrocytes = c("Sox9", "Etv4","Sall3","Grin2c","Kcng4"),
  RG_like=c("Tfap2c", "Hopx","Rhcg","Vnn1","Wnt8b"),
  nIPC=c("Mxd3","Ascl1","Neurog2","Rfc4","Lockd"),
  Neuroblast1=c("Eomes","Neurod4","Slc17a6","Mpped1","Calb2"),
Neuroblast2=c("Sox11","Bhlhe22","Neurod2","Gal","Igfbpl1"))
merge <- AddModuleScore(merge, features = marker_genes, name = "ModuleScore")
cluster_annotations <- merge@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(across(starts_with("ModuleScore"), mean))

predicted_cell_types <- cluster_annotations %>%
  rowwise() %>%
  mutate(predicted_type = names(marker_genes)[which.max(c_across(starts_with("ModuleScore")))])

# 将预测的类型加入到Seurat对象的元数据中
merge@meta.data$predicted_type <- factor(merge@meta.data$seurat_clusters, 
                                              levels = cluster_annotations$seurat_clusters,
                                              labels = predicted_cell_types$predicted_type)
											  
p1 <- DimPlot(merge, reduction = "umap", label = TRUE, repel = TRUE)
p2 <- DimPlot(merge, reduction = "umap", group.by = "predicted_type",label=TRUE)
DimPlot(merge, group.by = "predicted_type", label = TRUE, repel = TRUE) + NoLegend()
write.table(merge$seurat_clusters,file="merge_umapCluster.txt",quote=F,sep="\t",col.names=FALSE)
ggsave("merge_umap_predict_celltype.pdf", p1, width=12 ,height=10)
ggsave("merge_umap_predict_celltype_group.pdf", p2, width=12 ,height=10)


#执行差异表达分析

# 获取所有细胞群体
clusters <- unique(merge$manual_annotation)

# 定义组别
comparisons <- list(
  c("p7_l7_nestin", "p7_l7_wt"),
  c("p7_l7_nestin", "p7_l7_emxt"),
  c( "p7_l7_emxt", "p7_l7_wt")
)

# 存储结果
results_list <- list()

for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  for (cluster in clusters) {
    # 选择当前群体的细胞
    cluster_cells <- WhichCells(merge, expression = manual_annotation == cluster)
    
    # 进行差异表达分析
    degs <- FindMarkers(
      object = merge,
      ident.1 = group1,
      ident.2 = group2,
      cells = cluster_cells,
      group.by = "orig.ident", min.pct=0.1,logfc.threshold=0.25, only.pos=F
    )
    
    # 存储结果，命名格式为: "群体_比较组1_vs_比较组2"
    results_list[[paste(cluster, group1, "vs", group2, sep = "_")]] <- degs
  }
}
保存结果为CSV文件

将每个比较的结果保存为CSV文件。

for (name in names(results_list)) {
  output_file <- paste0(name, "_DEGs.csv")
  write.csv(results_list[[name]], file = output_file, row.names = TRUE)
}


