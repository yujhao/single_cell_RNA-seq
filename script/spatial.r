# jhyu
# 2025-02-27
# spatial_transcriptomics
suppressWarnings({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(clusterProfiler)
  library(fgsea)
  library(enrichplot)
  library(pheatmap)
  library(igraph)
  library(ggraph)
  library(SPOTlight)
  library(igraph)
  library(RColorBrewer)
  library(spacexr)
  library(doParallel)
  library(glmGamPoi)
  library(tidyverse)
})

#1 10×空间转录组数据读取 
dir.create("1.SpaceRanger_QC")
brain = readRDS("./SpaceRanger/brain.rds")
brain = PercentageFeatureSet(brain, pattern = "^MT-",
                             col.name = "percent.mito")
Idents(brain) = "sampleid"
for (i in c("nCount_Spatial","nFeature_Spatial","percent.mito")){
p1 = VlnPlot(brain,
             features = i,
             pt.size = 0.1) + NoLegend()
p2 = SpatialPlot(brain,
                 ncol = 2,
                 features=i,
                 pt.size.factor = 1.4,
                 alpha = 1)
plot = wrap_plots(p1,p2)
plot = plot
ggsave(paste0("1.SpaceRanger_QC/",i,".pdf"),
       plot = plot,
       device = "pdf",
       width = 16, 
       dpi = 300)
ggsave(paste0("1.SpaceRanger_QC/",i,".png"),
       plot = plot,
       device = "png",
       width = 16, 
       dpi = 300)
}

#2 标准化与降维聚类
dir.create("2.Clustering",recursive = T)
brain = Seurat::SCTransform(
  brain,
  assay = "Spatial",
  method = "glmGamPoi", 
  vars.to.regress =  "percent.mito",
  variable.features.n = 3000,
  verbose = TRUE,
  return.only.var.genes = FALSE)
brain = RunPCA(brain, 
               assay = "SCT", 
               features = VariableFeatures(brain))
p1 = DimPlot(brain, 
             reduction = "pca",
             group.by = "orig.ident")
p2 = ElbowPlot(brain, ndims = 50, 
               reduction = "pca")
p3 = wrap_plots(p1,p2)
p3
ggsave("2.Clustering/RunPCA.pdf", 
       plot = p3,
       device = "pdf", 
       width = 10, 
       height = 4, 
       dpi = 300)
brain = FindNeighbors(brain, 
                      reduction = "pca", 
                      dims = 1:30, 
                      features = VariableFeatures(brain))
brain = FindClusters(brain, 
                     verbose = TRUE,resolution = 0.4)
brain = RunUMAP(brain, 
                reduction = "pca", 
                dims = 1:30)
brain@meta.data$clusters = as.numeric(brain@meta.data$seurat_clusters)
brain@meta.data$clusters = as.factor(brain@meta.data$clusters) 
Idents(brain) = "clusters"
p1 = DimPlot(brain, 
             reduction = "umap", 
             label = TRUE)
p2 = SpatialPlot(brain, 
                 ncol = 2,
                 group.by = "clusters",
                 pt.size.factor = 1.6 )
p3 = wrap_plots(p1,p2)
p3
ggsave("2.Clustering/UMAP.pdf", 
       plot = p3, 
       device = "pdf",
       width = 16, 
       height = 4, 
       dpi = 300)

#3 marker基因鉴定
dir.create("3.Marker")
All_markers = FindAllMarkers(brain, 
                             only.pos = TRUE,
                             test.use = "wilcox")
top10 = All_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(All_markers, 
            "3.Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
write.table(top10, 
            "3.Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
p1 = DoHeatmap(brain, 
               features = top10$gene,
               group.by="clusters",
               size = 2)
ggsave("3.Marker/Heatmap.pdf", 
       plot = p1,
       device = "pdf", 
       width = 10, 
       height = 15, 
       dpi = 300)

#4 差异基因鉴定
dir.create("4.Diffexp")
contrasts = unlist(strsplit(c("group:H_cj1:H_cj2"), ":", perl = T))
numerator_subset = unlist(strsplit(contrasts[2], ",", perl = T))
denominator_subset = unlist(strsplit(contrasts[3], ",", perl = T))
cellmeta = Seurat::FetchData(brain, vars = contrasts[1]) %>% tibble::rownames_to_column(var = "barcode")
numerator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
    numerator_subset) %>% dplyr::pull(barcode)
denominator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
    denominator_subset) %>% dplyr::pull(barcode)
Diff_exp = FindMarkers(brain, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, 
                               only.pos = TRUE,
                               ident.1 = numerator_subset,
                               ident.2 = denominator_subset, 
                               group.by = contrasts[1],
                               test.use = "wilcox")
Diff_exp = Diff_exp[abs(Diff_exp$avg_log2FC) > 0.25,]
Diff_exp = Diff_exp %>% tibble::rownames_to_column(var = "gene") %>%
        dplyr::rename(pvalue = p_val, padj = p_val_adj)
numerator_means = Matrix::rowMeans(SeuratObject::GetAssayData(brain,
        slot = "data")[Diff_exp$gene, numerator])
denominator_means = Matrix::rowMeans(SeuratObject::GetAssayData(brain,
        slot = "data")[Diff_exp$gene, denominator])
Diff_exp1 = Diff_exp %>% dplyr::mutate(FoldChange = 2^avg_log2FC, baseMean = 1/2 *
        (log2(numerator_means) + log2(denominator_means))) %>%
        dplyr::rename(log2FoldChange = avg_log2FC) %>% dplyr::select(gene,
        everything()) %>% dplyr::select(-baseMean)
colnames(Diff_exp1) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange")
res_Significant = dplyr::filter(Diff_exp1, `p-value` < 0.05,
            abs(log2FoldChange) > log2(1.5))
res_Significant[which(res_Significant$log2FoldChange > 0),
        "Regulation"] <- "Up"
res_Significant[which(res_Significant$log2FoldChange < 0),
        "Regulation"] <- "Down"
colnames(res_Significant) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange","Regulation")
write.table(Diff_exp, 
            "4.Diffexp/all_diff.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(res_Significant, 
            "4.Diffexp/diff_p<0.05_FC>1.5.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

#5 富集分析
#5.1 GO富集分析
dir.create("5.enrichment")
genes_symbol <- as.character(res_Significant$gene)
eg = bitr(genes_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
ego <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
GO_dot = dotplot(ego,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free") 
GO_bar = barplot(ego,split = "ONTOLOGY")+ facet_grid(ONTOLOGY~.,scales = "free")
res_plot <- CombinePlots(list(GO_dot,GO_bar), nrow=1)
ggsave("5.enrichment/GO_results_all.pdf", plot=res_plot, width = 12,height = 10)
ggsave("5.enrichment/GO_results_all.png", plot=res_plot, width = 12,height = 10)

#5.2 KEGG富集分析
# https://davidbioinformatics.nih.gov/conversion.jsp
# http://bioinfo.org/kobas/genelist/

#6 RCTD
dir.create("6.RCTD")
cortex_sc <- readRDS("MNN.rds")
counts_ref = GetAssayData(cortex_sc, 
                          assay = "RNA",
                          slot = "counts")
cell_types_ref <- cortex_sc@meta.data %>%  dplyr::select("subclass") %>%
  dplyr::mutate(cell_type =str_replace_all(!!!rlang::syms("subclass"), "-", "_") %>%   
  as.factor()) %>% .$cell_type
cell_types_ref <- setNames(cell_types_ref,
                           cortex_sc@meta.data  %>% rownames)
nUMI_ref <- cortex_sc@meta.data %>% dplyr::select("nCount_RNA") %>% .$nCount_RNA
nUMI_ref <- setNames(nUMI_ref,
                     cortex_sc@meta.data %>% rownames)
reference = spacexr::Reference(counts_ref,
                               cell_types_ref,
                               nUMI_ref)
seurat_to_SpatialRNA = function(st_ob, st_assay) {
  images = Seurat::Images(st_ob)  
  coords = do.call(rbind, lapply(images, function(slice) {
    spatial_coord = data.frame(st_ob[[slice]]@coordinates)
    colnames(spatial_coord) = c("in_tissue", "y", "x", "pxl_col_in_fullres", "pxl_row_in_fullres")
    spatial_coord$y = -1 * (spatial_coord$y) + 78
    spatial_coord = spatial_coord %>% dplyr::select(c(x, y))
    return(spatial_coord)
  }))
  counts = Seurat::GetAssayData(st_ob, assay = st_assay, slot = "counts")
  puck = spacexr::SpatialRNA(coords, counts) 
  return(puck)
}
spatialRNA = seurat_to_SpatialRNA(brain, "Spatial")
RCTD_ob = spacexr::create.RCTD(
  spatialRNA=spatialRNA,      
  reference=reference,        
  test_mode = FALSE,          
  gene_cutoff = 1.25e-4,     
  fc_cutoff = 0.5,           
  gene_cutoff_reg = 2e-04,  
  fc_cutoff_reg = 0.75, 
  UMI_max = 2e+05,
  counts_MIN = 10,    
  UMI_min_sigma = 0,       
  class_df = NULL,           
  CELL_MIN_INSTANCE = 1 ,     
  cell_type_names = NULL,     
  cell_type_profiles = NULL, 
  MAX_MULTI_TYPES = 4,        
  keep_reference = F,        
  CONFIDENCE_THRESHOLD = 10, 
  DOUBLET_THRESHOLD = 25,    
  max_cores = 20 )
myRCTD = spacexr::run.RCTD(RCTD_ob,
                           doublet_mode = "full")
norm_weights = normalize_weights(myRCTD@results$weights)%>% as.data.frame
brain@misc$RCTD_result = norm_weights
maxn = function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
metadata = norm_weights %>% 
  dplyr::mutate(
    top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
    top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
  )
brain = Seurat::AddMetaData(brain, metadata = metadata)
mat=normalize_weights(myRCTD@results$weights)%>%as.data.frame()
p = SPOTlight::spatial_scatterpie(se_obj = brain,
                                  cell_types_all = colnames(mat),
                                  img_path = "tissue_lowres_image.png",
                                  pie_scale = 0.4)
ggsave("6.RCTD/Celltype_RCTD.pdf", 
       plot = p,
       device = "pdf", 
       width = 14, 
       height = 12,
       dpi = 300)
p = SpatialDimPlot(brain,group.by = "top1_celltype")+
  theme( legend.key = element_rect(fill = NA, color = NA),
         legend.text = element_text(size = 10),
         plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("6.RCTD/top1_celltype_RCTD.pdf",
       plot = p, 
       device = "pdf", 
       width = 14,
       height = 12,
       dpi = 300)


