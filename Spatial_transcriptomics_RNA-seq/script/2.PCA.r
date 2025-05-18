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