dir.create("10.SCENIC")
source("script/aux_rss.R")
source("script/dotHeatmap.R")
setwd("10.SCENIC")
dbs = c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(dbs) = c("500bp", "10kb") 
org="hgnc" 
dbDir="data/DATABASE"
scenicOptions <- initializeScenic(org = org ,dbDir = dbDir ,dbs = dbs ,nCores = 10)
scRNA_mnn = readRDS("data/Harmony.rds")
minCell4gene = round(0.01 * ncol(scRNA_mnn))
exprMat = scRNA_mnn@assays$RNA@counts 
genesKept <- geneFiltering(as.matrix(exprMat), 
		scenicOptions=scenicOptions, 
		minCountsPerGene=1,
		minSamples=minCell4gene)
exprMat_filtered <- exprMat[genesKept, ]
tf_names = getDbTfs( scenicOptions )
tf_names = CaseMatch(search=tf_names, match = rownames(scRNA_mnn))
arb.algo = import("arboreto.algo")
adjacencies = arb.algo$grnboost2(as.data.frame(t(as.matrix(exprMat_filtered))), tf_names=tf_names, verbose=T, seed=123L)
colnames(adjacencies) = c( "TF", "Target", "weight" )
saveRDS( adjacencies, file = getIntName(scenicOptions, "genie3ll") )
corrMat = bigcor(t(as.matrix(exprMat_filtered)),size = 2000, cores = 20, method = "spearman")
saveRDS(corrMat, file= getIntName(scenicOptions, "corrMat"))
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod = "top10perTarget")
scenicOptions <- initializeScenic(org = org ,dbDir = dbDir ,dbs = dbs ,nCores = 1)
runSCENIC_3_scoreCells(scenicOptions,log2(as.matrix(exprMat_filtered)+1))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC_mat = AUCell::getAUC(regulonAUC)
rownames(regulonAUC_mat) = gsub("_", "-", rownames(regulonAUC_mat))
regulonAUC_mat_out = regulonAUC_mat[-grep(pattern="-extended",rownames(regulonAUC_mat)),]
write.table(as.data.frame(regulonAUC_mat_out) %>% tibble::rownames_to_column(var = "regulon"),
            file.path("regulon_activity.xls"),
            sep = "\t", col.names =T, row.names =F)
scRNA_mnn[["SCENIC"]] = CreateAssayObject(counts = regulonAUC_mat)
scRNA_mnn = ScaleData(scRNA_mnn, assay = "SCENIC")
scRNA_mnn@tools$RunAUCell = regulonAUC
regulonTargetsInfo = loadInt(scenicOptions, "regulonTargetsInfo")
write.table(regulonTargetsInfo,
            file.path(output_dir, "0.1.TF_target_enrichment_annotation.xls"),
            sep = "\t", col.names =T, row.names =F, quote =F)
regulons <- loadInt(scenicOptions, "regulons")
sub_regulons = gsub(" .*","",rownames(regulonAUC_mat_out))
regulons = regulons[sub_regulons]
write.gmt(regulons, gmt_file = file.path(output_dir, "0.2.regulon_annotation.xls"))
Idents(object = scRNA_mnn) = scRNA_mnn@meta.data$clusters
cellInfo = scRNA_mnn@meta.data
col_anno = as.data.frame(scRNA_mnn@meta.data) %>% rownames_to_column(var="barcodes")
col_anno = col_anno[,c("barcodes","clusters")]
col_anno = col_anno %>% arrange(clusters) %>% column_to_rownames(var="barcodes")
regulonAUC_plotdata = regulonAUC_mat_out[,rownames(col_anno)]
bks <- unique(c(seq(-2.5,0, length=100),  seq(0,2.5, length=100)))
color_use = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde")
cluster_levels <- as.character(unique(col_anno$clusters))
names(color_use) <- cluster_levels
annotation_colors <- list(
  clusters = color_use
)
pdf("10.SCENIC/1.1.regulon_activity_heatmap_groupby_cells.pdf",height = 8) 
pheatmap::pheatmap( regulonAUC_plotdata,
                    scale = "row",
                    cluster_cols=F,
                    cluster_rows=F,
                    show_colnames= F,
                    color=colorRampPalette(c("#406AA8","white","#D91216"))(200),
                    annotation_col = col_anno,
                    annotation_colors = annotation_colors,
                    treeheight_col=10,
                    border_color=NA,breaks=bks,fontsize_row=6)
dev.off() 
regulonActivity_byclusters <- sapply(split(rownames(cellInfo), cellInfo$clusters), function(cells) 
				rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byclusters_Scaled <- t(scale(t(regulonActivity_byclusters), center = T, scale=T)) 
df = as.data.frame(regulonActivity_byclusters_Scaled) %>% tibble::rownames_to_column(var = "regulon")
write.table(df,file.path(output_dir, "1.2.centered_regulon_activity_groupby_design.xls"),
            sep = "\t", col.names =T, row.names =F, quote =F)
regulonAUC_plotdata <- regulonActivity_byclusters_Scaled[1:10, ]
pdf("10.SCENIC/1.3.regulon_activity_heatmap.pdf")
pheatmap::pheatmap( regulonAUC_plotdata,
                    cellwidth = 18,
                    cellheight = 18,
                    color=colorRampPalette(c("#406AA8","white","#D91216"))(299),
                    angle_col = 45,
                    treeheight_col=20, 
                    treeheight_row=20,
                    border_color=NA)
dev.off()
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "clusters"])
pdf('10.SCENIC/2.2.RSS_ranking_plot.pdf') 
setName = "1"
n = 5
rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  setName = setName
 p4 =  ggplot(thisRss, aes(x=rank, y=rss)) + 
    geom_point(color = "grey50", size = 1) + 
    ggtitle(setName) + 
    geom_point(data = subset(thisRss,rank < n+1),
               color = "red",size = 2)+
    geom_label_repel(aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     na.rm=TRUE) +
    theme_classic()
p4
dev.off()
write.table(thisRss[1:5,], file.path(output_dir, "2.1.regulon_RSS_annotation.xls"),
            sep = "\t", col.names = T, row.names =F, quote =F)
regulons_csi <- calculate_csi(regulonAUC,calc_extended = FALSE)
csi_csi_wide <- regulons_csi %>% tidyr::spread(regulon_2,CSI)
future_rownames <- csi_csi_wide$regulon_1
csi_csi_wide <- as.matrix(csi_csi_wide[,2:ncol(csi_csi_wide)])
rownames(csi_csi_wide) <- future_rownames
regulons_hclust <- hclust(dist(csi_csi_wide,method = "euclidean"))
nclust = 4
clusters <- cutree(regulons_hclust,k= nclust)
clusters_df <- data.frame("regulon" = names(clusters),"csi_module" = clusters)
cellinfo = scRNA_mnn@meta.data
csi_cluster_activity_wide <- calc_csi_module_activity(clusters_df,regulonAUC, cellinfo,cell_type_column = "seurat_clusters")
rownames(csi_cluster_activity_wide) = paste0("module",c(1:nclust))
plot = pheatmap::pheatmap(csi_cluster_activity_wide,
                    show_colnames = TRUE,
                    show_rownames = TRUE,
                    scale = "row",
                    color = viridis::viridis(n = 10),
                    cellwidth = 24,
                    cellheight = 24,
                    cluster_cols = TRUE,
                    cluster_rows = TRUE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean")
ggsave("10.SCENIC/3.3.csi_module_activity_heatmap.pdf",
    plot = plot,
    width = 8,
    height = 8,
    dpi = 1000,
    limitsize = F)
regulons_csi <- calculate_csi(regulonAUC,calc_extended = FALSE)
plot = plot_csi_modules(regulons_csi,nclust = nclust)
pdf("10.SCENIC/3.2.regulons_csi_correlation_heatmap.pdf") 
print(plot)
dev.off() 
write.table(clusters_df, "10.SCENIC/3.1.csi_module_annotation.xls", sep = "\t", col.names = T, row.names =F, quote =F)