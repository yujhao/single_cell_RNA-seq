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


#4.1 GO富集分析
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
ggsave("4.Diffexp/GO_results_all.pdf", plot=res_plot, width = 12,height = 10)
ggsave("4.Diffexp/GO_results_all.png", plot=res_plot, width = 12,height = 10)

#4.2 KEGG富集分析
# https://davidbioinformatics.nih.gov/conversion.jsp
# http://bioinfo.org/kobas/genelist/
