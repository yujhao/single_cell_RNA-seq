#7 邻域分析
dir.create("7.neighbor")
data_ob <- readRDS("rctd.rds")
if(length(Seurat::Images(data_ob))>1){
	futile.logger::flog.info("报错：该分析一次只针对一个样本")
	quit()
}
center <- "MalignantEpithelialCells"
locat <- data_ob@images[[1]]@coordinates[, c("row", "col")]
distest <- dist(locat, p = 2)
center_barcode <- rownames(data_ob@meta.data)[grep(center, unlist(data_ob[["top1_celltype"]]))]
distest2 <- as.data.frame(as.matrix(distest)) %>% .[rownames(data_ob@meta.data), center_barcode]
distest2$dist <- apply(distest2, 1, min)
center_distance <- paste0(center, "_distance")
data_ob[[center_distance]] <- distest2$dist
plot_data <- data.frame(celltype = factor(data_ob@meta.data[, "top1_celltype"],
											levels = unique(data_ob@meta.data[, "top1_celltype"])),
						dist = distest2$dist)
    
p1_slice <- SpatialPlot(data_ob,
										features = paste0(center, "_distance"),
										min.cutoff = 0.1,
										combine = FALSE,
										ncol = 2,
										cols="spectral",
										alpha =1,
										pt.size.factor = 1,
										image.alpha = 0.01,
										images = Seurat::Images(data_ob)[1],
										crop = FALSE
)
p1_density <- ggplot2::ggplot(plot_data) +
				ggplot2::stat_density(ggplot2::aes(x = dist, colour = celltype),
									geom = "line",
									position = "identity",
									size = 0.2) +
				ggplot2::scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02")) +
				ggplot2::theme_set(ggplot2::theme_bw()) +
				ggplot2::theme(panel.grid = ggplot2::element_blank()) +
				ggplot2::xlab(paste0("distance to ", center, " of ", Seurat::Images(data_ob)[1])) +
				ggplot2::scale_y_continuous(limits = c(0, 6))

layout <- glue::glue("{strrep('A',1)}{strrep('#',1)}\n{strrep('B',2)}")
p1_slice_density <- patchwork::wrap_plots(A = p1_slice[[1]],
                                              B = p1_density,
                                              widths=c(1+stringr::str_length(paste0(center, "_distance"))/8,2),
                                              design = layout)

ggsave("7.neighbor/distance_and_density.pdf", 
       plot = p1_slice_density,
       device = "pdf", 
       width = 14, 
       height = 12,
       dpi = 300)

#num取外围层数，推荐一层，两层
get_loc <- function(loc, barcode, num) {
    row <- loc[barcode,]$row
    col <- loc[barcode,]$col
    col_list <- list()
    for (i in 1:(num + 1)) {
      part_col <- loc[loc$row %in% c((row - (i - 1)):(row + (i - 1))) & loc$col %in% c((col - 2 * num + (i - 1)):(col + 2 * num - (i - 1))),]
        col_list[[i]] <- rownames(part_col)
      }
      col_list <- unique(unlist(col_list))
}
barcodes <- rownames(data_ob@meta.data[data_ob[["top1_celltype"]][,1] == center,])
allbarcode <- lapply(barcodes, get_loc, loc = locat, num = 2)
allbarcode <- unique(unlist(allbarcode))
spottype <- data.frame(celltype = data_ob[["top1_celltype"]])
spottype$spottype <- NA
spottype[rownames(spottype) %in% allbarcode & spottype[, 1] == center, 2] <- "intra_spot" #内点
spottype[rownames(spottype) %in% allbarcode & spottype[, 1] != center, 2] <- "inter_spot" #间点
spottype[is.na(spottype[, 2]), 2] <- "distal_spot" #外点
center_spot <- paste0(center, "_spot_type")
data_ob[[center_spot]] <- factor(spottype[, 2], levels = c("intra_spot", "inter_spot", "distal_spot"))
p2_celltype <- SpatialPlot(data_ob,
                            group.by ="top1_celltype",
                            pt.size.factor = 1,
                            image.alpha = 0,
                            crop = FALSE)
p2_spottype <- SpatialPlot(data_ob,
                            group.by = center_spot,
                            ncol = 2,
                            alpha = c(1.1, 1.5),
                            pt.size.factor = 1,
                            image.alpha = 0,
                            images = Seurat::Images(data_ob)[1],
                            crop = FALSE)
p2_celltype_spottype <- patchwork::wrap_plots(A = p2_celltype[[1]],
                                              B = p2_spottype[[1]],
                                              widths=c(4+(data_ob@meta.data[["top1_celltype"]] %>% unique %>%
                                                        stringr::str_length()%>% max%>%+5) /10,
                                                        4+stringr::str_length(center_spot)/10),
                                              design = "AB")
ggsave("7.neighbor/spottype_and_celltype.pdf", 
       plot = p2_celltype_spottype,
       device = "pdf", 
       width = 14, 
       height = 12,
       dpi = 300)

#当细胞类型是直接注释而不是 RCTD或 SpotLight 鉴定出时,metadata 中没有细胞类型的百分比信息,无法做 不同细胞类型在不同分层中的占比箱线图和不同层级中不同细胞类型的占比箱线图
if (center %in% colnames(data_ob@meta.data)){
      futile.logger::flog.info("step4: 分面boxplot图-celltype")
      celltype <- c(names(table(data_ob[["top1_celltype"]])), paste0(center, "_spot_type")) %>% data_ob[[.]]
      print(head(celltype))
      my_comparisons <- list(c("intra_spot", "inter_spot"), c("intra_spot", "distal_spot"), c("inter_spot", "distal_spot"))
      box_list <- list()
      for (i in 1:(dim(celltype)[2] - 1)) {
        p <- ggpubr::ggboxplot(celltype[, c(i, dim(celltype)[2])],
                              y = names(celltype)[i],
                              x = center_spot,
                              fill  = center_spot,
                              palette =c("#FF5733", "#ED9481", "grey"),
                              add = "jitter", shape = NA) +
            ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format")+
            ggplot2::theme(text = element_text(size = 10))+
            ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
            ggplot2::guides(fill = FALSE)
        box_list[[i]] <- p
      }
      p3 <- patchwork::wrap_plots(box_list)
      ggsave("7.neighbor/spot_type_of_different_celltype.pdf",
                                plot = p3,
                                width = 15,
                                height = 15,
                                dpi = 300, bg="white")
      write.table(celltype,"7.neighbor/spot_type_of_different_celltype.xls",sep = "\t",quote = F,row.names = F)
      futile.logger::flog.info("step5: 分面boxplot图-region")
      spots_type <-unique(data_ob@meta.data[,center_spot]) #names(table(data_ob[[center_spot]]))
    
      spot_celltype <- list()
      c_name <- unique(data_ob@meta.data[,"top1_celltype"])[!unique(data_ob@meta.data[,"top1_celltype"]) %in% c("MalignantEpithelialCells")] %>%as.character
      print(c_name)
      for (i in 1:length(spots_type)) {
        data_names = as.character(spots_type[i])
        data <- data_ob@meta.data %>%
                tibble::rownames_to_column("barcodes")%>%
                dplyr::filter(!!rlang::sym(center_spot)== spots_type[i]) %>% 
                dplyr::select(c("barcodes",c_name))  %>%
                tibble::column_to_rownames("barcodes") %>%
                reshape2::melt(, variable.name = 'celltype', value.name = 'ratio')
        write.table(data,"7.neighbor/celltype_of_different_spottype.xls",sep = "\t",quote = F,row.names = F)
        #colors<-colors[c_name]
        p <- ggpubr::ggboxplot(data,
                              y = "ratio",
                              x = "celltype",
                              fill = "celltype",
                              add = "jitter",
                              shape = NA,
                              title = spots_type[i]) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 10, hjust = 1))+
            ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
            ggplot2::guides(fill = FALSE)
        spot_celltype[[i]] <- p
      }
      p4 <- patchwork::wrap_plots(spot_celltype)
      ggsave("7.neighbor/celltype_of_different_spottype.pdf",
                                plot = p4,
                                width = 9,
                                height = 5,
                                dpi = 300, bg="white")
    }else{
      futile.logger::flog.info("由于数据中没有细胞类型百分比信息，无法提供不同细胞类型在不同分层中的占比箱线图和不同层级中不同细胞类型的占比箱线图")
    }