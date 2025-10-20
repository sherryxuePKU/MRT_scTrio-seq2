## load packages
library(Seurat)
library(ggplot2)
library(data.table)
library(ggtern)
library(ComplexHeatmap)
library(circlize)

## check packages' version
packageVersion("Seurat")

## state paths
count_dir <- "F:/Project/3P/data/RNA/umi_count/raw_data_count/"
save_count_dir <- "F:/Project/3P/data/RNA/umi_count/combined_count/MRT_UMI_count.txt"
infodir <- "F:/Project/3P/data/StatInfo/SampleInfo.txt"
save_seurat_dir <- paste0("F:/Project/3P/R_workspace/ST_", Sys.Date(), "_seurat.rda")
plot_outdir <- "F:/Project/3P/plot/"

## other parameter
min_gene_num <- 4000

## load colorlist
load("F:/Project/3P/R_workspace/colors.rda")

## marker list
TE_marker <- c("GATA3", "DAB2", "TFAP2C", "GATA2","CDX2","PTGES",
               "EMP2", "TGFBR3", "PDGFA", "KRT18","CLDN10", "PLAC8", "TRIML1")
PE_marker <- c("PDGFRA", "GATA4", "FOXA2","HNF1B", "COL4A1",  "FN1",
               "LINC00261","AMOTL1", "DPP4","SOX17")
EPI_marker <- c("SOX2",  "NANOG", "TDGF1","GDF3", "PRDM14",
                "NODAL","ARGFX", "DPPA2", "POU5F1")

## define function
merge_count_func <- function(indir, outdir, save_res=T){
  umilist <- list()
  fileName <- dir(indir)
  for (i in 1:length(fileName)) {
    umilist[[i]] <- fread(paste0(indir, fileName[i]), 
                          header = T, 
                          row.names = 1, 
                          stringsAsFactors = F)
  }
  umi_count <- do.call("cbind", umilist)
  if(save_res){
    write.table(umi_count, 
                file = outdir, 
                quote = F, sep = "\t", 
                col.names = T, row.names = T)
  }
  return(umi_count)
}

seurat_pipeline <- function(umi_count_dir=count_dir, 
                            info_dir=infodir, project="MRT", 
                            min_cell=3, min_gene=min_gene_num,
                            save.obj=T, save_path=save_seurat_dir){
  set.seed(0)
  data <- read.table(umi_count_dir,header = T, stringsAsFactors = F)
  info <- fread(info_dir, header = T, stringsAsFactors = F)
  info <- as.data.frame(info)
  rownames(info) <- info$Cell_id
  cell_vector <- intersect(info$Cell_id, colnames(data))
  seurat.obj <- CreateSeuratObject(counts = data[,cell_vector], min.cells = min_cell, min.features = min_gene, project = project)
  seurat.obj[["percent.mt"]]  <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  seurat.obj[["ID"]] <- info[colnames(seurat.obj), "ID"]
  seurat.obj[["Cell_id"]] <- rownames(seurat.obj@meta.data)
  seurat.obj[["Embryo"]] <- info[colnames(seurat.obj),"Embryo"]
  seurat.obj[["Grade"]] <- info[colnames(seurat.obj), "Grade"]
  seurat.obj[["Day"]] <- info[colnames(seurat.obj), "Day"]
  seurat.obj[["Batch"]] <- info[colnames(seurat.obj),"Batch"]
  seurat.obj[["Group"]] <- info[colnames(seurat.obj),"Group"]
  seurat.obj[["Index"]] <- info[colnames(seurat.obj), "Index"]
  seurat.obj[["Lineage"]] <- info[colnames(seurat.obj), "Lineage"]
  if(save.obj){
    if(is.null(save_path)){
      print("Please specify path to save!")
    } else {
      save(seurat.obj, file = save_path)
    }
  }
  return(seurat.obj)
}

myfeature_plot <- function(seurat.obj,markers,
                           min_expr=2,max_expr=4, 
                           col =c("lightgrey", "red"),
                           plot_dir="E:/OneDrive/3P/data/plot/"){
  mytheme <- theme_minimal()+theme(axis.ticks = element_blank(),
                                   axis.text = element_blank(),
                                   axis.title = element_blank(),
                                   panel.border = element_blank(),
                                   legend.key.size=unit(0.75,'cm'),
                                   legend.text = element_text(size = 15),
                                   legend.title = element_text(size = 20),
                                   legend.position = "none", 
                                   panel.grid = element_blank(),
                                   plot.title = element_text(hjust = 0.5, size = 25, face = "italic"))
  reduction_method <- "tsne"
  dim1 <- "tSNE_1"
  dim2 <- "tSNE_2"
  dim_used <- 1:2
  
  for (mk in markers) {
    coor_ratio <- (range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim1])[2] -
                     range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim1])[1])/
      (range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim2])[2] - 
         range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim2])[1])

    p <- FeaturePlot(seurat.obj, 
                     features = mk,cols = col, 
                     reduction = reduction_method,
                     order = T, dims = dim_used,
                     min.cutoff = min_expr,
                     max.cutoff = max_expr,
                     pt.size = 1.5)+
      coord_fixed(ratio = coor_ratio)+
      mytheme
    print(p)
  }
}

plot_lineage_ternary_func <- function(seurat.obj=ST.seurat,
                                      marker_dir,lineage=c("EPI", "PE", "TE"),
                                      group_by="Lineage",
                                      cols=MRT_colorlist[["Lineage_col"]]){
  Lineage_Marker <- read.table(marker_dir, header = T, stringsAsFactors = F)
  Lineage_Marker <- Lineage_Marker[,c("gene", "cluster")]
  colnames(Lineage_Marker) <- c("Gene_id", "Lineage")
  ggtern_df <- as.data.frame(t(FetchData(object = seurat.obj, slot = "data",
                                         vars = rownames(seurat.obj))))
  tmp <- list()
  for (i in c("EPI", "PE", "TE")) {
    tmp[[i]] <- colMeans(ggtern_df[Lineage_Marker[Lineage_Marker$Lineage==i,"Gene_id"],], 
                         na.rm = T)
  }
  Ternary_df <- as.data.frame(t(do.call("rbind", tmp)))
  Ternary_df$Var <- seurat.obj@meta.data[rownames(Ternary_df),group_by]
  
  ggtern(data=Ternary_df,aes(EPI,TE,PE, color=Var))+
    geom_point(size=1)+
    scale_color_manual(values = MRT_colorlist$Lineage_col)+
    theme_custom(col.T = cols["TE"], col.L = cols["EPI"],
                 col.R = cols["PE"],col.grid.minor = "grey90",
                 tern.panel.background = "white")+
    theme(legend.position = "right",
          legend.key = element_blank(),
          legend.title = element_blank())
  
}

## merge multiple umi_count.xls into one umi_count.xls
merge_count_func(indir = count_dir, outdir=save_count_dir, save_res = T)

## create seurat object and add metadata
ST.seurat <- seurat_pipeline(umi_count_dir = save_count_dir,
                             save_path = save_seurat_dir,
                             info_dir = infodir)
ST.seurat <-NormalizeData(ST.seurat, normalization.method = "LogNormalize", scale.factor = 1e+5, verbose = F)
ST.seurat <-FindVariableFeatures(ST.seurat, selection.method = "vst")
ST.seurat <- ScaleData(ST.seurat, features = rownames(ST.seurat))
ST.seurat <- RunPCA(ST.seurat,features = c(PE_marker, EPI_marker,TE_marker))

## fig.1.b
pdf(paste0(outdir, "MRT_", min_gene_num, ".reduce_dimension.pdf"))

coor_ratio <- calc_coord_ratio(seurat.obj = ST.seurat, reduction_method = reduction)
DimPlot(ST.seurat, reduction = "tsne", 
              pt.size = 2,
              cols = MRT_colorlist$Group_col,
              group.by = "Group") + 
  theme_bw(base_size = 25)+ 
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")+
  coord_equal(ratio = coor_ratio)

p2<- DimPlot(ST.seurat, reduction = "tsne", 
             pt.size = 2,
             cols = MRT_colorlist$Lineage_col,
             group.by = "Lineage") + 
  theme_bw(base_size = 25) + 
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")+
  coord_equal(ratio = coor_ratio)
print(p2)

dev.off()

## fig.1.c
Idents(ST.seurat) <- ST.seurat$Lineage
levels(ST.seurat)
new_lineage_marker <- FindAllMarkers(ST.seurat, test.use = "wilcox")
new_lineage_marker <- new_lineage_marker[new_lineage_marker$avg_log2FC >0,]
MRT_lineage_marker <- new_lineage_marker[!duplicated(new_lineage_marker$gene),]
write.table(MRT_lineage_marker,file = paste0(plot_outdir, "MRT_", min_gene_num, ".lineage_marker.txt"),
            quote = F, sep = "\t",row.names = F, col.names = T)

# MRT_lineage_marker <- read.table(file = paste0(outdir, "MRT_", min_gene_num, ".lineage_marker.txt"),
#                                  header = T, stringsAsFactors = F)

test <- MRT_lineage_marker %>% 
  filter(p_val_adj < 0.05 & pct.1 >.5 ) %>%
  select(avg_log2FC, cluster, gene) %>% 
  group_by(cluster) %>% 
  slice_max(n = 100, order_by=avg_log2FC) %>% 
  as.data.frame()

expr_scale <- ST.seurat[["RNA"]]@scale.data[test$gene,]
# 
# limit <- 1.5
# expr_scale[expr_scale> limit] <- limit
# expr_scale[expr_scale< -limit] <- (-limit)

anno_col <- data.frame(Lineage=ST.seurat$Lineage)
rownames(anno_col) <- rownames(ST.seurat@meta.data)

anno_colors <- list(Lineage=MRT_colorlist$Lineage_col)

col_anno <- HeatmapAnnotation(Lineage=ST.seurat$Lineage,Group=ST.seurat$Group,
                              col = list(Lineage=MRT_colorlist$Lineage_col,
                                         Group=MRT_colorlist$Group_col))
row_anno <- rowAnnotation(foo = anno_mark(at = which(test$gene %in% 
                                                       c(TE_marker, EPI_marker, PE_marker)), 
                                          labels = test$gene[test$gene %in% 
                                                               c(TE_marker, EPI_marker, PE_marker)],
                                          labels_gp = gpar(fontface="italic")))

colors_func <- colorRamp2(c(min(expr_scale, na.rm = T), mean(expr_scale, na.rm = T), max(expr_scale, na.rm = T)),
                          c("#91BFDB", "white", "#D73027"))

ht <- ComplexHeatmap::Heatmap(expr_scale,
                              show_row_names = F, show_column_names = F,
                              row_names_gp = gpar(fontface="italic"),
                              column_split = factor(anno_col$Lineage, 
                                                    levels = c("TE", "EPI", "PE")),
                              column_gap = unit(0, "mm"),
                              row_split = test$cluster, row_gap = unit(0, "mm"),
                              border = TRUE, border_gp = gpar(lty = 1, lwd = 1),
                              cluster_columns = T, cluster_column_slices = T,
                              column_dend_height = unit(8, "cm"), 
                              cluster_row = T,cluster_row_slices = T,
                              na_col = "lightgrey",
                              top_annotation = col_anno,
                              right_annotation = row_anno,
                              show_row_dend = T, show_column_dend = T,
                              col = colors_func,
                              name = "ht",
                              use_raster = T,
                              raster_device="CairoPNG",
                              raster_quality = 5)
pdf(paste0(outdir, "MRT_", min_gene_num, ".lineage_marker_enlarged.pdf"), width = 20,height = 15)
draw(ht)
dev.off()

## sup.fig.1.b 
VlnPlot(ST.seurat, features = c("nFeature_RNA"), ## nCount_RNA, percent.mt
        group.by = "Group",cols = MRT_colorlist$Group_col,
        adjust = 2,pt.size = 0)+
  theme_classic(base_size = 10)+
  labs(x="Group", y="nGene")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_blank())
ggsave(paste0(plot_outdir, "MRT_", min_gene_num, ".RNA_QC_nGene.pdf"),width = 3, height = 3)

## sup.fig.1.c
plot_lineage_ternary_func(seurat.obj = ST_seurat,
                          marker_dir = paste0(outdir, "MRT_", min_gene_num, ".lineage_marker.txt"))

## sup.fig.1.d
pdf(file=paste0(plot_outdir, "MRT_", min_gene_num, ".marker_featureplot.pdf"), width = 5, height = 5)
## EPI marker
myfeature_plot(ST.seurat,markers = "NANOG",min_expr = 1.5,max_expr = 3)
myfeature_plot(ST.seurat,markers = "SOX2",min_expr = 1,max_expr = 2)
myfeature_plot(ST.seurat,markers = "POU5F1",min_expr = 4,max_expr = 5)

## TE marker
myfeature_plot(ST.seurat,markers = "GATA3",min_expr = 2,max_expr = 3)
myfeature_plot(ST.seurat,markers = "GATA2",min_expr = 2,max_expr = 3)
myfeature_plot(ST.seurat,markers = "CDX2",min_expr = 2,max_expr = 3)

## PE marker
myfeature_plot(ST.seurat,markers = "GATA4",min_expr = 1,max_expr = 3)
myfeature_plot(ST.seurat,markers = "SOX17",min_expr = 2,max_expr = 3)
myfeature_plot(ST.seurat,markers = "FOXA2",min_expr = 1,max_expr = 3)
dev.off()

## sup.fig.1.e
DotPlot(ST.seurat,scale=T,
        cols = c("yellow", "red"),
        features = c(TE_marker,PE_marker, EPI_marker), 
        group.by = "Lineage")+
  theme_bw()+
  theme(axis.text.x = element_text(face = "italic", vjust = 1, 
                                   hjust = 1, angle = 45),
        axis.title = element_blank(),
        axis.text = element_text(size = 12.5))

ggsave(file=paste0(plot_outdir, "MRT_", min_gene_num, ".marker_dotplot.pdf"), width = 12, height=4)

