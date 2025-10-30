##-------------------------load packages-------------------------
library(Seurat)
library(ggplot2)
library(data.table)
library(ggtern)
library(ComplexHeatmap)
library(circlize)

## check packages' version
packageVersion("Seurat")

##-------------------------state paths-------------------------
count_dir <- "/path/to/your/umi_count/raw_data_count/"
save_count_dir <- "/path/to/your/umi_count/combined_count/MRT_UMI_count.txt"
infodir <- "/path/to/your/MRT_RNA.CellsPass_info.txt"
save_seurat_dir <- paste0("/path/to/your/ST_", Sys.Date(), "_seurat.rda")
plot_outdir <- "/path/to/your/plot/"

##-------------------------hyperparameters-------------------------
MIN_GENE_NUM <- 4000

##-------------------------load colorlist-------------------------
load("/path/to/your/metadata/colors.rda")

##-------------------------marker list-------------------------
TE_marker <- c("GATA3", "DAB2", "TFAP2C", "GATA2","CDX2","PTGES",
               "EMP2", "TGFBR3", "PDGFA", "KRT18","CLDN10", "PLAC8", "TRIML1")
PE_marker <- c("PDGFRA", "GATA4", "FOXA2","HNF1B", "COL4A1",  "FN1",
               "LINC00261","AMOTL1", "DPP4","SOX17")
EPI_marker <- c("SOX2",  "NANOG", "TDGF1","GDF3", "PRDM14",
                "NODAL","ARGFX", "DPPA2", "POU5F1")

##-------------------------define function-------------------------
source("analysis/RNA_transcriptome/helper_function.R")

##-------------------------create seurat object-------------------------
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

##-------------------------fig.1.b-------------------------
pdf(paste0(outdir, "MRT_", MIN_GENE_NUM, ".reduce_dimension.pdf"))

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

##-------------------------fig.1.c-------------------------
Idents(ST.seurat) <- ST.seurat$Lineage
levels(ST.seurat)
new_lineage_marker <- FindAllMarkers(ST.seurat, test.use = "wilcox")
new_lineage_marker <- new_lineage_marker[new_lineage_marker$avg_log2FC >0,]
MRT_lineage_marker <- new_lineage_marker[!duplicated(new_lineage_marker$gene),]
write.table(MRT_lineage_marker,file = paste0(plot_outdir, "MRT_", MIN_GENE_NUM, ".lineage_marker.txt"),
            quote = F, sep = "\t",row.names = F, col.names = T)

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

ht <- Heatmap(expr_scale,
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
pdf(paste0(outdir, "MRT_", MIN_GENE_NUM, ".lineage_marker_enlarged.pdf"), width = 20,height = 15)
draw(ht)
dev.off()

##-------------------------sup.fig.1.b-------------------------
VlnPlot(ST.seurat, features = c("nFeature_RNA"), ## nCount_RNA, percent.mt
        group.by = "Group",cols = MRT_colorlist$Group_col,
        adjust = 2,pt.size = 0)+
  theme_classic(base_size = 10)+
  labs(x="Group", y="nGene")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_blank())
ggsave(paste0(plot_outdir, "MRT_", MIN_GENE_NUM, ".RNA_QC_nGene.pdf"),width = 3, height = 3)

##-------------------------sup.fig.1.c-------------------------
plot_lineage_ternary_func(seurat.obj = ST_seurat,
                          marker_dir = paste0(outdir, "MRT_", MIN_GENE_NUM, ".lineage_marker.txt"))

##-------------------------sup.fig.1.d-------------------------
pdf(file=paste0(plot_outdir, "MRT_", MIN_GENE_NUM, ".marker_featureplot.pdf"), width = 5, height = 5)
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

##-------------------------sup.fig.1.e-------------------------
DotPlot(ST.seurat,scale=T,
        cols = c("yellow", "red"),
        features = c(TE_marker,PE_marker, EPI_marker), 
        group.by = "Lineage")+
  theme_bw()+
  theme(axis.text.x = element_text(face = "italic", vjust = 1, 
                                   hjust = 1, angle = 45),
        axis.title = element_blank(),
        axis.text = element_text(size = 12.5))

ggsave(file=paste0(plot_outdir, "MRT_", MIN_GENE_NUM, ".marker_dotplot.pdf"), width = 12, height=4)

