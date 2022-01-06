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
count_dir <- "F:/Project/3P/data/RNA/umi_count/raw_data_count/"
save_count_dir <- "F:/Project/3P/data/RNA/umi_count/combined_count/MRT_UMI_count.txt"
infodir <- "F:/Project/3P/data/StatInfo/SampleInfo.txt"
# save_seurat_dir <- paste0("F:/Project/3P/R_workspace/ST_", Sys.Date(), "_seurat.rda")
plot_outdir <- "F:/Project/3P/plot/MRT"

seurat_dir <- paste0(plot_outdir, "_", min_gene_num, ".seurat.rda")
load(seurat_dir)

##-------------------------other parameter-------------------------
min_gene_num <- 4000

##-------------------------load colors-------------------------
load("F:/Project/3P/R_workspace/colors.rda")
MRT_colorlist$Group_col["ST"] <- "#f5bd2a"
MRT_colorlist$Group_col["ICSI"] <- "#704d9c"

##-------------------------marker list-------------------------
TE_marker <- c("GATA3", "DAB2", "TFAP2C", "GATA2","CDX2","PTGES",
               "EMP2", "TGFBR3", "PDGFA", "KRT18","CLDN10", "PLAC8", "TRIML1")
PE_marker <- c("PDGFRA", "GATA4", "FOXA2","HNF1B", "COL4A1",  "FN1",
               "LINC00261","AMOTL1", "DPP4","SOX17")
EPI_marker <- c("SOX2",  "NANOG", "TDGF1","GDF3", "PRDM14",
                "NODAL","ARGFX", "DPPA2", "POU5F1")

##-------------------------define function-------------------------
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

seurat_pipeline <- function(
  umi_count_dir=count_dir, 
  info_dir=infodir, project="MRT", 
  min_cell=3, min_gene=min_gene_num,
  save.obj=T, save_path=save_seurat_dir
  ){
  set.seed(0)
  data <- read.table(umi_count_dir,header = T, stringsAsFactors = F)
  info <- fread(info_dir, header = T, stringsAsFactors = F)
  info <- as.data.frame(info)
  rownames(info) <- info$Cell_id
  cell_vector <- intersect(info$Cell_id, colnames(data))
  ST.seurat <- CreateSeuratObject(counts = data[,cell_vector], min.cells = min_cell, min.features = min_gene, project = project)
  ST.seurat[["percent.mt"]]  <- PercentageFeatureSet(ST.seurat, pattern = "^MT-")
  ST.seurat[["ID"]] <- info[colnames(ST.seurat), "ID"]
  ST.seurat[["Cell_id"]] <- rownames(ST.seurat@meta.data)
  ST.seurat[["Embryo"]] <- info[colnames(ST.seurat),"Embryo"]
  ST.seurat[["Grade"]] <- info[colnames(ST.seurat), "Grade"]
  ST.seurat[["Day"]] <- info[colnames(ST.seurat), "Day"]
  ST.seurat[["Batch"]] <- info[colnames(ST.seurat),"Batch"]
  ST.seurat[["Group"]] <- info[colnames(ST.seurat),"Group"]
  ST.seurat[["Index"]] <- info[colnames(ST.seurat), "Index"]
  ST.seurat[["Lineage"]] <- info[colnames(ST.seurat), "Lineage"]
  if(save.obj){
    if(is.null(save_path)){
      print("Please specify path to save!")
    } else {
      save(ST.seurat, file = save_path)
    }
  }
  return(ST.seurat)
}

myfeature_plot <- function(
  ST.seurat,markers,
  min_expr=2,max_expr=4, 
  col =c("lightgrey", "red"),
  plot_dir="E:/OneDrive/3P/data/plot/"
  ){
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
    
    p <- FeaturePlot(ST.seurat, 
                     features = mk,cols = col, 
                     reduction = reduction_method,
                     order = T, dims = dim_used,
                     min.cutoff = min_expr,
                     max.cutoff = max_expr,
                     pt.size = 1.5)+
      coord_fixed(ratio = coor_ratio)+
      mytheme
    # print(p)
    return(p)
  }
}

plot_lineage_ternary_func <- function(
  ST.seurat=ST.seurat,
  marker_dir,lineage=c("EPI", "PE", "TE"),
  group_by="Lineage",
  cols=MRT_colorlist[["Lineage_col"]]
  ){
  Lineage_Marker <- read.table(marker_dir, header = T, stringsAsFactors = F)
  Lineage_Marker <- Lineage_Marker[,c("gene", "cluster")]
  colnames(Lineage_Marker) <- c("Gene_id", "Lineage")
  ggtern_df <- as.data.frame(t(FetchData(object = ST.seurat, slot = "data",
                                         vars = rownames(ST.seurat))))
  tmp <- list()
  for (i in c("EPI", "PE", "TE")) {
    tmp[[i]] <- colMeans(ggtern_df[Lineage_Marker[Lineage_Marker$Lineage==i,"Gene_id"],], 
                         na.rm = T)
  }
  Ternary_df <- as.data.frame(t(do.call("rbind", tmp)))
  Ternary_df$Var <- ST.seurat@meta.data[rownames(Ternary_df),group_by]
  
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

calc_coord_ratio <- function(reduction_method="pca",dim=c(1,2),ST.seurat){
  #reduction_method <- "pca"
  if(reduction_method=="umap"){
    dim1 <- "UMAP_1"
    dim2 <- "UMAP_2"
  } else if(reduction_method=="tsne"){
    dim1 <- "tSNE_1"
    dim2 <- "tSNE_2"
  }else if(reduction_method=="pca"){
    dim1 <- paste0("PC_", dim[1])
    dim2 <- paste0("PC_", dim[2])
  }
  coor_ratio <- (range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim1])[2] -
                   range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim1])[1])/
    (range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim2])[2] - 
       range(ST.seurat@reductions[[reduction_method]]@cell.embeddings[,dim2])[1])
  return(coor_ratio)
}

circos_hclust <- function(info, hclust.out,levels= c("Group", "Lineage"), col_list){
  
  cir_df <- info[hclust.out$tree_col$labels[hclust.out$tree_col$order],levels]
  N <- length(levels)
  for(i in 1:N){
    t <- names(table(cir_df[,levels[i]]))
    M <- length(t)
    col_name <- paste0(levels[i], "_col")
    cir_df[,col_name] <- "white"
    for (j in 1:M) {
      cir_df[cir_df[,levels[i]]==t[j], col_name] <- col_list[[col_name]][t[j]]
    }
  }
  
  nc <- dim(cir_df)[1]
  
  dend  <- as.dendrogram(hclust.out$tree_col)
  dend_height  <- attr(dend, "height")
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0, start.degree = 0)
  circos.initialize("a", xlim = c(0,nc))
  circos.track(ylim = c(0, N), bg.border = NA,
               panel.fun = function(x, y) {
                 for (i in 1:N) {
                   circos.rect(1:nc-1, rep(N-i, nc),1:nc, rep(N-i+1, nc),
                               border = cir_df[,N+i],col = cir_df[,N+i])
                 }
               }
  )
  circos.track(ylim = c(0, dend_height), bg.border = NA, 
               track.height = 0.4, panel.fun = function(x, y) {circos.dendrogram(dend)})
  circos.clear()
}

stat_CNV <-function(
  x=infer_result_total,
  chromosome=paste0("chr",1:22),
  embryo_list=c("E10","E11","E12","E13","E14","E15"),
  cell_info=ST_info
){
  ## output for CNV of each cell
  CNV_Summary_Matrix <- as.data.frame(matrix(rep(0,ncol(x)*22),ncol=22))
  colnames(CNV_Summary_Matrix) <- paste("chr",c(1:22),sep="")
  rownames(CNV_Summary_Matrix)<-colnames(x)
  ## output for CNV of each embryo
  CNV_Summary_Matrix_embryo <- as.data.frame(matrix(rep(0,length(embryo_list)*22),ncol=22))
  colnames(CNV_Summary_Matrix_embryo) <- paste0("chr", c(1:22), sep= "")
  rownames(CNV_Summary_Matrix_embryo) <- embryo_list
  ## calc
  for (j in embryo_list){
    select_col <- cell_info[cell_info$Embryo==j,"Cell_id"]
    embryo_x <- x[,select_col]
    for (i in 1:length(chromosome)){
      chr <-chromosome[i]
      selected_row <- gene_bed[gene_bed$chr==chr,"gene_symbol"]
      sub_data<-embryo_x[selected_row,]
      #sub_data<- sub_data[,-ncol(sub_data)]
      #sub_data <- sub_data[,-c(1:2)]
      ## count CNV for each cell
      gain_sc <- apply(sub_data,2,function(x){ratio<-length(x[x>1.05])/length(x);return(ratio)})
      lost_sc <- apply(sub_data,2,function(x){ratio<-length(x[x<0.95 & x>=0])/length(x);return(ratio)}) 
      tmp_cell <- rep(0,length(gain_sc))
      for(k in 1:length(gain_sc)){
        if(i<23){
          if(gain_sc[k]>=0.6)tmp_cell[k]=1
          else if(lost_sc[k]>=0.6)tmp_cell[k]=-1
        }
        # else{
        #   if(gain[k]>0.50)tmp[k]=1
        #   else if(lost[k]>0.50)tmp[k]=-1
        # }
      }
      CNV_Summary_Matrix[select_col,i]<-tmp_cell
      
      ## count CNV for each embryo
      if(length(tmp_cell[tmp_cell<0])>=3 & length(tmp_cell[tmp_cell>0])<3){
        CNV_Summary_Matrix_embryo[j,i] <- -1
      }
      if(length(tmp_cell[tmp_cell>0])>=3 & length(tmp_cell[tmp_cell<0])<3){
        CNV_Summary_Matrix_embryo[j,i] <- 1
      }
      if(length(tmp_cell[tmp_cell>0])>=3 & length(tmp_cell[tmp_cell<0])>=3){
        CNV_Summary_Matrix_embryo[j,i] <- 2
      }
      # else if (length(tmp_cell[tmp_cell==-1])>=3){
      #   CNV_Summary_Matrix_embryo[j,i] <- -1
      # }
    }
  }
  output <- list(sc=CNV_Summary_Matrix,
                 embryo=CNV_Summary_Matrix_embryo)
  return(output)
}

##-------------------------merge umi_count.xls-------------------------
merge_count_func(indir = count_dir, outdir=save_count_dir, save_res = T)

##-------------------------create seurat object-------------------------
ST.seurat <- seurat_pipeline(umi_count_dir = save_count_dir,
                             save_path = save_seurat_dir,
                             info_dir = infodir)
ST.seurat <-NormalizeData(ST.seurat, normalization.method = "LogNormalize", scale.factor = 1e+5, verbose = F)
ST.seurat <-FindVariableFeatures(ST.seurat, selection.method = "vst")
ST.seurat <- ScaleData(ST.seurat, features = rownames(ST.seurat))
ST.seurat <- RunPCA(ST.seurat,features = c(PE_marker, EPI_marker,TE_marker))
ST.seurat <- RunTSNE(ST.seurat,dims = 1:5)

##-------------------------fig.1.b-------------------------
pdf(paste0(outdir, "MRT_", min_gene_num, ".reduce_dimension.png"))

reduction <- "tsne"
coor_ratio <- calc_coord_ratio(ST.seurat = ST.seurat, reduction_method = reduction)
p1 <- DimPlot(ST.seurat, reduction = "tsne", 
        pt.size = 2,
        cols = MRT_colorlist$Group_col,
        group.by = "Group") + 
  theme_bw(base_size = 25)+ 
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")+
  coord_equal(ratio = coor_ratio)

p2 <- DimPlot(ST.seurat, reduction = "tsne", 
             pt.size = 2,
             cols = MRT_colorlist$Lineage_col,
             group.by = "Lineage") + 
  theme_bw(base_size = 25) + 
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")+
  coord_equal(ratio = coor_ratio)
print(p1+p2)

dev.off()

##-------------------------fig.1.c-------------------------
Idents(ST.seurat) <- ST.seurat$Lineage
levels(ST.seurat)
new_lineage_marker <- FindAllMarkers(ST.seurat, test.use = "wilcox")
new_lineage_marker <- new_lineage_marker[new_lineage_marker$avg_log2FC >0,]
MRT_lineage_marker <- new_lineage_marker[!duplicated(new_lineage_marker$gene),]
write.table(MRT_lineage_marker,file = paste0(plot_outdir, "MRT_", min_gene_num, ".lineage_marker.txt"),
            quote = F, sep = "\t",row.names = F, col.names = T)

MRT_lineage_marker <- read.table(file = paste0(plot_outdir, "MRT_", min_gene_num, ".lineage_marker.txt"),
                                 header = T, stringsAsFactors = F)
library(dplyr)
test <- MRT_lineage_marker %>% 
  filter(p_val_adj < 0.05 & pct.1 >.5 ) %>%
  select(avg_log2FC, cluster, gene) %>% 
  group_by(cluster) %>% 
  slice_max(n = 100, order_by=avg_log2FC) %>% 
  as.data.frame()

library(ComplexHeatmap)
library(circlize)
expr_scale <- ST.seurat[["RNA"]]@scale.data[test$gene,]

limit <- 1.5
expr_scale[expr_scale> limit] <- limit
expr_scale[expr_scale< -limit] <- (-limit)

anno_col <- data.frame(Lineage=ST.seurat$Lineage)
rownames(anno_col) <- rownames(ST.seurat@meta.data)

anno_colors <- list(Lineage=MRT_colorlist$Lineage_col)

col_anno <- HeatmapAnnotation(
  Lineage=ST.seurat$Lineage,
  Group=ST.seurat$Group,
  col = list(Lineage=MRT_colorlist$Lineage_col,
             Group=MRT_colorlist$Group_col)
)
row_anno <- rowAnnotation(
  foo = anno_mark(at = which(test$gene %in% c(TE_marker, EPI_marker, PE_marker)), 
                  labels = test$gene[test$gene %in% c(TE_marker, EPI_marker, PE_marker)],
                  labels_gp = gpar(fontface="italic"))
)

colors_func <- colorRamp2(
  c(min(expr_scale, na.rm = T), 
    mean(expr_scale, na.rm = T), 
    max(expr_scale, na.rm = T)),
  c("#91BFDB", "white", "#D73027")
)

ht <- ComplexHeatmap::Heatmap(
  expr_scale,
  show_row_names = F, 
  show_column_names = F,
  row_names_gp = gpar(fontface="italic"),
  column_split = factor(anno_col$Lineage, 
                        levels = c("TE", "EPI", "PE")),
  column_gap = unit(0, "mm"),
  row_split = test$cluster, row_gap = unit(0, "mm"),
  border = TRUE, border_gp = gpar(lty = 1, lwd = 1),
  cluster_columns = T, cluster_column_slices = T,
  column_dend_height = unit(8, "cm"), 
  cluster_row = T,
  cluster_row_slices = T,
  na_col = "lightgrey",
  top_annotation = col_anno,
  right_annotation = row_anno,
  show_row_dend = T, 
  show_column_dend = F,
  col = colors_func,
  name = "ht",
  use_raster = T,
  raster_device="CairoPNG",
  raster_quality = 5
)
png(
  paste0(plot_outdir, "MRT_", 
         min_gene_num, 
         ".lineage_marker_enlarged.png"),
  width = 640, height = 480, units = "px"
)
draw(ht)
dev.off()

##-------------------------sup.fig.1.b-------------------------
VlnPlot(ST.seurat, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ## nCount_RNA, percent.mt
        group.by = "Group",cols = MRT_colorlist$Group_col,
        adjust = 2,pt.size = 0)
# +theme_classic(base_size = 10)+
#   labs(x="Group", y="nGene")+
#   theme(axis.title.x = element_blank(),
#         legend.position = "none",
#         plot.title = element_blank())
ggsave(paste0(plot_outdir, "MRT_", min_gene_num, ".RNA_QC.png"))

##-------------------------sup.fig.1.c-------------------------
plot_lineage_ternary_func(ST.seurat = ST_seurat,
                          marker_dir = paste0(outdir, "MRT_", min_gene_num, ".lineage_marker.txt"))

##-------------------------sup.fig.1.d-------------------------
# pdf(file=paste0(plot_outdir, "MRT_", min_gene_num, ".marker_featureplot.pdf"), width = 5, height = 5)
## EPI marker
p1 <- myfeature_plot(ST.seurat,markers = "NANOG",min_expr = 1.5,max_expr = 3)
p2 <- myfeature_plot(ST.seurat,markers = "SOX2",min_expr = 1,max_expr = 2)
p3 <- myfeature_plot(ST.seurat,markers = "POU5F1",min_expr = 4,max_expr = 5)
cowplot::plot_grid(p1,p2,p3,ncol = 3)
ggsave(paste0(plot_outdir, "MRT_", min_gene_num, ".EPI_marker_featureplot.png"))

## TE marker
p1 <- myfeature_plot(ST.seurat,markers = "GATA3",min_expr = 2,max_expr = 3)
p2 <- myfeature_plot(ST.seurat,markers = "GATA2",min_expr = 2,max_expr = 3)
p3 <- myfeature_plot(ST.seurat,markers = "CDX2",min_expr = 2,max_expr = 3)
cowplot::plot_grid(p1,p2,p3,ncol = 3)
ggsave(paste0(plot_outdir, "MRT_", min_gene_num, ".TE_marker_featureplot.png"))

## PE marker
p1 <- myfeature_plot(ST.seurat,markers = "GATA4",min_expr = 1,max_expr = 3)
p2 <- myfeature_plot(ST.seurat,markers = "SOX17",min_expr = 2,max_expr = 3)
p3 <- myfeature_plot(ST.seurat,markers = "FOXA2",min_expr = 1,max_expr = 3)
cowplot::plot_grid(p1,p2,p3,ncol = 3)
ggsave(paste0(plot_outdir, "MRT_", min_gene_num, ".PE_marker_featureplot.png"))
# dev.off()

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

ggsave(file=paste0(plot_outdir, "MRT_", min_gene_num, ".marker_dotplot.png"))

library(SCopeLoomR)
library(SCENIC)
library(Rtsne)
library(pheatmap)

##-------------------------fig.2.b-------------------------
load(paste0(plot_outdir, "MRT_", min_gene_num, ".seurat.rda"))
write.csv(t(as.matrix(ST.seurat@assays$RNA@counts)), 
          file = paste0(plot_outdir, "MRT_", min_gene_num, ".rawcount.csv"))
#以python代码
# import loompy as lp
# import numpy as np
# import scanpy as sc
# x=sc.read_csv("../../code/sample.csv")
# row_attrs = {
#   "Gene": np.array(x.var_names),
# }
# col_attrs = {
#   "CellID": np.array(x.obs_names)
# }
# lp.create("MRT_4000.loom",x.X.transpose(),row_attrs,
#           col_attrs)
scenicLoomPath <- paste0(plot_outdir, "MRT_", min_gene_num, ".auc_mtx.loom")
loom <- open_loom(scenicLoomPath)
regulonsAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulon_res <- getAUC(regulonsAUC)
set.seed(0)
tsne_out <- Rtsne(t(regulon_res),
                  dims = 2,pca =T,
                  perplexity = 10, 
                  theta = 0,
                  check_duplicates = FALSE)
rownames(tsne_out$Y) <- colnames(regulon_res)
ST.seurat[["scenic_tSNE_1"]] <- tsne_out$Y[,1]
ST.seurat[["scenic_tSNE_2"]] <- tsne_out$Y[,2]

## tsne
ratio <- diff(range(ST.seurat$scenic_tSNE_1))/diff(range(ST.seurat$scenic_tSNE_2))

## lineage
p1 <- ggplot(data = ST.seurat@meta.data, 
       aes(x=scenic_tSNE_1, y=scenic_tSNE_2,color=Lineage))+
  geom_point(size=2)+
  scale_color_manual(values = MRT_colorlist$Lineage_col)+
  coord_equal(ratio = ratio)+
  labs(x="tSNE_1", y="tSNE_2")+
  theme_bw(base_size = 20)+
  theme(legend.position = "top",
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave(file=paste0(plot_outdir, "MRT_", min_gene_num, ".pyscenic_tsne_lineage.png"))

## group
p2 <- ggplot(data = ST.seurat@meta.data, 
       aes(x=scenic_tSNE_1, y=scenic_tSNE_2,color=Group))+
  geom_point(size=2)+
  scale_color_manual(values = MRT_colorlist$Group_col)+
  coord_equal(ratio = ratio)+
  labs(x="tSNE_1", y="tSNE_2")+
  theme_bw(base_size = 20)+
  theme(legend.position = "top",
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave(file=paste0(plot_outdir, "MRT_", min_gene_num, ".pyscenic_tsne_group.pdf"))

cowplot::plot_grid(p1,p2, ncol = 2)

##-------------------------fig.2.c-------------------------
anno_col <- data.frame(Group=ST.seurat$Group,
                       Lineage=ST.seurat$Lineage)
rownames(anno_col) <- rownames(ST.seurat@meta.data)

anno_colors <- list(
  Lineage=MRT_colorlist$Lineage_col,
  Group=MRT_colorlist$Group_col
)

out.hclust <-pheatmap::pheatmap(regulon_res,
                                cluster_rows = F, cluster_cols = T,
                                annotation_col = anno_col,
                                show_colnames = F,show_rownames = F,
                                annotation_colors = anno_colors,
                                annotation_names_row=F,
                                annotation_names_col = T,
                                scale = "none",silent = F,
                                clustering_method = "ward.D2")

pdf(paste0(plot_outdir, "MRT_", min_gene_num, ".regulon_circos.pdf"),
    height = 7, width = 7)
circos_hclust(info = ST.seurat@meta.data, 
              hclust.out = out.hclust, 
              col_list = MRT_colorlist,
              levels= c("Group", "Lineage"))
dev.off()

##  infer gender
X_linked_genelist <- read.table("F:/Project/3P/reference/X_linked_genelist.txt", 
                                header = T, stringsAsFactors = F)
Y_linked_genelist <- read.table("F:/Project/3P/reference/Y_linked_genelist.txt", 
                                header = T, stringsAsFactors = F)
##-------------------------sup.fig.2.b-------------------------
## infer embryonic sex
ST.seurat[["Gender"]] <- "unkown"
ST.seurat[["Mean_X"]] <- 0
ST.seurat[["Mean_Y"]] <- 0
X_gene <- intersect(X_linked_genelist$Gene_id, rownames(ST.seurat))
Y_gene <- intersect(Y_linked_genelist$Gene_id, rownames(ST.seurat))
for (i in names(table(ST.seurat$Embryo))) {
  X_mean <- mean(apply(ST.seurat@assays$RNA@data[X_gene,colnames(ST.seurat)[ST.seurat$Embryo==i]], 2, mean), na.rm = T)
  Y_mean <- mean(apply(ST.seurat@assays$RNA@data[Y_gene,colnames(ST.seurat)[ST.seurat$Embryo==i]], 2, mean), na.rm = T)
  if(X_mean/Y_mean > 2){
    ST.seurat$Gender[ST.seurat$Embryo==i] <- "Female"
  } else {
    ST.seurat$Gender[ST.seurat$Embryo==i] <- "Male"
  }
  ST.seurat$Mean_X[ST.seurat$Embryo==i] <- X_mean
  ST.seurat$Mean_Y[ST.seurat$Embryo==i] <- Y_mean
}

# table(unique(ST.seurat@meta.data[,c("Gender", "Embryo")])[,c("Embryo", "Gender")])
# table(ST.seurat@meta.data[,c("Embryo", "Gender", 'Group')]) %>% as.data.frame.array() %>% t()
tmp <- unique(ST.seurat@meta.data[, c("Embryo", "Gender","Group", "Mean_X", "Mean_Y")])

library(dplyr)
tmp <- tmp %>% group_by(Gender)
tmp <- tmp %>% arrange(desc(Group),desc(Embryo), .by_group = T) %>% as.data.frame()
tmp$Embryo <- factor(tmp$Embryo, levels = rev(tmp$Embryo))

ggplot(data = tmp, aes(x=Embryo))+
  geom_point(aes(y=Mean_X), color=MRT_colorlist$Gender_col["Female"], size=3)+
  geom_point(aes(y=Mean_Y), color=MRT_colorlist$Gender_col["Male"], size=3)+
  labs(y="Expression level")+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file = paste0(plot_outdir, "_", min_gene_num, ".gender_infer.pdf"),
       width = 12,height = 4)


##-------------------------sup.fig.2.a-------------------------
## volcano plot
library(ggrepel)
# pdf(paste0(plot_outdir, "MRT_", min_gene_num, ".DEG_volcano.pdf"), width = 7,height = 7)
plot_list <- list()
for(lineage in c("TE", "EPI", "PE")){
  #gender <- "Female"
  for(gender in c("", "Female","Male")){
    cell_vector <- ST.seurat$Lineage==lineage & 
      ST.seurat$Gender !=gender
    if(gender==""){
      plot_name <- paste0(lineage,"_bothsex")
    } else if(gender=="Female"){
      plot_name <- paste0(lineage, "_Male")
    } else if(gender=="Male"){
      plot_name <- paste0(lineage, "_Female")
    }
    if(length(table((ST.seurat[,cell_vector]$Group)))==1) next
    group_deg <- FindMarkers(
      ST.seurat[,cell_vector],
      group.by = "Group",
      ident.1 = "ST",
      ident.2 = "ICSI", 
      min.pct = 0.25, 
      test.use = "wilcox"
    )
    write.table(group_deg, file = paste0(plot_outdir, "_", min_gene_num, ".DEG_",plot_name,".txt"),
                quote = F, sep = "\t", col.names = T, row.names = T)
    group_deg$`log10(adj.pvalue)` <- -log10(group_deg$p_val_adj)
    col_vec <- rep("grey", dim(group_deg)[1])
    col_vec[group_deg$avg_log2FC > 1 & group_deg$p_val_adj  < 0.01 ] <- "red"
    col_vec[group_deg$avg_log2FC < -1 & group_deg$p_val_adj  < 0.01 ] <- "blue"
    group_deg$color <- col_vec
    group_deg$label <- rownames(group_deg)
    group_deg[group_deg$color=="grey","label"] <- ""
    # ratio <- diff(range(group_deg$avg_log2FC))/diff(range(group_deg$`log10(adj.pvalue)`))
    
    ## volcano plot
    plot_list[[plot_name]] <- ggplot()+
      geom_point(data = group_deg, aes(x=avg_log2FC, y=`log10(adj.pvalue)`, color=color))+
      scale_color_manual(values = c(blue="blue", grey="grey", red="red"))+
      geom_hline(yintercept = 2, linetype="longdash")+
      geom_vline(xintercept = c(-1, 1), linetype="longdash")+
      # geom_text_repel(data = group_deg, aes(x=avg_log2FC, 
      #                                       y=`log10(adj.pvalue)`, 
      #                                       color=color,
      #                                       label=group_deg[,"label"]),
      #                 xlim = c(-max(abs(range(group_deg$avg_log2FC))),
      #                          max(abs(range(group_deg$avg_log2FC)))),
      #                 fontface="italic",
      #                 color="black",
      #                 max.overlaps = 10)+
      labs(x="Log2(fold change)", 
           title = plot_name,
           y="-Log10(ajusted p value)")+
      xlim(-max(abs(range(group_deg$avg_log2FC))),
           max(abs(range(group_deg$avg_log2FC))))+
      theme_bw(base_size = 20)+
      theme(#axis.text = element_text(size = 15),
        #axis.title = element_text(size = 15),
        legend.position = "none",
        plot.title = element_text(vjust = 1, hjust = .5))
    # print(p)
    plot_outdir <- "F:/Project/3P/plot/MRT"
    # ggsave(paste(plot_outdir, min_gene_num,plot_name, "DEG_volcano.png", sep = "."), width = 6, height = 7)
    # n <- n+1
    
    out_deg <- group_deg[group_deg$color!="grey",]
    if(nrow(out_deg)==0) next
    else out_deg$gene <- rownames(out_deg)
    # write.table(out_deg, file = paste0(plot_outdir,"_", min_gene_num, ".DEG_",plot_name,".txt"),
    #             quote = F, sep = "\t", col.names = T, row.names = F)
  }
}
# dev.off()

library(cowplot)
do.call("plot_grid", plot_list)

##-------------------------fig.2.a-------------------------
## linear regression
data <- as.data.frame(t(FetchData(object = ST.seurat, slot = "data",vars = rownames(ST.seurat))))

# pdf(file = paste0(plot_outdir, "MRT_", min_gene_num, ".regression_pearson.pdf"))

plot_list <- list()
for(lineage in c("TE", "EPI", "PE")){
  cor_df <- data.frame(ST=apply(data[,ST.seurat$Lineage==lineage & 
                                       ST.seurat$Group=="ST"], 1, mean),
                       ICSI=apply(data[,ST.seurat$Lineage==lineage & 
                                         ST.seurat$Group=="ICSI"], 1, mean))
  model.lm <- lm(ICSI ~ ST, data = cor_df)
  summary(model.lm)
  
  cor_df$Gene_id <- rownames(cor_df)
  deg_df <- read.table(
    paste0(plot_outdir, "_", min_gene_num, ".DEG_",lineage,"_bothsex.txt"), 
    header = T,
    row.names = 1, 
    stringsAsFactors = F
  )
  deg_df$Gene_id <- rownames(deg_df)
  # colnames(deg_df) <- sub("gene", "Gene_id", colnames(deg_df))
  cor_df <- merge(cor_df, deg_df, by = "Gene_id", all.x = T)
  
  col_vec <- rep("grey", dim(cor_df)[1])
  col_vec[cor_df$avg_log2FC > 1 & cor_df$p_val_adj  < 0.01 ] <- "red"
  col_vec[cor_df$avg_log2FC < -1 & cor_df$p_val_adj  < 0.01 ] <- "blue"
  cor_df$color <- col_vec
  cor_df$label <- cor_df$Gene_id
  cor_df[cor_df$color=="grey","label"] <- ""
  
  cor_ratio <- diff(range(cor_df$ICSI, na.rm = T))/diff(range(cor_df$ST, na.rm = T))
  plot_list[[lineage]] <- ggplot(data = cor_df, aes(x=ICSI, y=ST, color=color))+
    geom_point()+
    scale_colour_manual(values = c("grey"="black", "red"="red", "blue"="blue"))+
    annotate("text", label = paste0("Up:",as.numeric(table(cor_df$color)["red"])), 
             x = 1, y = 5, size = 10, colour = "red")+
    annotate("text", label = paste0("Down:",as.numeric(table(cor_df$color)["blue"])), 
             x = 5, y = 1, size = 10, colour = "blue")+
    coord_equal(ratio = cor_ratio)+
    labs(title = paste0(lineage, ":","R-squared=", 
                        format(summary(model.lm)$r.squared, digits = 4)))+
    theme_bw(base_size = 25)+
    theme(panel.border = element_rect(size = 1),
          #panel.grid = element_blank(),
          plot.title = element_text(hjust = .5),
          axis.ticks = element_line(size = 1),
          legend.key.size=unit(1,'cm'),
          legend.position = "none")
  # print(p)
  # ggsave(file = paste0(plot_outdir, "_",min_gene_num, ".",lineage,".regression_pearson.png"))
}
# dev.off()

do.call("plot_grid", c(plot_list, ncol=3))

##-------------------------fig.2.b-------------------------
write.csv(t(as.matrix(ST.seurat@assays$RNA@counts)), 
          file = paste0(plot_outdir, "MRT_", min_gene_num, ".rawcount.csv"))
#以python代码
# import loompy as lp
# import numpy as np
# import scanpy as sc
# x=sc.read_csv("../../code/sample.csv")
# row_attrs = {
#   "Gene": np.array(x.var_names),
# }
# col_attrs = {
#   "CellID": np.array(x.obs_names)
# }
# lp.create("MRT_4000.loom",x.X.transpose(),row_attrs,
#           col_attrs)
scenicLoomPath <- paste0(plot_outdir, "MRT_", min_gene_num, ".auc_mtx.loom")
loom <- open_loom(scenicLoomPath)
regulonsAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulon_res <- getAUC(regulonsAUC)
tsne_out <- Rtsne(t(regulon_res),dims = 2,pca =T,
                  perplexity = 10, theta = 0,check_duplicates = FALSE)
rownames(tsne_out$Y) <- colnames(regulon_res)
ST.seurat[["scenic_tSNE_1"]] <- tsne_out$Y[,1]
ST.seurat[["scenic_tSNE_2"]] <- tsne_out$Y[,2]

## tsne
ratio <- (range(ST.seurat$scenic_tSNE_1)[2]-range(ST.seurat$scenic_tSNE_1)[1])/
  (range(ST.seurat$scenic_tSNE_2)[2]-range(ST.seurat$scenic_tSNE_2)[1])

## lineage
ggplot(data = ST.seurat@meta.data, 
       aes(x=scenic_tSNE_1, y=scenic_tSNE_2,color=Lineage))+
  geom_point(size=2)+
  scale_color_manual(values = MRT_colorlist$Lineage_col)+
  coord_equal(ratio = ratio)+
  theme_bw(base_size = 25)+theme(legend.position = "top",
                                 axis.ticks = element_blank(),
                                 axis.text = element_blank())
ggsave(file=paste0(plot_outdir, "MRT_", min_gene_num, ".pyscenic_tsne_lineage.pdf"))

## group
ggplot(data = ST.seurat@meta.data, 
       aes(x=scenic_tSNE_1, y=scenic_tSNE_2,color=Group))+
  geom_point(size=2)+
  scale_color_manual(values = MRT_colorlist$Group_col)+
  coord_equal(ratio = ratio)+
  theme_bw(base_size = 25)+theme(legend.position = "top",
                                 axis.ticks = element_blank(),
                                 axis.text = element_blank())
ggsave(file=paste0(plot_outdir, "MRT_", min_gene_num, ".pyscenic_tsne_group.pdf"))

##-------------------------fig.2.c-------------------------
out.hclust <-pheatmap::pheatmap(regulon_res,
                                cluster_rows = F, cluster_cols = T,
                                annotation_col = anno_col,
                                show_colnames = F,show_rownames = F,
                                annotation_colors = anno_colors,
                                annotation_names_row=F,
                                annotation_names_col = T,
                                scale = "none",silent = F,
                                clustering_method = "ward.D2")

pdf(paste0(plot_outdir, "_", min_gene_num, ".regulon_circos.pdf"),
    height = 7, width = 7)
circos_hclust(info = ST.seurat@meta.data, 
              hclust.out = out.hclust, 
              col_list = MRT_colorlist,
              levels= c("Group", "Lineage"))
dev.off()


save(ST.seurat, file = paste0(plot_outdir, "_", min_gene_num, ".seurat.rda"))



##-------------------------fig.4.a-------------------------
## merge infercnv output
# dir <- "F:/Project/3P/data/plot/0421/"
# infer_result_ref <- read.table(paste0(dir, "infercnv_MRT_", min_gene_num,"/infercnv.references.txt"), 
#                                header = T, row.names = 1, stringsAsFactors = F)
# infer_result_obs <- read.table(paste0(dir, "infercnv_MRT_", min_gene_num,"/infercnv.observations.txt"), header = T, row.names = 1, stringsAsFactors = F)
# infer_result_total <- cbind(infer_result_obs, infer_result_ref)
# rm(infer_result_ref)
# rm(infer_result_obs)
infer_result_ref <- read.table("F:/Project/3P/plot/infercnv.references.txt",header = T, row.names = 1, stringsAsFactors = F)
infer_result_obs <- read.table("F:/Project/3P/plot/infercnv.observations.txt", header = T, row.names = 1, stringsAsFactors = F)
infer_result_total <- cbind(infer_result_obs, infer_result_ref)
rm(infer_result_ref)
rm(infer_result_obs)

## gene_metadata
setwd(dir)
gene_bed <- read.table("F:/Project/3P/plot/infercnv.gene_order.txt", header = F, stringsAsFactors = F)
colnames(gene_bed) <- c("gene_symbol", "chr", "start", "end")
rownames(gene_bed) <- gene_bed$gene_symbol
gene_bed <- gene_bed[rownames(infer_result_total),]
gene_bed <- gene_bed %>% group_by(chr)
gene_bed$chr <- factor(gene_bed$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrMT"))
gene_bed <- gene_bed %>% arrange(start, .by_group=T)
gene_bed <- as.data.frame(gene_bed)
gene_bed <- gene_bed[!gene_bed$chr %in% c("chrX", "chrY"),]
infer_result <- infer_result_total[gene_bed$gene_symbol,]

## cell_metadata
embryo_list <- unique(ST.seurat$Embryo)

inferCNV_stat <- stat_CNV(
  infer_result_total, 
  embryo_list = embryo_list, 
  cell_info = ST.seurat@meta.data
)
# inferCNV_stat$embryo["E11","chr8"]<- -1
# inferCNV_stat$embryo["E45",c("chr6", "chr8", "chr15")]<- 1
# inferCNV_stat$sc[info$Cell_id[info$Embryo=="E11"],"chr8"] <- -1
# inferCNV_stat$sc[info$Cell_id[info$Embryo=="E45"],c("chr6", "chr8", "chr15")] <- 1
inferCNV_stat$sc[info$Cell_id[info$Embryo=="E6"],"chr4"] <- -1
inferCNV_stat$sc[info$Cell_id[info$Embryo=="E6"],"chr9"] <- 1

## arrange chromosome by cell
inferCNV_sc <- inferCNV_stat$sc

inferCNV_sc[inferCNV_sc==0] <- 2
inferCNV_sc[inferCNV_sc==-1] <- 1.5
inferCNV_sc$Order <- apply(inferCNV_sc, 1, function(x){min(c(1:22)[which(x!=2)])})
inferCNV_sc <- inferCNV_sc %>% arrange(Order)
info_tmp <- unique(ST.seurat@meta.data[,c("Cell_id", "Lineage", "Group")])
inferCNV_sc$Cell_id <- rownames(inferCNV_sc)
inferCNV_sc <- merge(inferCNV_sc, info_tmp, by = "Cell_id")
inferCNV_sc <- inferCNV_sc %>% group_by(Group, Lineage) %>% arrange(Order, .by_group=T) %>% as.data.frame()
inferCNV_sc$Lineage <- factor(inferCNV_sc$Lineage, levels = c("EPI", "PE", "TE"))
inferCNV_sc$Group <- factor(inferCNV_sc$Group, levels = c("ICSI", "ST"))


## plot heatmap
chr_col <- rep(c("black", "white"),11)
names(chr_col) <- paste0("chr", 1:22)
col_anno <- HeatmapAnnotation(
  chr=paste0("chr", 1:22),
  col = list(chr=chr_col),
  show_annotation_name = F,
  border = T,
  show_legend = F
  )
row_anno <- rowAnnotation(
  Lineage=inferCNV_sc$Lineage,
  Group=inferCNV_sc$Group,
  col = list(Lineage=MRT_colorlist$Lineage_col,
             Group=MRT_colorlist$Group_col),
  show_annotation_name=F
)
ht <- Heatmap(
  as.matrix(inferCNV_sc[,2:23]),
  name="CNV",
  show_row_names = F,
  show_column_names = F,
  left_annotation = row_anno,
  top_annotation = col_anno,
  row_split = inferCNV_sc[,c("Lineage", "Group")],
  row_gap = unit(0, "mm"),
  row_title_rot = 0,
  cluster_rows=F,
  border = T,
  show_row_dend = F,
  cluster_columns=FALSE,
  column_split = factor(paste0("chr", 1:22), levels = paste0("chr", 1:22)),
  column_gap = unit(0, "mm"),
  column_title_rot = 90,
  cluster_row_slices = F,
  col=c("1"="red","2"="white", "1.5"="blue")
)
draw(ht)


##-------------------------fig.4.b-------------------------
ST.seurat[["inferCNV_CNV_count"]] <- 0
ST.seurat@meta.data[rownames(inferCNV_stat$sc),"inferCNV_CNV_count"] <- rowSums(abs(inferCNV_stat$sc))
ST.seurat[["inferCNV_ploidy_sc"]] <- "euploidy"
ST.seurat$inferCNV_ploidy_sc[ST.seurat$inferCNV_CNV_count!=0] <- "aneuploidy"

ploidy_count_embryo_per <- ST.seurat@meta.data %>% 
  dplyr::select(Group, Embryo,inferCNV_ploidy_sc) %>% 
  dplyr::group_by(Group, Embryo, inferCNV_ploidy_sc) %>%
  dplyr::summarise(Ploidy=n()) %>% 
  dplyr::mutate(total_cell=sum(Ploidy),CNV_freq=round(Ploidy/sum(Ploidy)*100,2)) %>%
  dplyr::filter(total_cell >=30) %>%
  as.data.frame()

write.table(ploidy_count_embryo_per, 
            file = paste0(outdir, "MRT_", min_gene_num, "ploidy_count_embryo_per.txt"),
            quote = F, sep = "\t",
            col.names = T, row.names = F)

ggplot(data = ploidy_count_embryo_per %>% 
         filter(inferCNV_ploidy_sc=="aneuploidy"), 
       aes(x=Group, y=CNV_freq, fill=Group))+
  geom_boxplot(width=.6, outlier.shape = NA)+
  labs(y="CNV frequency")+
  scale_fill_manual(values = MRT_colorlist$Group_col)+
  geom_signif(comparisons = list(c("ST", "ICSI")),
              test = "t.test",
              textsize = 10,
              y_position = 65)+
  theme_bw(base_size = 25)+
  theme(axis.title.x = element_blank())




##-------------------------fig.4.b-------------------------
