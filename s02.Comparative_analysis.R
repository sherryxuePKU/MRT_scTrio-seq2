## load packages
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(SCopeLoomR)
library(SCENIC)
library(Rtsne)
library(pheatmap)

## check packages' version
packageVersion("Seurat")

## state paths
save_seurat_dir <- paste0("F:/Project/3P/R_workspace/ST_", Sys.Date(), "_seurat.rda")
plot_outdir <- "F:/Project/3P/plot/"

## define function
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

## load input
load(save_seurat_dir)
X_linked_genelist <- read.table("F:/Project/3P/reference/X_linked_genelist.txt", 
                                header = T, stringsAsFactors = F)
Y_linked_genelist <- read.table("F:/Project/3P/reference/data/Y_linked_genelist.txt", 
                                header = T, stringsAsFactors = F)

## other parameter
min_gene_num <- 4000

## sup.fig.2.a
## volcano plot
pdf(paste0(plot_outdir, "MRT_", min_gene_num, ".DEG_volcano.pdf"), width = 7,height = 7)
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
    group_deg <- FindMarkers(ST.seurat[,cell_vector], 
                             ident.1 = "ST",ident.2 = "ICSI", 
                             min.pct = 0.25, test.use = "wilcox")
    write.table(group_deg, file = paste0(outdir, "MRT_", min_gene_num, ".DEG_",plot_name,".txt"), 
                quote = F, sep = "\t", col.names = T, row.names = T)
    group_deg$`log10(adj.pvalue)` <- -log10(group_deg$p_val_adj)
    col_vec <- rep("grey", dim(group_deg)[1])
    col_vec[group_deg$avg_log2FC > 1 & group_deg$p_val_adj  < 0.01 ] <- "red"
    col_vec[group_deg$avg_log2FC < -1 & group_deg$p_val_adj  < 0.01 ] <- "blue"
    group_deg$color <- col_vec
    group_deg$label <- rownames(group_deg)
    group_deg[group_deg$color=="grey","label"] <- ""
    ratio <- (range(group_deg$avg_log2FC)[2]-range(group_deg$avg_log2FC)[1])/
      (range(group_deg$`log10(adj.pvalue)`)[2]-range(group_deg$`log10(adj.pvalue)`)[1])
    
    ## volcano plot
    p<- ggplot(data = group_deg, aes(x=avg_log2FC, y=`log10(adj.pvalue)`, color=color))+
      geom_point()+
      scale_color_manual(values = c(blue="blue", grey="grey", red="red"))+
      geom_hline(yintercept = 2, linetype="longdash")+
      geom_vline(xintercept = c(-1, 1), linetype="longdash")+
      geom_text_repel(aes(label=group_deg[,"label"]),
                      xlim = c(-max(abs(range(group_deg$avg_log2FC))),
                               max(abs(range(group_deg$avg_log2FC)))),
                      fontface="italic",
                      color="black",
                      max.overlaps = 10)+
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
    #print(p)
  }
}
dev.off()

## sup.fig.2.c
## heatmap
pdf(paste0(plot_outdir, "MRT_", min_gene_num, ".DEG_heatmap.pdf"), width = 8,height = 5)
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
    group_deg <- FindMarkers(ST.seurat[,cell_vector], 
                             ident.1 = "ST",ident.2 = "ICSI", 
                             min.pct = 0.25, test.use = "wilcox")
    write.table(group_deg, file = paste0(outdir, "MRT_", min_gene_num, ".",lineage,"_DEG_bothsex.txt"), 
                quote = F, sep = "\t", col.names = T, row.names = T)
    group_deg$`log10(adj.pvalue)` <- -log10(group_deg$p_val_adj)
    col_vec <- rep("grey", dim(group_deg)[1])
    col_vec[group_deg$avg_log2FC > 1 & group_deg$p_val_adj  < 0.01 ] <- "red"
    col_vec[group_deg$avg_log2FC < -1 & group_deg$p_val_adj  < 0.01 ] <- "blue"
    group_deg$color <- col_vec
    group_deg$label <- rownames(group_deg)
    group_deg[group_deg$color=="grey","label"] <- ""
    ratio <- (range(group_deg$avg_log2FC)[2]-range(group_deg$avg_log2FC)[1])/
      (range(group_deg$`log10(adj.pvalue)`)[2]-range(group_deg$`log10(adj.pvalue)`)[1])
    
    
    ## heatmap
    selected_col <- colnames(ST.seurat)[cell_vector]
    #selected_row <- rownames(group_deg)[group_deg$color!="grey"]
    selected_row <- group_deg[group_deg$color!="grey","label"]
    selected_row <- c("FMR1NB","HIST1H1E","HIST1H1B","ACTG2","MAGEB2","VCX2",
                      "MRPL2","HSPA1B", "HSPA1A", "RPS4Y1", "EIF1AY")
    if(length(selected_row)==0) next()
    expr_scale <- ST.seurat[["RNA"]]@scale.data[selected_row,selected_col]
    
    limit <- 1.5
    expr_scale[expr_scale> limit] <- limit
    expr_scale[expr_scale< -limit] <- (-limit)
    
    anno_col <- data.frame(Gender=ST.seurat$Gender,
                           Group=ST.seurat$Group,
                           #Embryo=ST.seurat$Embryo, 
                           Lineage=ST.seurat$Lineage)
    rownames(anno_col) <- rownames(ST.seurat@meta.data)
    
    anno_colors <- list(
      Lineage=MRT_colorlist$Lineage_col,
      Group=MRT_colorlist$Group_col,
      Gender=MRT_colorlist$Gender_col
    )
    
    
    anno_col2 <- anno_col[selected_col,]
    anno_col2$cell_id <- rownames(anno_col2)
    anno_col2 <- anno_col2 %>% group_by(Group)
    anno_col2 <- anno_col2 %>% arrange(Gender, .by_group=T)
    
    col_anno <- HeatmapAnnotation(Group=anno_col2$Group,
                                  Gender = anno_col2$Gender,
                                  col = list(Group=MRT_colorlist$Group_col,
                                             Gender=MRT_colorlist$Gender_col),
                                  border = F,
                                  show_annotation_name = T,
                                  show_legend = T)
    
    colors_func <- colorRamp2(c(min(expr_scale, na.rm = T),
                                mean(expr_scale, na.rm = T),
                                max(expr_scale, na.rm = T)),
                              c("#91BFDB", "white", "#D73027"))
    
    ht <- ComplexHeatmap::Heatmap(expr_scale[,anno_col2$cell_id],
                                  show_row_names = T,
                                  show_column_names = F,
                                  row_names_gp = gpar(fontface="italic"),
                                  column_split = anno_col2[,c("Group","Gender")],
                                  column_gap = unit(0, "mm"),
                                  row_split = group_deg[group_deg$color!="grey","color"],
                                  row_gap = unit(0, "mm"),
                                  #column_title_rot = 0,
                                  #column_title_gp = gpar(fontsize = 10),
                                  border = TRUE,
                                  border_gp = gpar(lty = 1, lwd = 1),
                                  cluster_columns = T,
                                  cluster_column_slices = F,
                                  cluster_row = F,
                                  cluster_row_slices = F,
                                  row_order = selected_row,
                                  na_col = "lightgrey",
                                  top_annotation = col_anno,
                                  show_row_dend = T,
                                  show_column_dend = F,
                                  clustering_method_rows = "ward.D2",
                                  clustering_distance_rows = "euclidean",
                                  col = colors_func,
                                  name = "ht",
                                  column_title=plot_name,
                                  #row_title = "Differentially expressed gene",
                                  use_raster = T,
                                  raster_device="CairoPNG",
                                  raster_quality = 5)
    draw(ht)
  }
}
dev.off()

## fig.2.a
## linear regression
data <- as.data.frame(t(FetchData(object = ST.seurat, slot = "data",vars = rownames(ST.seurat))))

pdf(file = paste0(plot_outdir, "MRT_", min_gene_num, ".regression_pearson.pdf"))
for(lineage in c("TE", "EPI", "PE")){
  cor_df <- data.frame(ST=apply(data[,ST.seurat$Lineage==lineage & 
                                       ST.seurat$Group=="ST"], 1, mean),
                       ICSI=apply(data[,ST.seurat$Lineage==lineage & 
                                         ST.seurat$Group=="ICSI"], 1, mean))
  model.lm <- lm(ICSI ~ ST, data = cor_df)
  summary(model.lm)
  
  cor_df$Gene_id <- rownames(cor_df)
  deg_df <- read.table(paste0(outdir, "MRT_", min_gene_num, ".DEG_",lineage,"_bothsex.txt"), 
                       row.names = 1,header = T, stringsAsFactors = F)
  deg_df$Gene_id <- rownames(deg_df)
  cor_df <- merge(cor_df, deg_df, by = "Gene_id", all.x = T)
  
  col_vec <- rep("grey", dim(cor_df)[1])
  col_vec[cor_df$avg_log2FC > 1 & cor_df$p_val_adj  < 0.01 ] <- "red"
  col_vec[cor_df$avg_log2FC < -1 & cor_df$p_val_adj  < 0.01 ] <- "blue"
  cor_df$color <- col_vec
  cor_df$label <- cor_df$Gene_id
  cor_df[cor_df$color=="grey","label"] <- ""
  
  cor_ratio <- (range(cor_df$ICSI)[2]-range(cor_df$ICSI)[1])/
    (range(cor_df$ST)[2]-range(cor_df$ST)[1])
  p <- ggplot(data = cor_df, aes(x=ICSI, y=ST, color=color))+
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
  print(p)
}
dev.off()

## sup.fig.2.b
## infer embryonic sex
seurat.obj <- ST.seurat
seurat.obj[["Gender"]] <- "unkown"
seurat.obj[["Mean_X"]] <- 0
seurat.obj[["Mean_Y"]] <- 0
X_gene <- intersect(X_linked_genelist$Gene_id, rownames(seurat.obj))
Y_gene <- intersect(Y_linked_genelist$Gene_id, rownames(seurat.obj))
for (i in names(table(seurat.obj$Embryo))) {
  X_mean <- mean(apply(seurat.obj@assays$RNA@data[X_gene,colnames(seurat.obj)[seurat.obj$Embryo==i]], 2, mean), na.rm = T)
  Y_mean <- mean(apply(seurat.obj@assays$RNA@data[Y_gene,colnames(seurat.obj)[seurat.obj$Embryo==i]], 2, mean), na.rm = T)
  if(X_mean/Y_mean > 2){
    seurat.obj$Gender[seurat.obj$Embryo==i] <- "Female"
  } else {
    seurat.obj$Gender[seurat.obj$Embryo==i] <- "Male"
  }
  seurat.obj$Mean_X[seurat.obj$Embryo==i] <- X_mean
  seurat.obj$Mean_Y[seurat.obj$Embryo==i] <- Y_mean
}

# table(unique(ST.seurat@meta.data[,c("Gender", "Group", "Embryo")])[,c("Gender", "Group")])
# tmp <- unique(ST.seurat@meta.data[, c("Embryo", "Gender","Group", "Mean_X", "Mean_Y")])

tmp <- tmp %>% group_by(Gender)
tmp <- tmp %>% arrange(desc(Group),desc(Embryo), .by_group = T) %>% as.data.frame()
tmp$Embryo <- factor(tmp$Embryo, levels = rev(tmp$Embryo))

ggplot(data = tmp, aes(x=Embryo))+
  geom_point(aes(y=Mean_X), color=MRT_colorlist$Gender_col["Female"], size=3)+
  geom_point(aes(y=Mean_Y), color=MRT_colorlist$Gender_col["Male"], size=3)+
  labs(y="Expression level")+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file = paste0(plot_outdir, "MRT_", min_gene_num, ".gender_infer.pdf"),
       width = 12,height = 4)

## fig.2.b
write.csv(t(as.matrix(ST.seurat@assays$RNA@counts)), 
          file = paste0(plot_outdir, "MRT_", min_gene_num, ".rawcount.csv"))
#以下为python代码
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
scenicLoomPath <- paste0(outdir, "MRT_", min_gene_num, ".auc_mtx.loom")
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

## fig.2.c
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
