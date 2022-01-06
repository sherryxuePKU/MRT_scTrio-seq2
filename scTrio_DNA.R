##-------------------------load packages-------------------------
library(Seurat)
library(ggplot2)
library(data.table)
library(ggtern)
library(ComplexHeatmap)
library(circlize)

## check packages' version
packageVersion("Seurat")

##-------------------------other parameter-------------------------
min_gene_num <- 4000

##-------------------------state paths-------------------------
plot_outdir <- "F:/Project/3P/plot/MRT"
seurat_dir <- paste0(plot_outdir, "_", min_gene_num, ".seurat.rda")
info_dir <- paste0(plot_outdir, "Meth_Total_SampleInfo.txt")
load(seurat_dir)

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

myfeature_plot <- function(ST.seurat,markers,
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

plot_lineage_ternary_func <- function(ST.seurat=ST.seurat,
                                      marker_dir,lineage=c("EPI", "PE", "TE"),
                                      group_by="Lineage",
                                      cols=MRT_colorlist[["Lineage_col"]]){
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

calc_coord_ratio <- function(reduction_method="pca",dim=c(1,2),seurat.obj){
  ST.seurat <- seurat.obj
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

##-------------------------load colorlist-------------------------
load("F:/Project/3P/R_workspace/colors.rda")
MRT_colorlist$Group_col["ST"] <- "#f5bd2a"
MRT_colorlist$Group_col["ICSI"] <- "#704d9c"

##-------------------------quality control-------------------------
info_df <- read.table(info_dir, header = T, stringsAsFactors = F)
attach(info_df)
info_df<- info_df[
  CpG_TotalCpG.1X. > 2e+6 & 
    Lambda_percent < 25 &
    #Stage %in% c(4,5) & 
    Conversion_ratio > 99 &
    lambda_percent==10,]
detach()

##-------------------------fig.3.a------------------------
meth.seurat <- subset(ST.seurat, subset = ID %in% info_df$Sample)
meth.seurat <-NormalizeData(meth.seurat, normalization.method = "LogNormalize", scale.factor = 1e+5, verbose = F)
meth.seurat <-FindVariableFeatures(meth.seurat, selection.method = "vst")
meth.seurat <- ScaleData(meth.seurat, features = rownames(meth.seurat), verbose = F)
meth.seurat <- RunPCA(meth.seurat,features = c(PE_marker, EPI_marker,TE_marker))
# DimHeatmap(meth.seurat, dims = 1:10, balanced = TRUE)
# ElbowPlot(meth.seurat)
meth.seurat <- RunTSNE(meth.seurat, dims = 1:5)

reduction <- "tsne"
plot_list <- list()
for (i in c("Group", "Lineage")) {
  p <- DimPlot(meth.seurat, reduction = reduction, 
                pt.size = 2,
                # cols = MRT_colorlist[paste0(i, "_col")],
                group.by = i) + 
    scale_color_manual(values = MRT_colorlist[[paste0(i, "_col")]])+
    theme_bw(base_size = 25) + 
    theme(axis.line = element_line(size = 1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(hjust = .5),
          legend.position = "right")
  plot_list[[i]] <- p + coord_equal(ratio = diff(range(p$data$tSNE_1))/diff(range(p$data$tSNE_2)))
}

library(cowplot)
do.call("plot_grid", plot_list)

##-------------------------fig.3.b------------------------
library(ggpubr)

group_var <- "Lineage"
output_type <- "Ratio"

plot_list <- list()
for(j in c("CpG","CHH", "CHG")){
  for(i in c(1)){
    if(output_type=="Ratio"){
      y <- paste0(j, "_MethRatio.", i, "X.")
      ylab <- "DNA methylation level (%)"
    } else if(output_type=="Total"){
      y <- paste0(j, "_Total",j,".",i,"X.")
      ylab <- "CpG site covered (1X)"
    }
    plot_list[[j]] <- ggboxplot(
       info_df, 
       x = group_var,y = y,
       fill = "Group", 
       # width = .8,
       # outlier.shape = NA,
       palette = MRT_colorlist$Group_col
      )+
      stat_compare_means(
         aes(group = Group), 
         label="p.signif",
         size=5,vjust = 1
        )+
      labs(y=ylab,title=j)+
      theme_bw(base_size = 15)+
      theme(axis.title.x = element_blank(),
            plot.title = element_text(hjust = .5),
            strip.background = element_blank())
    
    # print(p)
  }
}

do.call("plot_grid", c(plot_list, ncol=3))
##-------------------------fig.3.e------------------------
library(scales)
library(ggridges)

## load ref info
ref_info <- read.table("F:/Project/3P/data/REF/2017_NG/2017_NG_StatInfo.txt", 
                       header = T, stringsAsFactors = F)
ref_info$Group <- "Ref"
ref_info$Embryo <- "Ref"
ref_info$Stage <- ""
var_cols <- c("Sample", "Lineage", "Group", "Embryo","Stage")
merge_info <- rbind(merge_df[,var_cols], ref_info[,var_cols])

## load pca input 
pca_input <- read.table("F:/Project/3P/Meth/GC_Merge_PCA_Result_300bp_per0.3.txt",
                        header = T, stringsAsFactors = F, row.names = 1)
pca_input$Sample <- rownames(pca_input)
pca_input <- merge(pca_input, merge_info, by = "Sample")
pca_input <- pca_input[pca_input$Lineage!="Blastocyst",]
pca_input <- pca_input[pca_input$Embryo !="E15",]
pca_input$Pseudotime <- rescale(pca_input$PC1, to=c(0,100)) ## rescale PC1 to 0~100

## rename group to Ref_8-cell, Ref_Morula, Ref_Blastcyst, ICSI, ST
pca_input$Group_1 <- pca_input$Group
selector <- pca_input$Lineage %in% c("8-cell", "Morula")
pca_input$Group_1[selector] <- paste0(pca_input$Group[selector], "_", pca_input$Lineage[selector])
pca_input$Group_1[pca_input$Group_1=="Ref"] <- paste0(pca_input$Group[pca_input$Group_1=="Ref"],"_Blastocyst") 
ref_group <- paste0("Ref_", c("8-cell", "Morula", "Blastocyst"))
pca_input$Group_1 <- factor(pca_input$Group_1, levels = rev(c(ref_group, "ICSI", "ST")))

## ridge plot
ggplot(data = pca_input[pca_input$Embryo !="E15",], 
       aes(x=Pseudotime, y=Group_1, fill=Group_1))+
  geom_density_ridges()+
  scale_fill_manual(values = brewer.pal(8,"Set3")[4:8])+
  theme_bw(base_size = 20)+
  theme(axis.title.y = element_blank(),
        legend.position = "none")



##-------------------------sup.fig.2.b------------------------
meth_summary <- function(x, cells,group, type="ratio", C_type="CpG"){
  if (type=="site") {
    colume <- switch(C_type,
                     CpG=c(1,seq(3,18,5)),
                     CHH=c(1,seq(23,38,5)),
                     CHG=c(1,seq(43,58,5)))
    
  } else if (type=="ratio"){
    colume <- switch(C_type,
                     CpG=c(1,seq(5,20,5)),
                     CHH=c(1,seq(25,40,5)),
                     CHG=c(1,seq(45,60,5)))
  } else {
    stop("Unrecognized type!")
  }
  x <- na.omit(x)
  rownames(x) <- x$Sample
  tmp <- x[cells,colume]
  tmp <- na.omit(tmp)
  #if(class(tmp[,2])=="integer"){
  if(type=="site"){
    tmp[,2:5] <- round(tmp[, 2:5]/1e+06,2)
  }
  tmp_summary <- as.data.frame(t(apply(tmp[,2:5], 2, summary)))
  tmp_summary$Depth <- sapply(strsplit(rownames(tmp_summary), "[.]"), "[[", 2)
  rownames(tmp_summary) <- 1:dim(tmp_summary)[1]
  tmp_summary$Var <- group
  return(tmp_summary)
}

output_type <- "ratio"
plot_list <- list()
# pdf(paste0("Meth_", output_type, "_Group_Depth_lambdaPercent.pdf"),width = 8, height = 7)
for (i in c("CpG","CHH","CHG")) {
  depth_summary_list <- list()
  for (j in c("ICSI", "ST")) {
    for (k in c(10)){
      # selector <- merge_df$Group==j & 
      #   merge_df$Lambda_percent<50 & 
      #   merge_df$lambda_percent==k
      # cells <- merge_df$Sample[selector]
      cells <- merge_df$Sample
      depth_summary_list[[paste0(j,"_",k)]] <- meth_summary(
        stat_meth, 
        cells = cells,
        group = j,
        type=output_type,
        C_type = i
      )
      depth_summary_list[[paste0(j,"_",k)]]$Var <- paste0(depth_summary_list[[paste0(j,"_",k)]]$Var, "_", k)
    }
  }
  
  depth_summary_df<- do.call("rbind", depth_summary_list)
  depth_summary_df$lambda_percent <- sapply(strsplit(depth_summary_df$Var, "_"), "[[", 2)
  depth_summary_df$Var <- sapply(strsplit(depth_summary_df$Var, "_"), "[[", 1)
  depth_summary_df$Depth <- factor(depth_summary_df$Depth, levels = c("1X", "3X", "5X", "10X"))
  
  if(output_type=="ratio"){
    y_title <- "DNA methylation level (%)"
  } else if(output_type=="site"){
    y_title <- "Median of covered sites(M)"
  }
  plot_list[[i]]<- ggplot(data = depth_summary_df, 
                          aes_string(x="Depth",y="Median",color="Var", group="Var"))+
    geom_errorbar(
      aes(ymin=`1st Qu.`, ymax=`3rd Qu.`), 
      width=.05, size=1)+
    geom_line(size=1)+
    scale_color_manual(values = MRT_colorlist$Group_col)+
    scale_x_discrete(label=c(">=1X", ">=3X", ">=5X", ">=10X"))+
    labs(y=y_title,title = i)+
    theme_bw(base_size = 15)+
    theme(legend.title = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(vjust = 1, hjust = .5)
    )
  # print(p)
}
# dev.off()

do.call("plot_grid", c(plot_list, nrow=1))

##-------------------------sup.fig.2.f------------------------
library(reshape2)

anno_list <- list()
indir <- "F:/Project/3P/data/DNA/CpG_site_anno/"
for(i in dir(indir)){
  fileName <- paste0(indir, i)
  anno_list[[i]] <- read.table(fileName,header = T, stringsAsFactors = F)
}
# anno_df<- do.call("rbind", anno_list)
# rownames(anno_df) <- 1:nrow(anno_df)
# write.table(
#   anno_df,
#   file = paste0("F:/Project/3P/plot/MRT_4000.CpG_site_anno.txt"),
#   col.names = T, row.names = F, quote = F, sep = "\t"
# )

info_tmp <- merge_df[,c("Sample", "Group", "Lineage", "Embryo", "Stage", "lambda_percent")]
anno_df<- merge(anno_df, info_tmp, by.x = "sample", by.y="Sample")

anno_long <- reshape2::melt(
  anno_df, 
  id.vars = c("sample", "Group","Lineage", 
              "Embryo", "Stage", "lambda_percent")
)
anno_long$variable <- factor(
  anno_long$variable, 
  levels = c("CGI","Promoter", "HCP", "ICP", "LCP", 
             "Exon", "Intron","Intergenic", "Intragenic",
             "LINE", "LTR", "L1", "L2",
             "SINE", "SVA", "MIR", "ALR", "Alu",
             "ERV1", "ERVK", "ERVL.MaLR","ERVL")
)

library(ggpubr)
library(RColorBrewer)
# pdf("DNA_element_Group_Lineage.pdf", height = 10, width = 15)
ggboxplot(
  anno_long, 
  x = "variable", 
  y = "value",
  fill = "Group",
  outlier.shape = NA, 
  #add = "jitter",
  palette = MRT_colorlist$Group_col
  )+
  facet_grid(Lineage~.)+
  stat_compare_means(
    aes(group = Group), 
       label="p.signif", 
       size=5, vjust = .5
    )+
  labs(y="DNA methylation level(%)")+
  theme_bw(base_size = 15)+
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# dev.off()
##-------------------------sup.fig.3.4------------------------
## prepare input of profile plot
# region  <- "genebody"
for (region in c("genebody", "CGI")) {
  if(region=="CGI"){
    x_label <- c("-15kb","Start", "End", "+15kb")
  }
  if(region=="genebody"){
    x_label <- c("-15kb","TSS", "TES", "+15kb")
  }
  # genebody_meth_list <- list()
  # indir <- paste0("F:/Project/3P/data/DNA/",region, "_profile_mtx/")
  # for(i in dir(indir)){
  #   fileName <- paste0(indir, i)
  #   genebody_meth_list[[i]] <- read.table(fileName,header = T, stringsAsFactors = F)
  # }
  # genebody_meth_merge_df<- do.call("rbind", genebody_meth_list)
  # rownames(genebody_meth_merge_df) <- 1:nrow(genebody_meth_merge_df)
  # write.table(
  #   genebody_meth_merge_df,
  #   file = paste0("F:/Project/3P/plot/MRT_4000.",region, "_meth_merge.txt"),
  #   col.names = T, row.names = F, quote = F, sep = "\t"
  # )
  
  genebody_meth_merge_df <- read.table(
    paste0("F:/Project/3P/plot/MRT_4000.",region, "_meth_merge.txt"), 
    header = T, stringsAsFactors = F
  )
  info_tmp <- info_df[,c("Sample", "Group", "Lineage", "Embryo", "Stage", "lambda_percent")]
  genebody_meth_merge_df<- merge(genebody_meth_merge_df,info_tmp, by="Sample")
  
  ## plot for each embryo, quality control of embryo
  for(i in c("ICSI", "ST")){
    # pdf(paste0(region,"Profile_Embryo_", i, ".pdf"), height = 6, width = 8)
    p <- ggplot(data = genebody_meth_merge_df[genebody_meth_merge_df$Group==i,], 
                aes(x=Coord, y=data, group=Sample, color=Group))+
      geom_line()+scale_x_continuous(breaks = c(0,50, 250,300), 
                                     labels = x_label)+
      labs(y="DNA methylation level (%)")+
      facet_wrap(.~Embryo, ncol = 4)+
      coord_fixed(ratio = 4)+
      #facet_grid(Lineage~Stage)+
      scale_color_manual(values =  MRT_colorlist$Group_col[i])+
      theme_bw(base_size = 15)+
      theme(panel.grid.minor.x = element_blank(), 
            axis.title.x = element_blank(), 
            legend.position = "right",
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
            #strip.text = element_text(size = 15),
            plot.title = element_text(vjust = 1, hjust = .5))
    print(p)
    # dev.off()
  }
}
##-------------------------fig.4.f------------------------
library(dplyr)

genebody_meth_group <- genebody_meth_merge_df %>% 
  select(c("data", "Coord", "Group")) %>% 
  group_by(Group, Coord) %>%
  dplyr::summarise(
    Median=median(data), 
    Q1=quantile(data,0.25), 
    Q3=quantile(data,0.75)
  ) %>%
  mutate(Pos=Coord) 

ggplot(data = genebody_meth_group,aes(x=Coord))+
  geom_ribbon(aes(ymin=Q1, ymax=Q3, fill=Group),
              alpha=.3)+
  geom_line(aes(y=Median, group=Group, color=Group), 
            size=1.5)+
  scale_fill_manual(values = MRT_colorlist$Group_col)+
  scale_color_manual(values = MRT_colorlist$Group_col)+
  scale_x_continuous(breaks = c(0,50, 250,300), labels = x_label)+
  labs(y="DNA methylation level (%)")+
  ylim(c(0,50))+
  coord_fixed(ratio = 3.5)+
  theme_bw(base_size = 25)+
  theme(panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "right",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = .5,vjust = 1),
        #strip.text = element_text(size = 15),
        plot.title = element_text(vjust = 1, hjust = .5))



##-------------------------fig.4.d------------------------
DMR_dir <- "F:/Project//3P/data/DNA/Meth/DMR/data/noCNV/"
regions <- c("CGI","HCP", "ICP", "LCP", 
             "Exon", "Intron","Intergenic", "Intragenic",
             "LINE", "LTR", "L1", "L2",
             "SINE", "SVA", "MIR", "ALR", "Alu",
             "ERV1", "ERVK", "ERVL-MaLR","ERVL", 
             "gene", "promoter_gene")
file.list <- list()
hyper_dmr_count_mtx <- matrix(ncol = 2, nrow = length(regions))
rownames(hyper_dmr_count_mtx) <- regions
hypo_dmr_count_mtx <- matrix(ncol = 2, nrow = length(regions))
rownames(hypo_dmr_count_mtx) <- regions
hyper_DMR_summary <- list()
hypo_DMR_summary <- list()
for (lineage in c("TE", "PE", "EPI")) {
  #lineage <- "PE"
  count_DMR <- TRUE
  # pdf(paste0(DMR_dir, "/plot/",lineage, "_merged.","CpG_correlation.pdf"), width = 7.73, height = 6.41)
  for(i in regions){
    fileName <- paste0(DMR_dir, lineage,"/DMR_merged.", i,"_CpG.bed")# TE_merged.gene_CpG
    file.list[[i]] <- as.data.frame(fread(fileName, header = F, stringsAsFactors = F))
    colnames(file.list[[i]]) <- c("Chr", "Start", "End", "Element", "ST_Meth", "NC_Meth")
    
    mtx <- file.list[[i]]
    total <- nrow(mtx)
    hypo_dmr <- nrow(mtx[(mtx$ST_Meth<10 & mtx$NC_Meth>40),])
    hyper_dmr <- nrow(mtx[(mtx$ST_Meth>40 & mtx$NC_Meth<10),])
    hyper_dmr_count_mtx[i,] <- c(total, hyper_dmr/total*100)
    hypo_dmr_count_mtx[i,] <- c(total, hypo_dmr/total*100)
    
    hyper_DMR_summary[[lineage]] <- as.data.frame(hyper_dmr_count_mtx)
    colnames(hyper_DMR_summary[[lineage]]) <- c("Count", "Freq")
    hypo_DMR_summary[[lineage]] <- as.data.frame(hypo_dmr_count_mtx)
    colnames(hypo_DMR_summary[[lineage]]) <- c("Count", "Freq")
}

## count hyper or hypo DMR 
DMR_summary <- hypo_DMR_summary
DMR_summary <- do.call("rbind", DMR_summary)
DMR_summary <- na.omit(DMR_summary)
DMR_summary$Lineage <- sapply(strsplit(rownames(DMR_summary), "[.]"), "[[", 1)
DMR_summary$Lineage <- factor(DMR_summary$Lineage, levels = c("TE", "EPI", "PE"))
DMR_summary$Element <- sapply(strsplit(rownames(DMR_summary), "[.]"), "[[", 2)
DMR_summary$Element <- factor(DMR_summary$Element, 
                              levels = c("CGI","gene","promoter_gene", "HCP", "ICP", "LCP", 
                                          "Exon", "Intron","Intergenic", "Intragenic",
                                          "LINE", "LTR", "L1", "L2",
                                          "SINE", "SVA", "MIR", "ALR", "Alu",
                                          "ERV1", "ERVK", "ERVL-MaLR","ERVL"
                                         )
                              )

ggplot(data = DMR_summary, aes(x=Element, y=Freq, fill=Lineage))+
  geom_bar(stat = "identity", position = "dodge")+
  coord_flip()+
  ylim(c(0,45))+
  #scale_y_continuous(breaks = seq(0,50,10),labels = seq(0,50,10), limits = seq(0,50,10))+
  scale_fill_manual(values = MRT_colorlist$Lineage_col)+
  theme_bw(base_size = 20)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 1))
}

##-------------------------sup.fig.5.b------------------------
# indir <- "F:/Project/3P/data/DNA/CNV/data/1M/"
# fileNames <- dir(indir)
# fileNames <- fileNames[grep("read", fileNames)]
# data_list <- list()
# for(i in fileNames){
#   data_list[[i]] <- read.table(file = paste0(indir,i), header = T, stringsAsFactors = F)
#   colnames(data_list[[i]]) <- gsub("E34_D7_ST", "E34_D7_Allo_ST", colnames(data_list[[i]]))
#   colnames(data_list[[i]]) <- gsub("E35_D7_ST", "E35_D7_Allo_ST", colnames(data_list[[i]]))
# }
# 
# library(purrr)
# merge2 <- function(x,y)base::merge(x,y,by=c("Chr", "Start"))
# data <- data_list %>% purrr::reduce(merge2)
# write.table(data, file = "F:/Project/3P/data/DNA/CNV/data/MRT.PBAT_FreeC_read.merge.txt", 
#             col.names = T, row.names = F, sep = "\t", quote = F)
data <- read.table("F:/Project/3P/data/DNA/CNV/data/MRT.PBAT_FreeC_read.merge.txt",
                   header = T, stringsAsFactors = F)
data <- data[,c("Chr", "Start", merge_df$Sample)]



