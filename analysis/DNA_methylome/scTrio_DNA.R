##-------------------------load packages-------------------------
library(Seurat)
library(ggplot2)
library(data.table)
library(ggtern)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(scales)
library(ggridges)
library(reshape2)
library(dplyr)

## check packages' version
packageVersion("Seurat")

##-------------------------hyperparameter-------------------------
MIN_GENE_NUM <- 4000

##-------------------------state paths-------------------------
input_dir <- "/path/to/your/data/"  # specify your input path
plot_outdir <- "/path/to/your/plot/"  # specify your output path
seurat_dir <- paste0(input_dir, "seurat/MRT_", MIN_GENE_NUM, ".seurat.rda")
info_dir <- paste0(input_dir, "metadata/MRT_Meth.CellsPass_info.txt")
color_dir <- paste0(input_dir, "metadata/colors.rda")

##-------------------------marker list-------------------------
TE_marker <- c("GATA3", "DAB2", "TFAP2C", "GATA2","CDX2","PTGES",
               "EMP2", "TGFBR3", "PDGFA", "KRT18","CLDN10", "PLAC8", "TRIML1")
PE_marker <- c("PDGFRA", "GATA4", "FOXA2","HNF1B", "COL4A1",  "FN1",
               "LINC00261","AMOTL1", "DPP4","SOX17")
EPI_marker <- c("SOX2",  "NANOG", "TDGF1","GDF3", "PRDM14",
                "NODAL","ARGFX", "DPPA2", "POU5F1")

##-------------------------load colorlist-------------------------
load(color_dir)
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
load(seurat_dir)
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


do.call("plot_grid", plot_list)

##-------------------------fig.3.b------------------------
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
anno_list <- list()
indir <- "/path/to/your/CpG_site_anno/"
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

require(ggpubr)
require(RColorBrewer)
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
  
  genebody_meth_merge_df <- read.table(
    paste0("/path/to/your/MRT_", MIN_GENE_NUM, ".",region, "_meth_merge.txt"), 
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
DMR_dir <- "/path/to/your/DMRs/"
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
                              levels = regions[-c(19,20)]
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
data <- read.table("/path/to/your/MRT.PBAT_FreeC_read.merge.txt",
                   header = T, stringsAsFactors = F)
data <- data[,c("Chr", "Start", merge_df$Sample)]



