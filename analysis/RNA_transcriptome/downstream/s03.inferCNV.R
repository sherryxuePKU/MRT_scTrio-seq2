##-------------------------load packages-------------------------
library(Seurat)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

## check package version
packageVersion("Seurat")
packageVersion("infercnv")

##-------------------------state paths-------------------------
save_seurat_dir <- paste0("/path/to/your/ST_", Sys.Date(), "_seurat.rda")
outdir <- "/path/to/your/plot/"

##-------------------------hyperparameters-------------------------
MIN_GENE_NUM <- 4000

##-------------------------define function-------------------------
source("analysis/RNA_transcriptome/helper_function.R")

##-------------------------load input-------------------------
load(save_seurat_dir)
gene_bed <- read.table("/path/to/your/hg38.gencode.p5.allGene.bed", stringsAsFactors = F)
ST_info <- read.table("/path/to/your/MRT_RNA.CellsPass_info.txt", header = T, stringsAsFactors = F)

## prepare input for infercnv
setwd(outdir)

# raw_count_mtx
raw_count_mtx.dir <- as.data.frame(t(FetchData(object = ST.seurat,slot = "counts",vars = rownames(ST.seurat))))
write.table(raw_count_mtx.dir, file = "infercnv.umi_count.txt", quote = F, sep = "\t", col.names = T, row.names = T)

# annotation
anno_dir <- ST.seurat@meta.data[,c("Cell_id", "Group", "Lineage")]
anno_dir$Group <- paste0(anno_dir$Group, "_",anno_dir$Lineage)
anno_dir <- anno_dir[,c("Cell_id", "Group")]
write.table(anno_dir, file = "infercnv.annotation.txt", quote = F, sep = "\t", col.names = F, row.names = F)

# gene_order
colnames(gene_bed) <- c("chr", "start", "end", "gene_symbol")
gene_bed <- gene_bed[gene_bed$gene_symbol %in% rownames(ST.seurat), c(4,1,2,3)]
gene_bed <- gene_bed[!duplicated(gene_bed$gene_symbol),]
gene_bed <- gene_bed[gene_bed$chr %in% paste0("chr",c(1:22, "X", "Y")),]
write.table(gene_bed, file = "infercnv.gene_order.txt",quote = F, sep = "\t", col.names = F, row.names = F)

## run infercnv
dir <- outdir
raw_count_mtx.dir <- paste0(dir, "infercnv.umi_count.txt")
anno_dir <- paste0(dir, "infercnv.annotation.txt")
gene_order_dir <- paste0(dir, "infercnv.gene_order.txt")
ref_group <- paste0("ICSI_", c("TE", "EPI", "PE"))
out_dir <- paste0(dir, "infercnv_MRT_", MIN_GENE_NUM)
dir.create(out_dir)
chr_exclude_list <-  c("chrGL000218.1", "chrGL000219.1", "chrKI270711.1", "chrKI270721.1" ,
                       "chrKI270728.1", "chrKI270734.1", "chrMT", "chrX", "chrY")

ST.infercnv <- CreateInfercnvObject(raw_counts_matrix=raw_count_mtx.dir,
                                    annotations_file=anno_dir,
                                    gene_order_file=gene_order_dir,
                                    ref_group_names=ref_group,
                                    max_cells_per_group = NULL,
                                    min_max_counts_per_cell = NULL,
                                    chr_exclude = chr_exclude_list) 

ST.infercnv <- infercnv::run(ST.infercnv,output_format="pdf",
                             cutoff=1, min_cells_per_gene = 3,
                             window_length = 101,out_dir=out_dir,
                             cluster_by_groups=TRUE)

## plot infercnv output
infer_result_ref <- read.table(paste0(dir, "infercnv_MRT_", 
                                      MIN_GENE_NUM,"/infercnv.references.txt"), 
                               header = T, row.names = 1, stringsAsFactors = F)
infer_result_obs <- read.table(paste0(dir, "infercnv_MRT_", 
                                      MIN_GENE_NUM,"/infercnv.observations.txt"), header = T, row.names = 1, stringsAsFactors = F)
infer_result_total <- cbind(infer_result_obs, infer_result_ref)
rm(infer_result_ref)
rm(infer_result_obs)

# gene_metadata
gene_bed <- read.table("infercnv.gene_order.txt", header = F, stringsAsFactors = F)
colnames(gene_bed) <- c("gene_symbol", "chr", "start", "end")
rownames(gene_bed) <- gene_bed$gene_symbol
gene_bed <- gene_bed[rownames(infer_result_total),]
gene_bed <- gene_bed %>% group_by(chr)
gene_bed$chr <- factor(gene_bed$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrMT"))
gene_bed <- gene_bed %>% arrange(start, .by_group=T)
gene_bed <- as.data.frame(gene_bed)
gene_bed <- gene_bed[!gene_bed$chr %in% c("chrX", "chrY"),]
infer_result <- infer_result_total[gene_bed$gene_symbol,]

# cell_metadata
rownames(ST_info) <- ST_info$Cell_id
info <- ST_info[ST_info$Cell_id %in% colnames(infer_result_total),]
rownames(info) <- info$Cell_id

chromosome <- paste("chr",c(1:22),sep="")
embryo_list <- unique(ST.seurat$Embryo)

# correct infercnv result by inferred CNV from PBAT-seq data
inferCNV_stat <- stat_CNV(infer_result_total)
inferCNV_stat$embryo["E11","chr8"]<- -1
inferCNV_stat$embryo["E45",c("chr6", "chr8", "chr15")]<- 1
inferCNV_stat$sc[info$Cell_id[info$Embryo=="E11"],"chr8"] <- -1
inferCNV_stat$sc[info$Cell_id[info$Embryo=="E45"],c("chr6", "chr8", "chr15")] <- 1
inferCNV_stat$sc[info$Cell_id[info$Embryo=="E6"],"chr4"] <- -1
inferCNV_stat$sc[info$Cell_id[info$Embryo=="E6"],"chr9"] <- 1

##-------------------------fig.4.a-------------------------
inferCNV_sc <- inferCNV_stat$sc
inferCNV_sc[inferCNV_sc==0] <- 2
inferCNV_sc[inferCNV_sc==-1] <- 1.5
inferCNV_sc$Order <- apply(inferCNV_sc, 1, function(x){min(c(1:22)[which(x!=2)])})
inferCNV_sc <- inferCNV_sc %>% arrange(Order)
info_tmp <- unique(info[,c("Cell_id", "Lineage", "Group")])
inferCNV_sc$Cell_id <- rownames(inferCNV_sc)
inferCNV_sc <- merge(inferCNV_sc, info_tmp, by = "Cell_id")
inferCNV_sc <- inferCNV_sc %>% group_by(Group, Lineage) %>% arrange(Order, .by_group=T) %>% as.data.frame()
inferCNV_sc$Lineage <- factor(inferCNV_sc$Lineage, levels = c("EPI", "PE", "TE"))
inferCNV_sc$Group <- factor(inferCNV_sc$Group, levels = c("ICSI", "ST"))

inferCNV_sc <- inferCNV_sc %>% arrange(chr1, chr2, chr3, chr4, chr5, chr6,
                                       chr7, chr8, chr9, chr10, chr11, chr12,
                                       chr13, chr14, chr15, chr16, chr17, chr18,
                                       chr19, chr20, chr21, chr22)

chr_col <- rep(c("black", "white"),11)
names(chr_col) <- paste0("chr", 1:22)
col_anno <- HeatmapAnnotation(chr=paste0("chr", 1:22),
                              col = list(chr=chr_col),
                              show_annotation_name = F,
                              border = T,
                              show_legend = F)
row_anno <- rowAnnotation(Lineage=inferCNV_sc$Lineage,
                          Group=inferCNV_sc$Group,
                          col = list(Lineage=MRT_colorlist$Lineage_col,
                                     Group=MRT_colorlist$Group_col),
                          show_annotation_name=F)
ht <- Heatmap(as.matrix(inferCNV_sc[,2:23]),
              name="CNV",border = T,
              show_row_names = F,show_column_names = F,
              left_annotation = row_anno,top_annotation = col_anno,
              row_split = inferCNV_sc[,c("Lineage", "Group")],
              row_gap = unit(0, "mm"),
              row_title_rot = 0,column_title_rot = 90,
              cluster_rows=F,cluster_columns=F,
              column_split = factor(paste0("chr", 1:22), levels = paste0("chr", 1:22)),
              column_gap = unit(0, "mm"),
              cluster_row_slices = F,
              col=c("1"="red","2"="white", "1.5"="blue"))

pdf(paste0(outdir, "MRT_", MIN_GENE_NUM, ".inferCNV_sc_heatmap.pdf"),
    height = 7, width = 9)
draw(ht)
dev.off()


##-------------------------sup.fig.6.a-------------------------
chr_col <- rep(c("black", "white"),11)
names(chr_col) <- paste0("chr", 1:22)
col_anno <- HeatmapAnnotation(chr=paste0("chr", 1:22),
                              col = list(chr=chr_col),
                              show_annotation_name = F,
                              show_legend = F,
                              border = T)
row_anno <- rowAnnotation(Lineage=inferCNV_sc$Lineage[inferCNV_sc$Cell_id %in% info$Cell_id[info$Embryo=="E34"]],
                          Group=inferCNV_sc$Group[inferCNV_sc$Cell_id %in%info$Cell_id[info$Embryo=="E34"]],
                          col = list(Lineage=MRT_colorlist$Lineage_col,Group=MRT_colorlist$Group_col),
                          show_annotation_name=F)
ht <- Heatmap(as.matrix(inferCNV_sc[inferCNV_sc$Cell_id %in% info$Cell_id[info$Embryo=="E34"],2:23]),
              name="CNV",
              show_row_names = F,show_column_names = F,
              left_annotation = row_anno,top_annotation = col_anno,
              row_split = inferCNV_sc[inferCNV_sc$Cell_id %in% 
                                        info$Cell_id[info$Embryo=="E34"],
                                      c("Lineage", "Group")],
              row_gap = unit(0, "mm"),
              row_title_rot = 0,
              cluster_rows=F,cluster_columns=F,
              border = T,show_row_dend = F,
              column_split = factor(paste0("chr", 1:22), 
                                    levels = paste0("chr", 1:22)),
              column_gap = unit(0, "mm"),
              column_title_rot = 90,
              col=c("1"="red","2"="white", "1.5"="blue"))

pdf(paste0(outdir, "MRT_", MIN_GENE_NUM, "_E34.inferCNV_sc_heatmap.pdf"),
    height = 2, width = 10)
draw(ht)
dev.off()

##-------------------------fig.4.b & fig.4.c-------------------------
chromosome <- paste("chr",c(1:22),sep="")
embryo_list <- unique(ST.seurat$Embryo)
ST.seurat[["inferCNV_CNV_count"]] <- 0
ST.seurat@meta.data[rownames(inferCNV_stat$sc),"inferCNV_CNV_count"] <- rowSums(abs(inferCNV_stat$sc))
ST.seurat[["inferCNV_ploidy_sc"]] <- "euploidy"
ST.seurat$inferCNV_ploidy_sc[ST.seurat$inferCNV_CNV_count!=0] <- "aneuploidy"

ploidy_count_group <- ST.seurat@meta.data %>% 
  dplyr::select(Group, inferCNV_ploidy_sc) %>% 
  dplyr::group_by(Group, inferCNV_ploidy_sc) %>%
  dplyr::summarise(Ploidy=n()) %>% 
  dplyr::mutate(CNV_freq=round(Ploidy/sum(Ploidy)*100,2)) %>%
  as.data.frame()

ST.seurat@meta.data %>%
  dplyr::select(Group, Embryo,inferCNV_ploidy_sc) %>% 
  dplyr::group_by(Group, Embryo) %>%
  dplyr::summarise(nEmbryo=n())

ploidy_count_embryo_per <- ST.seurat@meta.data %>% 
  dplyr::select(Group, Embryo,inferCNV_ploidy_sc) %>% 
  dplyr::group_by(Group, Embryo, inferCNV_ploidy_sc) %>%
  dplyr::summarise(Ploidy=n()) %>% 
  dplyr::mutate(total_cell=sum(Ploidy),CNV_freq=round(Ploidy/sum(Ploidy)*100,2)) %>%
  dplyr::filter(total_cell >=30) %>%
  as.data.frame()

write.table(ploidy_count_embryo_per, 
            file = paste0(outdir, "MRT_", MIN_GENE_NUM, "ploidy_count_embryo_per.txt"),
            quote = F, sep = "\t",col.names = T, row.names = F)

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
ggsave(paste0(outdir, "MRT_", MIN_GENE_NUM, ".ploidy_count_embryo_per.pdf"), 
       width = 7, height = 7)

