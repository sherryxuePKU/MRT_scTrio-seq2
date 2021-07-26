#!/datc/houyu/software/miniconda3/envs/R4.0/bin/Rscript
arg <- commandArgs(T)


#indir <- "/datd/houyu/project/MRT/Trio_DNA/03.methylkit"
#sample <- "/datd/houyu/project/MRT/Trio_DNA/ST0103_sample.xls"
#outdir <- "/datd/houyu/project/MRT/Trio_DNA/04.Info/ST1115/ST1115"
#region <- "CGI"
dir <- arg[1]
batch <- arg[2]
#outdir <- arg[3]
region <- arg[3]

# .libPaths("/datc/houyu/software/miniconda3/envs/R4.0/lib/R/library")
# library(ggplot2)
#Group_color <- c(ST="#E74632",NC="#5068A8", Ref="#FFB923")

indir <- paste0(dir, "/03.methylkit")
outdir <- paste0(dir, "/04.info/",batch, "/", batch, ".meth_", region, ".absoluteDistance.avg.txt") #ST1227.meth_CGI.absoluteDistance.avg.txt
sample <- paste0(dir, "/", batch, "_sample.xls")

sample <- read.table(sample, header = F, stringsAsFactors = F)
sample <- sample$V1
#info_df <- read.table(info_dir, header = T, stringsAsFactors = F)

file_list <- list()
for(i in sample){
  filename <- paste0(indir, "/", i, "/", i, "_CpG.", region,".absoluteDistance.txt")
  file_list[[i]] <- read.table(file = filename, header = F, stringsAsFactors = F)
  file_list[[i]] <- unique(file_list[[i]])
  file_list[[i]] <- file_list[[i]][,-1]
  file_list[[i]] <- colMeans(na.rm = T, file_list[[i]])
  file_list[[i]] <- data.frame(file_list[[i]])
  file_list[[i]]$Coord <- 1:nrow(file_list[[i]])
  file_list[[i]]$Pos <- as.character(1:nrow(file_list[[i]]))
  file_list[[i]]$Pos <- factor(file_list[[i]]$Pos, levels = file_list[[i]]$Pos)
  file_list[[i]]$Sample <- i
}

merge_file <- do.call("rbind", file_list)
#merge_file <- merge(tmp, info_df[info_df$Lineage=="TE",c("Sample", "Group")], by = "Sample")
#merge_file$Color <- "white"
#merge_file[merge_file$Group=="ST","Color"] <- Group_color["ST"]
#merge_file[merge_file$Group=="NC","Color"] <- Group_color["NC"]
#cols <- merge_file$Color
#names(cols) <- merge_file$Sample
#colnames(merge_file) <- c("Sample","data", "Coord", "Pos","Group", "Color")
colnames(merge_file) <- c("data", "Coord", "Pos","Sample")
write.table(merge_file, file = outdir, quote = F, sep = "\t", col.names = T, row.names = F)


# pdf(paste0(outdir, ".meth_",region,".pdf"))
# #ggplot(data = merge_file, aes(x=Coord, y=data, color=Sample))+geom_line()+scale_x_continuous(breaks = c(25, 125), labels = c("TSS", "TES"))+scale_color_manual(values = cols)+theme_bw()+theme(panel.grid.minor.x = element_blank(), axis.title = element_blank(), legend.position = "none")
# if(region=="genebody"){
#   ggplot(data = merge_file, aes(x=Coord, y=data, color=Sample))+geom_line()+scale_x_continuous(breaks = c(0,25, 125,150), labels = c("-15kb","TSS", "TES", "+15kb"))+ labs(y="DNA methylation level (%)")+theme_bw()+theme(panel.grid.minor.x = element_blank(), axis.title = element_blank(), legend.position = "none")
# } else if(region=="CGI"){
#   ggplot(data = merge_file, aes(x=Coord, y=data, color=Sample))+geom_line()+scale_x_continuous(breaks = c(0,25, 125,150), labels = c("-15kb","Start", "End", "+15kb"))+ labs(y="DNA methylation level (%)")+theme_bw()+theme(panel.grid.minor.x = element_blank(), axis.title = element_blank(), legend.position = "none")
# }
# dev.off()

