#!/datc/houyu/software/miniconda3/envs/R4.0/bin/Rscript
arg <- commandArgs(T)

dir <- arg[1]
batch <- arg[2]
#outdir <- arg[3]
region <- arg[3]

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
colnames(merge_file) <- c("data", "Coord", "Pos","Sample")
write.table(merge_file, file = outdir, quote = F, sep = "\t", col.names = T, row.names = F)

