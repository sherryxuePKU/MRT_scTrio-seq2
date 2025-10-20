arg <- commandArgs(T)

dir <- arg[1]
batch <-arg[2]

C_type <- "CpG"
#project <- arg[6]
input <- paste0(dir, "/03.methylkit")
output <- paste0(dir,"/04.info/",batch,"/",batch,".CpG.site_anno.xls")

sample <- paste0(dir, "/", batch, "_sample.xls")
sample <- read.table(sample, stringsAsFactors = F)
sample <- sample$V1

subgroup_dir <- "~/tangfuchou_coe/xuexiaohui/database/hg38/Annotation/subgroup.txt"
subgroup <- read.table(subgroup_dir, stringsAsFactors = F)
subgroup <- sapply(strsplit(subgroup$V1, "[.]"), "[[",2)
#sample <- dir(path = input)

N <- length(sample)
M <- length(subgroup)
fileName <- list()
#annolist <- list()
annodf_CG <- list()
#annodf_CH <- list()
#annodf_C <- list()
for (i in 1:N) {
  annodf_CG[[i]] <- data.frame(sample=sample[i])
  #annodf_CH[[i]] <- data.frame(sample=sample[i])
  #annodf_C[[i]] <- data.frame(sample=sample[i])
  for (j in 1:M) {
    fileName <- paste0(input, "/",sample[i], "/",sample[i],".",subgroup[j],"_",C_type, ".bed")
    #print(fileName)
    anno_tmp<- read.table(fileName, header = F, stringsAsFactors = F)
    colnames(anno_tmp) <-  c("Chr","Start","End","Strand","CpG_total","Freq_C")
    annodf_CG[[i]][,subgroup[j]] <- sum(anno_tmp[,"CpG_total"]*anno_tmp[,"Freq_C"])/sum(anno_tmp[,"CpG_total"])
    #annodf_CH[[i]][,subgroup[j]] <- sum(anno_tmp[, "CH_meth"])/sum(anno_tmp[, "CH_total"])
    #annodf_C[[i]][,subgroup[j]] <- sum(anno_tmp[, c("CpG_meth","CH_meth")])/sum(anno_tmp[, c("CpG_total","CH_total")])
  }
}
merged_CG <- do.call("rbind", annodf_CG)
#merged_CH <- do.call("rbind", annodf_CH)
#merged_C <- do.call("rbind", annodf_C)

write.table(merged_CG, file = output, row.names = F, col.names = T, quote = F, sep = "\t")

