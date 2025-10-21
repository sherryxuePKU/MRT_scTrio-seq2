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
                            min_cell=3, min_gene=MIN_GENE_NUM,
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

stat_CNV <-function(x){
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
    select_col <- ST.seurat$Cell_id[ST.seurat$Embryo==j]
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