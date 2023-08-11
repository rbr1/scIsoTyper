#### Run on Oxford compute
srun -p short --cpus-per-task 5 --pty bash
module purge
module load HDF5/1.10.5-gompi-2019a
module load umap-learn/0.3.10-foss-2019a-Python-3.7.2
module load Seurat/3.1.2-foss-2019a-R-3.6.0
module load Harmony/1.0.0-foss-2019a-R-3.6.0
R

library("Seurat")
library('harmony') 
library(ggplot2)
library(pryr)
library(future) 

library(ggplot2)
library(BiocParallel)
library(org.Hs.eg.db)

.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/users/immune-rep/kvi236/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos='http://cran.ma.imperial.ac.uk/')

concat = function(v) {
  res = ""
  for (i in 1:length(v)){res = paste0(res,v[i])}
  res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha)) }

a = 1
if(a==1){
  file="/users/immune-rep/mfj169/SINGLE_CELL_RNA_SEQ/10X_PIPELINE/Samples_PDAC150K_2.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,"To_use_in_PDAC150K"]=="Yes"),]
  sample_id = as.character(p[,"Sample_Name"])
  sample_output_id = as.character(p[,"Sample_Name"])
  BCR.location = as.character(p[,"Location_of_BCR"])
  BCR.location = gsub("all_contig","filtered_contig",BCR.location)
  TCR.location = as.character(p[,"Location_of_TCR"])
  TCR.location = gsub("all_contig","filtered_contig",TCR.location)
  Overall_sample_group = as.character(p[,"Patient"])
  Site = as.character(p[,"Sample_type"])
  batch = "PDAC150Ka"
  
}
PDAC150K_dir = "/gpfs3/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/PLOTS/"
Peng_dir = "/well/immune-rep/shared/10X_GENOMICS/PDAC_CHINESE_DATA/PLOTS/"
Steel_dir = "/well/immune-rep/shared/10X_GENOMICS/PDAC_GSE155698/Plots/"
Peng_dir2 = "/well/immune-rep/users/mfj169/10X_GENOMICS/2019_PUBLISHED_PDAC/PLOTS/"
out_dir = "/gpfs2/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/COMBINED_WITH_PUBLISHED/"
batch = "COMBINED_PDAC150K"

## First get the datasets that you would like to merge for label transfer. 
## We recommend merging related cell types only (such as myeloids, B cells, T cells, non-immune cells separately) to ensure best batch correction. 
## Here, we use harmony to perform that batch correction between datasets, and this is exemplified on the myeloid cells in this subroutine. 

Get_datasets_and_merge <- function(out_dir, PDAC150K_dir, Peng_dir, Steel_dir, Peng_dir2, batch){
  file = concat(c(PDAC150K_dir,"Overall_UMAP_annotations_PDAC150Ka_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[grep("Myeloid", p[,"cell_types"]),]
  
  # Peng all old annotations
  file = concat(c(Peng_dir, "Overall_UMAP_annotations_2019_PUBLISHED_PDAC_Myeloid_cells.txt"))
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  # Steele all old annotations
  file = concat(c(Steel_dir, "Overall_UMAP_annotations_PDAC_GSE155698_all.txt"))
  p2 <- as.matrix(read.csv(file, head=T, sep="\t"))
  p2=p2[which( p2[,3] %in% c("DC", "MGK","granulocyte","mast","momac","monocyte")),]
  
  ###################################### PDAC150K
  type = "Myeloid_cells"
  analysis = "Merge_PDAC150K_Peng_Steele"
  data_1 = readRDS(file = concat(c(PDAC150K_dir,"Seurat_harmonised_Myeloid_cells_clusters_classified_PDAC150Ka_23_12_2021_inc_ILCs.pbmc")))
  head(data_1@meta.data)
  table(data_1@meta.data$filtered_cell_type)
  
  length(intersect(rownames(p), rownames(data_1@meta.data)))
  length(rownames(p))
  length(rownames(data_1@meta.data))
  
  inter1 = intersect(rownames(p), rownames(data_1@meta.data))
  data_1 <- subset(x = data_1, cells = inter1 )
  
  cell_ids =rownames(data_1@meta.data)
  Patient = data_1@meta.data$Patient
  Sample.Type = data_1@meta.data$Sample.Type
  overall_cell_type = data_1@meta.data$filtered_cell_type
  names(overall_cell_type) = cell_ids
  names(Patient) = cell_ids
  names(Sample.Type) = cell_ids
  count_data1= data_1@assays$ RNA@ counts
  
  #################################### Peng
  data_2 = readRDS(file=concat(c(Peng_dir2,"/Seurat_UMAP_2019_PUBLISHED_PDAC_21_09_2021_myeloid_cell.rds")))
  head(data_2@meta.data)
  
  length(intersect(rownames(p1), rownames(data_2@meta.data)))
  length(rownames(p1))
  length(rownames(data_2@meta.data))
  
  cell_ids =rownames(data_2@meta.data)
  Patient2 = data_2@meta.data$orig.ident
  Sample.Type2 = data_2@meta.data$Histology
  names(Patient2) = cell_ids
  names(Sample.Type2) = cell_ids
  count_data2= data_2@assays$ RNA@ counts
  
  #################################### Steele
  data_3= readRDS(file=concat(c(Steel_dir, "Seurat_harmony_Myeloid_annotions2_PDAC_GSE155698.pbmc")))
  head(data_3@meta.data)
  
  length(intersect(rownames(p2), rownames(data_3@meta.data)))
  length(rownames(p2))
  length(rownames(data_3@meta.data))
  
  inter3 = intersect(rownames(p2), rownames(data_3@meta.data))
  data_3 <- subset(x = data_3, cells = inter3 )
  
  cell_ids =rownames(data_3@meta.data)
  Patient3 = data_3@meta.data$orig.ident
  Sample.Type3 = Patient3
  Sample.Type3[grep("TISSUE", Sample.Type3)] = "Tumour"
  Sample.Type3[grep("PBMC", Sample.Type3)] = "PBMC"
  Sample.Type3[grep("AdjNorm_TISSUE", Sample.Type3)] = "Normal"
  names(Patient3) = cell_ids
  names(Sample.Type3) = cell_ids
  count_data3= data_3@assays$ RNA@ counts
  
  #################################### Merge and combine
  PLOTS = ""
  nz1 = apply(count_data1, 1, sum)
  nz2 = apply(count_data2, 1, sum)
  nz3 = apply(count_data3, 1, sum)
  intersect_genes = intersect(intersect(names(which(nz1>0)), names (which(nz2>0))), names(which(nz3>0)))
  
  data1 <- CreateSeuratObject(count_data1[intersect_genes,], project = "PDAC150K")
  data2 <- CreateSeuratObject(count_data2[intersect_genes,], project = "PENG")
  data3 <- CreateSeuratObject(count_data3[intersect_genes,], project = "STEELE")
  
  data1@meta.data$orig.ident = Patient
  data1@meta.data$Sample.Type = Sample.Type
  data1@meta.data$cell_type = overall_cell_type
  data1@meta.data$source = rep("PDAC150K", length(overall_cell_type))
  head(data1@meta.data)
  
  data2@meta.data$orig.ident = Patient2
  data2@meta.data$Sample.Type = Sample.Type2
  data2@meta.data$cell_type = rep("UNKNOWN", length(Patient2))
  data2@meta.data$source = rep("PENG", length(Sample.Type2))
  head(data2@meta.data)
  
  data3@meta.data$orig.ident = Patient3
  data3@meta.data$Sample.Type = Sample.Type3
  data3@meta.data$cell_type = rep("UNKNOWN", length(Patient3))
  data3@meta.data$source = rep("STEELE", length(Sample.Type3))
  head(data3@meta.data)
  
  
  merge<- merge(data1,y = data2)
  merge<- merge(merge,y = data3)
  
  head(merge@meta.data)
  table(merge@meta.data$orig.ident)
  table(merge@meta.data$Sample.Type)
  table(merge@meta.data$cell_type)
  table(merge@meta.data$source)
  
  Sample.Type = merge@meta.data$Sample.Type
  Sample.Type[which(Sample.Type=="biopsy")] = "Tumour"
  Sample.Type[which(Sample.Type=="blood")] = "PBMC"
  
  merge@meta.data$Sample.Type = Sample.Type
  cell_type =merge@meta.data$cell_type
  source=merge@meta.data$source
  names(Sample.Type) = rownames(merge@meta.data)
  names(source) = rownames(merge@meta.data)
  names(cell_type) = rownames(merge@meta.data)
  merge1 = merge
  #################################### Dim red and harmony
  merge1 <- NormalizeData(merge1)
  merge1 <- ScaleData(merge1)
  merge1 <- FindVariableFeatures(merge1, selection.method = "vst", nfeatures = 2000)
  merge1 <- RunPCA(merge1, verbose = FALSE)
  merge1 <- RunHarmony(merge1, c("orig.ident","Sample.Type", "source"),plot_convergence = TRUE, nclust = 50, max.iter.cluster = 100, max.iter.harmony = 10)
  
  merge1 = RunUMAP(merge1 ,reduction = "harmony", dims = 1:30,n.components = 10)
  merge1 = FindNeighbors(merge1 ,reduction = "umap", dims = 1:2)
  merge1 = FindClusters(merge1,reduction = "umap" ,resolution = 1)
  
  head(merge1@meta.data)
  type = "initial_merging"
  library(cowplot)
  p1 <- DimPlot(object = merge1, reduction = "harmony", pt.size = .1, group.by = "source", do.return = TRUE)
  p2 <- VlnPlot(object = merge1, features = "harmony_1", group.by = "source", pt.size = .1)
  p3 <- DimPlot(object = merge1, reduction = "harmony", pt.size = .1, group.by = "cell_type", do.return = TRUE)
  p4 <- DimPlot(object = merge1, reduction = "harmony", pt.size = .1, group.by = "seurat_clusters", do.return = TRUE)
  
  fileout1=concat(c(out_dir,PLOTS,"Seurat_", batch,"_",type,"_1.pdf"))
  w=6
  pdf(file=fileout1, height=w*2, width=w*2)
  par(mfrow= c(1,1), mar = c(4,4,4,4))
  plot_grid(p1,p2,p3,p4)
  dev.off()
  
  fileout1=concat(c(out_dir,PLOTS,"Seurat_", batch,"_",type,"_2.pdf"))
  w=6
  pdf(file=fileout1, height=w*2, width=w*2.5)
  par(mfrow= c(1,1), mar = c(4,4,4,4))
  p1 <- DimPlot(merge1, reduction = "umap", group.by = "Sample.Type", pt.size = .1, do.return = TRUE)
  p2 <- DimPlot(merge1, reduction = "umap", group.by = "cell_type", pt.size = .1, do.return = TRUE)
  p3 <- DimPlot(merge1, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, do.return = TRUE,label = T)
  p4 <- DimPlot(merge1, reduction = "umap", group.by = "source", pt.size = .1, do.return = TRUE)
  plot_grid(p1,p2,p3,p4)
  dev.off()
  
  
  fileout1=concat(c(out_dir,PLOTS,"Seurat_", batch,"_",type,"_3.pdf"))
  w=6
  pdf(file=fileout1, height=w*1, width=w*15)
  par(mfrow= c(1,1), mar = c(4,4,4,4))
  DimPlot(merge1, reduction = "umap", group.by = "Sample.Type", pt.size = .1, split.by = "cell_type")
  dev.off()
  saveRDS(file = concat(c(out_dir,PLOTS,"Seurat_merged_", batch,"_",type,".rds")), merge1) 
  print (concat(c("RDS file saved at: ",out_dir,PLOTS,"Seurat_merged_", batch,"_",type,".rds")))
}

## open Seurat object with merged datasets
type = "Myeloid_cells" ## this will be included in all the output file names
merge1 = readRDS(file = concat(c(out_dir,PLOTS,"Seurat_merged_", batch,"_",type,".rds"))) 

cell_type_label = "cell_type" # this is the column in the Seurat metadata which includes the reference cell annotations, and the query annotations should be labelled with "UNKOWN"
subsampling_level = 300 # this is the maximum number of cells subsampled per cell type for the classifier to be trained on (where more cells are available)
dimred_use = "harmony" # this is the dimensionality reduction to use, either "harmony" or "umap"

Run_SVMCellTransfer<-function(merge1, out_dir,PLOTS, batch,type, cell_type_label = "cell_type", subsampling_level = 300, dimred_use = "harmony"){
  library(e1071)
  cell_ids = rownames(merge1@meta.data)
  cell_type = merge1@meta.data[,cell_type_label]
  cells_pre_annotated = cell_ids[which(cell_type!="UNKNOWN")]
  cells_of_interest_query =setdiff(cell_ids, cells_pre_annotated)
  
  ct = merge1@meta.data[,cell_type_label]
  
  names(ct) = rownames(merge1@meta.data)
  table(ct[cells_of_interest_query])
  ref_metadata=ct[cells_pre_annotated]
  ref_exp_full=merge1@ reductions$ harmony@ cell.embeddings[cells_pre_annotated,]
  dat = data.frame(ref_exp_full, y = as.factor(ref_metadata))
  
  ## only select X cells per classification
  include = NULL
  classes = sort(unique(ref_metadata))
  for(c in c(1:length(classes))){
    w = which(ref_metadata==classes[c])
    if(length(w)>subsampling_level){w = sample(w,subsampling_level)}
    include = c(include, w)
  }
  table(ref_metadata[include])
  dat = dat[include,]
  
  # fit SVM model
  svmfit = svm(y ~ ., data = dat, kernel = "radial", cost = 10, scale = FALSE)
  print("SVM model fitted")
  print(svmfit)
  
  # predict data
  if(dimred_use == "harmony"){
    query_exp_full=merge1@ reductions$ harmony@ cell.embeddings[cells_of_interest_query,]
    exp_full=merge1@ reductions$ harmony@ cell.embeddings
  }else{
    query_exp_full=merge1@ reductions$ umap@ cell.embeddings[cells_of_interest_query,]
    exp_full=merge1@ reductions$ umap@ cell.embeddings
  }
  #query_pred= predict(svmfit, query_exp_full)
  ref_pred= predict(svmfit, ref_exp_full)
  perc_correct = length(which(ref_pred==ref_metadata))*100/length(ref_metadata)
  perc_correct
  table(ref_pred,ref_metadata)
  cell_typ_pred= predict(svmfit, exp_full)
  table(cell_typ_pred)
  
  table(cell_typ_pred, merge1@meta.data$source)
  
  cell_type_overall = cell_type
  names(cell_type_overall) = cell_ids
  cell_type_overall[cells_of_interest_query] = as.character(cell_typ_pred[cells_of_interest_query])
  table(cell_type_overall)
  merge1@meta.data$cell_type_overall = cell_type_overall
  
  ## original cell labels returned
  names(cell_type) = cell_ids
  names(cell_type_overall) = cell_ids
  
  # check 
  table(cell_type_overall[cells_pre_annotated],cell_type[cells_pre_annotated])
  
  merge1@meta.data$refined_annotations = cell_type_overall
  return(merge1)
}
merge1 = Run_SVMCellTransfer(merge1, out_dir,PLOTS, batch,type, cell_type_label = "cell_type", subsampling_level = 300, dimred_use = "harmony")
  

Plot_results<-function(merge1, out_dir,PLOTS, batch,type){

  library(cowplot)
  fileout1=concat(c(out_dir,PLOTS,"Seurat_filtering_", batch,"_",type,"_3.pdf"))
  w=4
  pdf(file=fileout1, height=w*2, width=w*4)
  par(mfrow= c(1,1), mar = c(4,4,4,4))
  p1 <- DimPlot(object =merge1, reduction = "umap", group.by = "source", pt.size = .1, do.return = TRUE,label = T)
  p2 <- DimPlot(object =merge1, reduction = "umap", group.by = "cell_type", pt.size = .1, do.return = TRUE)
  p3 <- DimPlot(object =merge1, reduction = "umap", group.by = "cell_type", pt.size = .1, do.return = TRUE)
  p4 <- DimPlot(object =merge1, reduction = "umap", group.by = "refined_annotations", pt.size = .1, do.return = TRUE, label = F)
  
  plot_grid(p1,p3,p4)
  dev.off()
  
  ########## plot expression of key genes 
  ### get gene lists
  file = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/MYELOID/cluster_annotation_gene_signatures.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  sign = p[,"annot"]
  gene.signature = strsplit(gsub(" ", "",p[,"gene.signature"],fixed = T),",")
  names(gene.signature) = sign
  gene.signature1 = gene.signature
  sign_genes = sort(unique(unlist(gene.signature)))
  genes = rownames(merge1@ assays$RNA@ counts)
  
  file = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/MYELOID/MYELOID_BLOOD/cluster_annotation_gene_signatures.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  sign = p[,"annot"]
  gene.signature = strsplit(gsub(" ", "",p[,"gene.list"],fixed = T),",")
  names(gene.signature) = sign
  gene.signature2 = gene.signature
  sign_genes1 = sort(unique(unlist(gene.signature)))
  sign_genes=sort(intersect(unique(sign_genes, sign_genes1), genes))
  
  genes =rownames(merge1@assays$ RNA@ data)	
  
  genes_of_interest =sort(unique(intersect(genes, sign_genes)))
  exp = merge1@ assays$ RNA@ data
  exp1 = as.matrix(exp[genes_of_interest,])
  cell_types = sort(unique(merge1@meta.data$refined_annotations))
  
  mat_mean = matrix(data=0, nrow = length(cell_types), ncol = length(genes_of_interest), dimnames = c(list(cell_types),list(genes_of_interest))) 
  for(i in c(1:length(cell_types))){
    w = rownames(merge1@meta.data)[which(merge1@meta.data$refined_annotations == cell_types[i])]
    mat_mean [cell_types[i],] = apply(exp1[,w], 1, mean)
  }
  mat_mean1= mat_mean[,which(apply(mat_mean, 2, sd)!=0)]
  mat_mean1  = mat_mean1 [,c(which(apply(mat_mean1, 2, max)>1,which(colnames(mat_mean1) %in% c("AICDA"))))]
  
  hc = hclust(dist(t(mat_mean1)))
  labels = hc$ labels[hc$order] 
  hc = hclust(dist(mat_mean1))
  labels1 = hc$ labels[hc$order] 
  
  Idents(object = merge1) = merge1@meta.data$refined_annotations
  levels(x = merge1) = sort(unique(merge1@meta.data$refined_annotations))
  library(RColorBrewer)
  cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
  fileout1=concat(c(out_dir,PLOTS,"Seurat_filtering_", batch,"_",type,"_4.pdf"))
  w=4
  pdf(file=fileout1, height=w*1.9, width=w*8)
  par(mfrow= c(1,1), mar = c(5,5,3,3))
  DotPlot(merge1, assay = "RNA", features = labels) +
    scale_color_viridis_c() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
  
  xy = merge1@ reductions$ umap@ cell.embeddings[,c(1:2)]
  cell_type = merge1@meta.data$refined_annotations
  table(merge1@meta.data$refined_annotations, merge1@meta.data$orig.ident)
  
  orig.ident = merge1@meta.data$orig.ident
  w = grep("CD45p", orig.ident)
  orig.ident[w] = apply(cbind(merge1@meta.data$orig.ident[w],merge1@meta.data$Sample.Type[w]),1,paste,collapse = "_")
  orig.ident=gsub("_CD45p2_Tumour","_biopsy",gsub("_CD45p1_Tumour","_biopsy",orig.ident))
  orig.ident=gsub("_CD45p2_PBMC","_blood",gsub("_CD45p1_PBMC","_blood",orig.ident))
  sort(unique(orig.ident))
  
  out = cbind(xy, apply(cbind("Myeloid" , cell_type), 1, paste, collapse = " "),orig.ident, merge1@meta.data$source, merge1@meta.data$Sample.Type)
  colnames(out) = c("UMAP1", "UMAP2","cell type","sample","source","sample_type")
  out_file_table=concat(c(out_dir,PLOTS,"Overall_UMAP_", batch,"_",type,".txt"))
  write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  
  t = table(merge1@meta.data$refined_annotations, merge1@meta.data$orig.ident)
  out_file_table = concat(c(out_dir,PLOTS,"/Overall_cell_counts_", batch,"_",type,".txt"))
  write.table(t(t), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")

}
Plot_results(merge1, out_dir,PLOTS, batch,type)

saveRDS(file = concat(c(out_dir,PLOTS,"Seurat_merged_SVMCellTranfer_annotated_", batch,"_",type,".rds")), merge1) 

