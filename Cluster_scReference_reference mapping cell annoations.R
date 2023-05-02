srun -p short --cpus-per-task 10 --pty bash
module purge
module load HDF5/1.10.5-gompi-2019a
module load umap-learn/0.3.10-foss-2019a-Python-3.7.2
module load Seurat/3.1.2-foss-2019a-R-3.6.0
module load Harmony/1.0.0-foss-2019a-R-3.6.0

R
############### 
library("Seurat")
library('harmony') 
library(ggplot2)
library(pryr)
library(future) 

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



################# get an example dataset (Peng et al myeloid cells) that you can infer cell annotation labels from a pre-labelled dataset (PancrImmune):
# example data is myeloid data but here I have given you some outputs that you can use for the other cell type objects
cell_type_obejct = "myeloid"

Get_example_data<-function(){
  batch = "PDAC150K_COMBINED"
  PLOTS = "/"
  out_dir = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/COMBINED_WITH_PUBLISHED/"
  out_dir_raw = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/COMBINED_WITH_PUBLISHED/"
  
  file = concat(c("/gpfs3/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/PLOTS/Overall_UMAP_annotations_PDAC150Ka_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[grep("Myeloid", p[,"cell_types"]),]
  
  # Peng all old annotations
  file = concat(c("/well/immune-rep/shared/10X_GENOMICS/PDAC_CHINESE_DATA/PLOTS/Overall_UMAP_annotations_2019_PUBLISHED_PDAC_Myeloid_cells.txt"))
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  ###################################### PDAC150K
  type = "Myeloid_cells"
  analysis = "Merge_PDAC150K_Peng_Steele"
  data_1 = readRDS(file = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/PLOTS//Seurat_harmonised_Myeloid_cells_clusters_classified_PDAC150Ka_23_12_2021_inc_ILCs.pbmc")
  head(data_1@meta.data)
  table(data_1@meta.data$filtered_cell_type)
  
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
  out_dir = "/well/immune-rep/users/mfj169/10X_GENOMICS/2019_PUBLISHED_PDAC/"
  PLOTS = "PLOTS/"
  batch = "2019_PUBLISHED_PDAC"
  type = "Myeloid_cells"
  data_2 = readRDS(file=concat(c(out_dir,PLOTS,"/Seurat_UMAP_2019_PUBLISHED_PDAC_21_09_2021_myeloid_cell.rds")))
  head(data_2@meta.data)
  
  cell_ids =rownames(data_2@meta.data)
  Patient2 = data_2@meta.data$orig.ident
  Sample.Type2 = data_2@meta.data$Histology
  names(Patient2) = cell_ids
  names(Sample.Type2) = cell_ids
  count_data2= data_2@assays$ RNA@ counts
  
  #################################### Merge and combine
  out_dir = "/gpfs2/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/COMBINED_WITH_PUBLISHED/"
  batch = "COMBINED_PDAC150K"
  PLOTS = ""
  nz1 = apply(count_data1, 1, sum)
  nz2 = apply(count_data2, 1, sum)
  intersect_genes = intersect(names(which(nz1>0)), names (which(nz2>0)))
  
  data1 <- CreateSeuratObject(count_data1[intersect_genes,], project = "PDAC150K")
  data2 <- CreateSeuratObject(count_data2[intersect_genes,], project = "PENG")

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
  
  merge1<- merge(data1,y = data2)

  head(merge1@meta.data)
  table(merge1@meta.data$orig.ident)
  table(merge1@meta.data$Sample.Type)
  table(merge1@meta.data$cell_type)
  table(merge1@meta.data$source)
  
  #################################### Dimensionality reduction and batch correction via harmony
  merge1 <- NormalizeData(merge1)
  merge1 <- ScaleData(merge1)
  merge1 <- FindVariableFeatures(merge1, selection.method = "vst", nfeatures = 2000)
  merge1 <- RunPCA(merge1, verbose = FALSE)
  merge1 <- RunHarmony(merge1, c("orig.ident","Sample.Type", "source"),plot_convergence = TRUE, nclust = 50, max.iter.cluster = 100, max.iter.harmony = 10)
  
  merge1 = RunUMAP(merge1 ,reduction = "harmony", dims = 1:30,n.components = 10)
  merge1 = FindNeighbors(merge1 ,reduction = "umap", dims = 1:2)
  merge1 = FindClusters(merge1,reduction = "umap" ,resolution = 1)
  
  head(merge1@meta.data)
  return(merge1)
}

object<-Get_example_data()

################# plot the output of merging datasets together to ensure that this is done properly:
## location for output (prefix): please edit to what you need

Plot_merged_output_data<-function(object, fileout_prefix){
  library(cowplot)
  p1 <- DimPlot(object = merge1, reduction = "harmony", pt.size = .1, group.by = "source", do.return = TRUE)
  p2 <- VlnPlot(object = merge1, features = "harmony_1", group.by = "source", pt.size = .1)
  p3 <- DimPlot(object = merge1, reduction = "harmony", pt.size = .1, group.by = "cell_type", do.return = TRUE)
  p4 <- DimPlot(object = merge1, reduction = "harmony", pt.size = .1, group.by = "seurat_clusters", do.return = TRUE)
  
  fileout1=concat(c(fileout_prefix,"_1.pdf"))
  w=6
  pdf(file=fileout1, height=w*2, width=w*2)
  par(mfrow= c(1,1), mar = c(4,4,4,4))
  plot_grid(p1,p2,p3,p4)
  dev.off()
  
  p1 <- DimPlot(merge1, reduction = "umap", group.by = "Sample.Type", pt.size = .1, do.return = TRUE)
  p2 <- DimPlot(merge1, reduction = "umap", group.by = "cell_type", pt.size = .1, do.return = TRUE)
  p3 <- DimPlot(merge1, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, do.return = TRUE,label = T)
  p4 <- DimPlot(merge1, reduction = "umap", group.by = "source", pt.size = .1, do.return = TRUE)
  
  fileout1=concat(c(fileout_prefix,"_2.pdf"))
  w=6
  pdf(file=fileout1, height=w*2, width=w*2.5)
  par(mfrow= c(1,1), mar = c(4,4,4,4))
  plot_grid(p1,p2,p3,p4)
  dev.off()
}

fileout_prefix=concat(c("Example_Seurat_filtering_", batch,"_",cell_type_object,"_"))
Plot_merged_output_data(object, fileout_prefix)
  
################# Reference cell type annotation from pre-batch-corrected data including reference and query cells

# The Seurat object needs to have a column in meta data 
# called "cell_type" in which the reference cell annotations are provided
# and the query cell annotations are given the term "UNKNOWN"

# Batch correction should be done via harmony with the 
# object@ reductions$ harmony@ cell.embeddings as the embeddings object name

# number_cells_sample: the maximum number of reference cells per cell subtype 
# to sample on which to build the model (default = 300)

# The output column will be labelled "scReference_label"

scReference<-function(object, number_cells_sample = 300){
  ### build a classifier on old cells
  library(e1071)
  cell_ids = rownames(object@meta.data)
  cell_type = object@meta.data$cell_type
  cells_pre_annotated = cell_ids[which(cell_type!="UNKNOWN")]
  cells_of_interest_query =setdiff(cell_ids, cells_pre_annotated)
  
  ct = object@meta.data[,"cell_type"]
  
  names(ct) = rownames(object@meta.data)
  table(ct[cells_of_interest_query])
  ref_metadata=ct[cells_pre_annotated]
  ref_exp_full=object@ reductions$ harmony@ cell.embeddings[cells_pre_annotated,]
  dat = data.frame(ref_exp_full, y = as.factor(ref_metadata))
  
  ## only select X cells per classification
  n = 300
  include = NULL
  classes = sort(unique(ref_metadata))
  for(c in c(1:length(classes))){
    w = which(ref_metadata==classes[c])
    if(length(w)>n){w = sample(w,n)}
    include = c(include, w)
  }
  print("subsample levels per cell type:")
  print(table(ref_metadata[include]))
  dat = dat[include,]
  
  svmfit = svm(y ~ ., data = dat, kernel = "radial", cost = 10, scale = FALSE)
  print(svmfit)
  
  # predict data
  query_exp_full=merge1@ reductions$ harmony@ cell.embeddings[cells_of_interest_query,]
  exp_full=merge1@ reductions$ harmony@ cell.embeddings
  cell_typ_pred= predict(svmfit, exp_full)
  
  # check how good the prediction was on the reference data
  library(caret)
  conf = confusionMatrix(dat$y, cell_typ_pred[rownames(dat)])
  print(conf$ overall )
  
  print(table(cell_typ_pred, object@meta.data$cell_type))
  
  ## create output column with predicted labels
  scReference_label = cell_type
  names(scReference_label) = cell_ids
  scReference_label[cells_of_interest_query] = as.character(cell_typ_pred[cells_of_interest_query])
  table(scReference_label)
  object@meta.data$scReference_label = scReference_label
  return(object)
}

scReference(object, number_cells_sample = 300)

saveRDS(file=concat(c("Seurat_scReference_", batch,"_",cell_type_object,".seurat")), object)

################# check output (standard plots)

fileout_prefix=concat(c("Example_Seurat_scReference_", batch,"_",cell_type_object,"_"))

Check_scReference_output<-function(object, fileout_prefix){
  object@meta.data$scReference_label = scReference_label
  ## original cell labels returned
  names(cell_type) = cell_ids
  names(cell_type_overall) = cell_ids
  table(cell_type_overall[cells_pre_annotated],cell_type[cells_pre_annotated])
  
  merge1@meta.data$refined_annotations = cell_type_overall
  
  table(cell_type_overall, merge1@meta.data$source)
  table(cell_type_overall, merge1@meta.data$Sample.Type)
  table(cell_type_overall, merge1@meta.data$seurat_clusters)
  
  library(cowplot)
  fileout1=concat(c(fileout_prefix,"_1.pdf"))
  w=4
  pdf(file=fileout1, height=w*2, width=w*4)
  par(mfrow= c(1,1), mar = c(4,4,4,4))
  p1 <- DimPlot(object =merge1, reduction = "umap", group.by = "source", pt.size = .1, do.return = TRUE,label = T)
  p2 <- DimPlot(object =merge1, reduction = "umap", group.by = "cell_type", pt.size = .1, do.return = TRUE)
  p3 <- DimPlot(object =merge1, reduction = "umap", group.by = "scReference_label", pt.size = .1, do.return = TRUE)
  plot_grid(p1,p2,p3)
  dev.off()
}

Check_scReference_output(object, fileout_prefix)

Myeloid_plot_gene_expression<-function(object, fileout_prefix){
  file = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/MYELOID/cluster_annotation_gene_signatures.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  sign = p[,"annot"]
  gene.signature = strsplit(gsub(" ", "",p[,"gene.signature"],fixed = T),",")
  names(gene.signature) = sign
  gene.signature1 = gene.signature
  sign_genes = sort(unique(unlist(gene.signature)))
  genes = rownames(object@ assays$RNA@ counts)
  
  file = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/MYELOID/MYELOID_BLOOD/cluster_annotation_gene_signatures.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  sign = p[,"annot"]
  gene.signature = strsplit(gsub(" ", "",p[,"gene.list"],fixed = T),",")
  names(gene.signature) = sign
  gene.signature2 = gene.signature
  sign_genes1 = sort(unique(unlist(gene.signature)))
  sign_genes=sort(intersect(unique(sign_genes, sign_genes1), genes))
  
  genes =rownames(object@assays$ RNA@ data)	
  
  genes_of_interest =sort(unique(intersect(genes, sign_genes)))
  exp = object@ assays$ RNA@ data
  exp1 = as.matrix(exp[genes_of_interest,])
  cell_types = sort(unique(object@meta.data$refined_annotations))
  
  mat_mean = matrix(data=0, nrow = length(cell_types), ncol = length(genes_of_interest), dimnames = c(list(cell_types),list(genes_of_interest))) 
  for(i in c(1:length(cell_types))){
    w = rownames(object@meta.data)[which(object@meta.data$refined_annotations == cell_types[i])]
    mat_mean [cell_types[i],] = apply(exp1[,w], 1, mean)
  }
  mat_mean1= mat_mean[,which(apply(mat_mean, 2, sd)!=0)]
  mat_mean1  = mat_mean1 [,c(which(apply(mat_mean1, 2, max)>1))]
  
  hc = hclust(dist(t(mat_mean1)))
  labels = hc$ labels[hc$order] 
  hc = hclust(dist(mat_mean1))
  labels1 = hc$ labels[hc$order] 
  
  Idents(object = object) = object@meta.data$refined_annotations
  levels(x = object) = sort(unique(object@meta.data$refined_annotations))
  library(RColorBrewer)
  cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
  fileout1=concat(c(fileout_prefix,"_2.pdf"))
  w=4
  pdf(file=fileout1, height=w*1.9, width=w*8)
  par(mfrow= c(1,1), mar = c(5,5,3,3))
  DotPlot(object, assay = "RNA", features = labels) +
    scale_color_viridis_c() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
  
  xy = object@ reductions$ umap@ cell.embeddings[,c(1:2)]
  cell_type = object@meta.data$refined_annotations
  table(object@meta.data$refined_annotations, object@meta.data$orig.ident)
  
  orig.ident = object@meta.data$orig.ident
  w = grep("CD45p", orig.ident)
  orig.ident[w] = apply(cbind(object@meta.data$orig.ident[w],object@meta.data$Sample.Type[w]),1,paste,collapse = "_")
  orig.ident=gsub("_CD45p2_Tumour","_biopsy",gsub("_CD45p1_Tumour","_biopsy",orig.ident))
  orig.ident=gsub("_CD45p2_PBMC","_blood",gsub("_CD45p1_PBMC","_blood",orig.ident))
  sort(unique(orig.ident))
  
  out_file_table = concat(c(fileout_prefix,"_Overall_UMAP_", batch,".txt"))
  out = cbind(xy, apply(cbind("Myeloid" , cell_type), 1, paste, collapse = " "),orig.ident, object@meta.data$source, object@meta.data$Sample.Type)
  colnames(out) = c("UMAP1", "UMAP2","cell type","sample","source","sample_type")
  write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  t = table(object@meta.data$refined_annotations, object@meta.data$orig.ident)
  out_file_table = concat(c(fileout_prefix,"_Overall_cell_counts_", batch,".txt"))
  write.table(t(t), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}

B_cell_plot_gene_expression<-function(object, fileout_prefix){
  ########## plot full list of labels
  ref_cell_data = "/well/immune-rep/users/mfj169/10X_GENOMICS/REFERENCE_FILES/B_subset_marker_genes.txt"
  genes = rownames(pbmc@ assays$ RNA@ counts)
  p <- as.matrix(read.csv(ref_cell_data, head=T, sep="\t"))
  list_phenotypes = NULL
  phenotypes = NULL
  for(i in c(1:length(p[,1]))){
    g = strsplit(p[i,2], ",")[[1]]
    g = gsub(" ","", g)
    g = intersect(g, genes)
    if(length(g)>=3){
      phenotypes = c(phenotypes, p[i,1])
      list_phenotypes = c(list_phenotypes, list(g))
    }}
  names(list_phenotypes) = phenotypes
  
  genes =rownames(merge1@assays$ RNA@ data)	
  
  genes_of_interest = c("TXNIP FCER2 FCMR SELL BANK1 EGR1 CD69 CXCR4 DUSP2 JUN IL6 FCRL3 FCRL2 LINC01857 SAMSN1 SIGLEC10 RASSF6 FRZB HOPX BTNL9 FGFR1 IGHD IGHM IGHG1 CD27 CD38 SDC1 JCHAIN PRDM1 XBP1 MZB1 SSR4 GPR183 CD44 KLF2 TNFRSF13B VIM PLAC8 FCRL4 CCR1 ITGAX HMGB2 TUBA1B MKI67 AURKB CD19 MAPK1 CD81 CD40 CD22 CD24 CR2 CD3E CD7 HCST GZMA ZAP70 CD160")
  genes_of_interest = unique(c(unlist(strsplit(genes_of_interest," ", fixed = T)), unlist(list_phenotypes)))
  genes_of_interest = c(genes_of_interest, "IGHD","IGHM","IGHG1","IGHA1")
  genes_of_interest = sort(unique(intersect(genes_of_interest, genes)))
  genes_of_interest = setdiff(genes_of_interest,c( "CD3E"   ))

  genes =rownames(object@assays$ RNA@ data)	
  
  genes_of_interest =sort(unique(intersect(genes, c(genes_of_interest, sign_genes))))
  exp = object@ assays$ RNA@ data
  exp1 = as.matrix(exp[genes_of_interest,])
  cell_types = sort(unique(object@meta.data$refined_annotations))
  
  mat_mean = matrix(data=0, nrow = length(cell_types), ncol = length(genes_of_interest), dimnames = c(list(cell_types),list(genes_of_interest))) 
  for(i in c(1:length(cell_types))){
    w = rownames(object@meta.data)[which(object@meta.data$refined_annotations == cell_types[i])]
    mat_mean [cell_types[i],] = apply(exp1[,w], 1, mean)
  }
  mat_mean1= mat_mean[,which(apply(mat_mean, 2, sd)!=0)]
  mat_mean1  = mat_mean1 [,c(which(apply(mat_mean1, 2, max)>1))]
  
  hc = hclust(dist(t(mat_mean1)))
  labels = hc$ labels[hc$order] 
  hc = hclust(dist(mat_mean1))
  labels1 = hc$ labels[hc$order] 
  
  Idents(object = object) = object@meta.data$refined_annotations
  levels(x = object) = sort(unique(object@meta.data$refined_annotations))
  library(RColorBrewer)
  cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
  fileout1=concat(c(fileout_prefix,"_2.pdf"))
  w=4
  pdf(file=fileout1, height=w*1.9, width=w*8)
  par(mfrow= c(1,1), mar = c(5,5,3,3))
  DotPlot(object, assay = "RNA", features = labels) +
    scale_color_viridis_c() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
  
  xy = object@ reductions$ umap@ cell.embeddings[,c(1:2)]
  cell_type = object@meta.data$refined_annotations
  table(object@meta.data$refined_annotations, object@meta.data$orig.ident)
  
  orig.ident = object@meta.data$orig.ident
  w = grep("CD45p", orig.ident)
  orig.ident[w] = apply(cbind(object@meta.data$orig.ident[w],object@meta.data$Sample.Type[w]),1,paste,collapse = "_")
  orig.ident=gsub("_CD45p2_Tumour","_biopsy",gsub("_CD45p1_Tumour","_biopsy",orig.ident))
  orig.ident=gsub("_CD45p2_PBMC","_blood",gsub("_CD45p1_PBMC","_blood",orig.ident))
  sort(unique(orig.ident))
  
  out_file_table = concat(c(fileout_prefix,"_Overall_UMAP_", batch,".txt"))
  out = cbind(xy, apply(cbind("Myeloid" , cell_type), 1, paste, collapse = " "),orig.ident, object@meta.data$source, object@meta.data$Sample.Type)
  colnames(out) = c("UMAP1", "UMAP2","cell type","sample","source","sample_type")
  write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  t = table(object@meta.data$refined_annotations, object@meta.data$orig.ident)
  out_file_table = concat(c(fileout_prefix,"_Overall_cell_counts_", batch,".txt"))
  write.table(t(t), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}

T_cell_plot_gene_expression<-function(object, fileout_prefix){
  file = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/Ordered_gene_list_T_cells.txt"
  p <- as.matrix(read.csv(file, head=F, sep="\t"))
  T_cell_gene_set = p[,1]
  names(T_cell_gene_set) = p[,2]
  
  file = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/Ordered_gene_list_NK_cells.txt"
  p <- as.matrix(read.csv(file, head=F, sep="\t"))
  NK_cell_gene_set = p[,1]
  names(NK_cell_gene_set) = p[,2]
  
  genes =rownames(merge1@assays$ RNA@ data)
  labels = intersect(T_cell_gene_set ,genes)
  
  genes =rownames(object@assays$ RNA@ data)	
  
  genes_of_interest =sort(unique(intersect(genes, labels)))
  exp = object@ assays$ RNA@ data
  exp1 = as.matrix(exp[genes_of_interest,])
  cell_types = sort(unique(object@meta.data$refined_annotations))
  
  mat_mean = matrix(data=0, nrow = length(cell_types), ncol = length(genes_of_interest), dimnames = c(list(cell_types),list(genes_of_interest))) 
  for(i in c(1:length(cell_types))){
    w = rownames(object@meta.data)[which(object@meta.data$refined_annotations == cell_types[i])]
    mat_mean [cell_types[i],] = apply(exp1[,w], 1, mean)
  }
  mat_mean1= mat_mean[,which(apply(mat_mean, 2, sd)!=0)]
  mat_mean1  = mat_mean1 [,c(which(apply(mat_mean1, 2, max)>1))]
  
  hc = hclust(dist(t(mat_mean1)))
  labels = hc$ labels[hc$order] 
  hc = hclust(dist(mat_mean1))
  labels1 = hc$ labels[hc$order] 
  
  Idents(object = object) = object@meta.data$refined_annotations
  levels(x = object) = sort(unique(object@meta.data$refined_annotations))
  library(RColorBrewer)
  cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
  fileout1=concat(c(fileout_prefix,"_2.pdf"))
  w=4
  pdf(file=fileout1, height=w*1.9, width=w*8)
  par(mfrow= c(1,1), mar = c(5,5,3,3))
  DotPlot(object, assay = "RNA", features = labels) +
    scale_color_viridis_c() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
  
  xy = object@ reductions$ umap@ cell.embeddings[,c(1:2)]
  cell_type = object@meta.data$refined_annotations
  table(object@meta.data$refined_annotations, object@meta.data$orig.ident)
  
  orig.ident = object@meta.data$orig.ident
  w = grep("CD45p", orig.ident)
  orig.ident[w] = apply(cbind(object@meta.data$orig.ident[w],object@meta.data$Sample.Type[w]),1,paste,collapse = "_")
  orig.ident=gsub("_CD45p2_Tumour","_biopsy",gsub("_CD45p1_Tumour","_biopsy",orig.ident))
  orig.ident=gsub("_CD45p2_PBMC","_blood",gsub("_CD45p1_PBMC","_blood",orig.ident))
  sort(unique(orig.ident))
  
  out_file_table = concat(c(fileout_prefix,"_Overall_UMAP_", batch,".txt"))
  out = cbind(xy, apply(cbind("Myeloid" , cell_type), 1, paste, collapse = " "),orig.ident, object@meta.data$source, object@meta.data$Sample.Type)
  colnames(out) = c("UMAP1", "UMAP2","cell type","sample","source","sample_type")
  write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  t = table(object@meta.data$refined_annotations, object@meta.data$orig.ident)
  out_file_table = concat(c(fileout_prefix,"_Overall_cell_counts_", batch,".txt"))
  write.table(t(t), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}

if(cell_type_obejct == "myeloid"){Myeloid_plot_gene_expression(object, fileout_prefix)}
if(cell_type_obejct == "B_cell"){B_cell_plot_gene_expression(object, fileout_prefix)}
if(cell_type_obejct == "T_cell"){T_cell_plot_gene_expression(object, fileout_prefix)}

###################### Get pseudobulk output (to be completed soon!)
Get_pseudobulk_expression<-function(object){}
  
