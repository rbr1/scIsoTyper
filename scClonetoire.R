Draw_box_plot<-function(box,x,width,c,lwd,line_col){
	segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
	segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
	segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
	rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
	segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}
	
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha)) }

concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}

######################
batch = "PDAC150Ka"
analysis = "VDJ_clonality"

######################
##### get meta data (cell type) with VDJ information
input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"

##### get meta data (cell type) with VDJ information
Get_BCR_information<-function(input_directory_groups){
  type = "BCR"
  type1 = "B_cells"
  file = concat(c(input_directory_groups,"VDJ_information_",type,"_PDAC150Ka.txt"))
  VDJ <- as.matrix(read.csv(file, head=T, sep="\t"))
  clone = VDJ[,"clone1"]
  
  file = concat(c(input_directory_groups,"Cell_annotation_ALL_PDAC150Ka.txt"))
  ### this file can be replaced with the metadata information from your Seurat object
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  p1 <- p1[grep("B cell", p1[,"cell_type"]), ]
  id1 = p1[,"barcode"]
  sample = p1[,"sample"]
  id = apply(cbind(id1, sample), 1, function(x){
    x1 = gsub(concat(c(x[2], "_")), "", x[1])
    return(concat(c(x1,"||", x[2])))})
  id = gsub("-1","", id, fixed = T)
  length(intersect(id, names(clone)))
  cell_type = p1[,"cell_type"]
  Patient = p1[,"sample1"]
  Sample.Type= p1[,"sample"]
  Patient = gsub("_CD45p1","", Patient)
  Patient = gsub("_CD45p2","", Patient)
  Sample.Type[grep("biopsy", Sample.Type)] = "biopsy"
  Sample.Type[grep("blood", Sample.Type)] = "blood"
  names(cell_type) = id
  sample = apply(cbind(Patient, Sample.Type), 1, paste, collapse = "-")
  table(sample)
  samples = sort(unique(sample))
  cell_type1 = cell_type
  cell_type[grep("memory", cell_type)] = "B cell memory"
  cell_type[grep("activated", cell_type)] = "B cell activated"
  cell_type[grep("plasma", cell_type)] = "PB/PC"
  cell_type[grep("GC", cell_type)] = "B cell GC"
  
  names(cell_type) = id
  names(cell_type1) = id
  names(sample) = id
  
  list_return = c(list(VDJ), list(cell_type), list(cell_type1), list(sample), list(Patient), list(Sample.Type))
  names(list_return) = c("VDJ_object", "cell_type","cell_type_broad","sample","Patient","Sample.Type")
  return(list_return)
}

Get_TCR_information<-function(input_directory_groups){
  type = "TCR"
  type1 = "T_cells"
  file = concat(c(input_directory_groups,"VDJ_information_",type,"_PDAC150Ka.txt"))
  VDJ <- as.matrix(read.csv(file, head=T, sep="\t"))
  clone = VDJ[,"clone1"]
  
  file = concat(c(input_directory_groups,"Cell_annotation_ALL_PDAC150Ka.txt"))
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  p1 <- p1[grep("T cell", p1[,"cell_type"]), ]
  id1 = p1[,"barcode"]
  sample = p1[,"sample"]
  id = apply(cbind(id1, sample), 1, function(x){
    x1 = gsub(concat(c(x[2], "_")), "", x[1])
    return(concat(c(x1,"||", x[2])))})
  id = gsub("-1","", id, fixed = T)
  names(id) = id1
  length(intersect(id, names(clone)))
  cell_type = p1[,"cell_type"]
  Patient = p1[,"sample1"]
  Sample.Type= p1[,"sample"]
  Patient = gsub("_CD45p1","", Patient)
  Patient = gsub("_CD45p2","", Patient)
  Sample.Type[grep("biopsy", Sample.Type)] = "biopsy"
  Sample.Type[grep("blood", Sample.Type)] = "blood"
  sample = apply(cbind(Patient, Sample.Type), 1, paste, collapse = "-")
  
  file = concat(c(input_directory_groups,"Cell_names_broad.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,"Full.annotation"]%in% cell_type),]
  p=p[which(p[,"Very.broad.annotation"]%in% c("ILC","T cell gdT","T cell MAIT")==F),]
  Very.broad.annotation = p[,"Very.broad.annotation"]
  Broad.annotation = p[,"Broad.annotation"]
  names(Very.broad.annotation) = p[,"Full.annotation"]
  names(Broad.annotation) = p[,"Full.annotation"]
  
  cell_type1 = Broad.annotation[cell_type]
  
  names(cell_type) = id
  names(cell_type1) = id
  names(sample) = id
  names(Patient) = id
  
  list_return = c(list(VDJ), list(cell_type), list(cell_type1), list(sample), list(Patient), list(Sample.Type))
  names(list_return) = c("VDJ_object", "cell_type","cell_type_broad","sample","Patient","Sample.Type")
  return(list_return)
}

VDJ_list_BCR = Get_BCR_information(input_directory_groups)
VDJ_list_TCR = Get_TCR_information(input_directory_groups)


### BCR
BCR_subsampled_intra_subset_clonality<-function(VDJ_list_BCR, output_directory, batch, threshold_sample_depth = 5){
  VDJ_list = VDJ_list_BCR
  type = "BCR"
  type1 = "B_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  
  ### get total cell numbers for total and cell subsets
  count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
  totals = rowSums(count_per_sample)
  count_per_sample = cbind(count_per_sample, totals)
  
  print("Calculating duplicate rates")
  m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  table(sample[w_clone], cell_type_broad[w_clone])
  repeats = 1000
  thresholds1 = threshold_sample_depth
  for(c in c(1:length(cell_types_broad))){
    print (concat(c(c, " of ",length(cell_types_broad)," done" )))
    for(s in c(1:length(samples))){
      w2 = intersect(which(sample == samples[s]), which(cell_type_broad== cell_types_broad[c]))
      id_sub = intersect((names(sample)[w2]), w_clone)
      clon = clone[id_sub]
      clon = clon[which(clon!='-')]
      if(length(clon)>= thresholds1){
        n_dupl = NULL
        if(max(table(clon))==1){n_dupl=0
        }else{
          for(r in c(1:repeats)){
            rand = sample(clon, thresholds1)
            t = table(rand)
            n_dupl = c(n_dupl, length(which(t>1)))
          }}
        m_n_duplicate_clones[samples[s], cell_types_broad[c]] = mean(n_dupl)
      }}
  }
  all = c(cell_types )
  non_naive = setdiff(cell_types,"naive B cells")
  
  groups_cell_type = c(list(non_naive), list(all))
  groups_cell_type_names = c("non-naive","all")
  
  m_n_duplicate_clones_sub = matrix(data = -1,nrow = length(samples), ncol = length(groups_cell_type_names), dimnames = c(list(samples), list(groups_cell_type_names)))
  
  repeats = 1000
  thresholds1 = threshold_sample_depth
  for(c in c(1:length(groups_cell_type_names))){
    print (concat(c(c, " of ",length(groups_cell_type_names)," done" )))
    for(s in c(1:length(samples))){
      w2 = intersect(which(sample == samples[s]), which(cell_type %in% groups_cell_type[[c]]))
      id_sub = intersect((names(sample)[w2]), w_clone)
      clon = clone[id_sub]
      clon = clon[which(clon!='-')]
      if(length(clon)>= thresholds1){
        n_dupl = NULL
        if(max(table(clon))==1){n_dupl=0
        }else{
          for(r in c(1:repeats)){
            rand = sample(clon, thresholds1)
            t = table(rand)
            n_dupl = c(n_dupl, length(which(t>1)))
          }}
        m_n_duplicate_clones_sub[samples[s], groups_cell_type_names[c]] = mean(n_dupl)
      }}
  }
  m_n_duplicate_clones1 = cbind(m_n_duplicate_clones, m_n_duplicate_clones_sub)
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_intra.txt"))
  write.table(m_n_duplicate_clones1, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}
BCR_subsampled_intra_subset_clonality(VDJ_list_BCR, output_directory,batch,threshold_sample_depth = 5)

BCR_subsampled_inter_subset_clonality<-function(VDJ_list_BCR, output_directory, batch){
  VDJ_list = VDJ_list_BCR
  type = "BCR"
  type1 = "B_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  
  ### get total cell numbers for total and cell subsets
  count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
  totals = rowSums(count_per_sample)
  count_per_sample = cbind(count_per_sample, totals)
  
  print("Calculating inter duplicate rates")
  m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  repeats = 1000
  thresholds1 = min(c(50,quantile(totals, 0.05)))
  for(s in c(1:length(samples))){
    print (concat(c(s, " of ",length(samples)," done" )))
    w1 = intersect(names(sample)[which(sample == samples[s])], w_clone)
    repeated_measures = matrix(data = 0,nrow = length(c(1:repeats)), ncol = length(cell_types_broad), dimnames = c(list(c(1:repeats)), list(cell_types_broad)))
    if(length(w1)>= thresholds1){
      for(r in c(1:repeats)){
        rand = sample(w1, thresholds1)
        t = table(clone[rand], cell_type_broad[rand])
        expanded = which(rowSums(t)>1)
        if(length(expanded)==1){
          repeated_measures[r,] = repeated_measures[r,]+t[expanded,]
        }else{
          if(length(expanded)>1){
            repeated_measures[r,] = repeated_measures[r,]+colSums(t[expanded,])
          }}}
      if(sum(repeated_measures)>=5){
        m_n_duplicate_clones[samples[s],]= colSums(repeated_measures)*100/sum(repeated_measures)
      }}
    }
      
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_inter.txt"))
  write.table(m_n_duplicate_clones, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}
BCR_subsampled_inter_subset_clonality(VDJ_list_BCR, output_directory,batch)

BCR_renyi_gini_d5_subsampled<-function(VDJ_list_BCR, output_directory){
  type = "BCR"
  type1 = "B_cells"
  VDJ_list = VDJ_list_BCR
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  
  count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
  totals = rowSums(count_per_sample)
  count_per_sample = cbind(count_per_sample, totals)
  
  #### clonal expansion from subsampled cells. 
  m_gini = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  m_renyi = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  m_d5 = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  
  print("Calculating Gini, Renyi and D5")
  thresholds = 15
  repeats = 1000
  library(ineq)
  library(edgeR)
  for(c in c(1:length(cell_types_broad))){
    print (concat(c(c, " of ",length(cell_types_broad)," done" )))
    g1_gini = NULL
    for(s in c(1:length(samples))){
      w2 = intersect(which(sample == samples[s]), which(cell_type_broad== cell_types_broad[c]))
      id_sub = intersect((names(sample)[w2]), w_clone)
      clon = clone[id_sub]
      clon = clon[which(clon!='-')]
      if(length(clon)>= thresholds){
        ginis = NULL
        reynis = NULL
        d5s = NULL
        for(r in c(1:repeats)){
          rand = sample(clon, thresholds)
          dist = table(rand)
          g = ineq::Gini(dist)
          ginis = c(ginis, g)
          dist = dist/sum(dist)
          r = sum(-(dist*log(dist)))
          renyis = c(reynis, r)
          d5 = rev(sort(table(rand)))
          if(length(d5)>=3){d5 = sum(d5[c(1:3)])*100/sum(d5)
          }else{d5 = 100}
          d5s = c(d5s, d5)
        }
        g1_gini = c(g1_gini,list(ginis))
        m_gini[samples[s], cell_types_broad[c]] = mean(ginis)
        m_renyi[samples[s], cell_types_broad[c]] = mean(renyis)
        m_d5[samples[s], cell_types_broad[c]] = mean(d5s)
      }}
  }
  
  gini_all = rep(-1,length(samples))
  names(gini_all)= samples
  renyi_all = gini_all
  d5_all = gini_all
  thresh = floor(min(totals)*0.95)
  for(s in c(1:length(samples))){
    w2 =which(sample == samples[s])
    id_sub = names(sample)[w2]
    clon = clone[id_sub]
    clon = clon[which(clon!='-')]
    if(length(clon)>= thresh){
      ginis = NULL
      reynis = NULL
      for(r in c(1:repeats)){
        rand = sample(clon, thresh)
        dist = table(rand)
        g = ineq::Gini(dist)
        ginis = c(ginis, g)
        dist = dist/sum(dist)
        r = sum(-(dist*log(dist)))
        renyis = c(reynis, r)
        d5 = rev(sort(table(rand)))
        if(length(d5)>5){d5 = sum(d5[c(1:5)])*100/sum(d5)
        }else{d5 = 100}
        d5s = c(d5s, d5)
      }
      gini_all[samples[s]] = mean(ginis)
      renyi_all[samples[s]] = mean(renyis)
      d5_all[samples[s]] = mean(d5s)
    }}	
  
  m_gini = cbind(m_gini, gini_all)
  m_renyi = cbind(m_renyi, renyi_all)
  m_d5 = cbind(m_d5, d5_all)
  colnames(m_d5)[length(colnames(m_d5))] = "all"
  colnames(m_gini)[length(colnames(m_gini))] = "all"
  colnames(m_renyi)[length(colnames(m_renyi))] = "all"
  
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_gini.txt"))
  write.table(m_gini, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_renyi.txt"))
  write.table(m_renyi, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_d5.txt"))
  write.table(m_d5, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}
BCR_renyi_gini_d5_subsampled(VDJ_list_BCR, output_directory)

BCR_isotype_SHM_VJ<-function(VDJ_list_BCR, output_directory){
  VDJ_list = VDJ_list_BCR
  type = "BCR"
  type1 = "B_cells"
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  isotype = VDJ[,"constant_region1"]
  IGHKL = VDJ[,"chain2"]
  SHM1 = VDJ[,"V_mm1"]
  SHM2 = VDJ[,"V_mm2"]
  V_gene1 = VDJ[,"V_gene1"]
  V_gene2 = VDJ[,"V_gene2"]
  VJ_gene1 = apply(cbind(VDJ[,"V_gene1"], VDJ[,"J_gene1"]), 1, paste, collapse = ":")
  VJ_gene2 = apply(cbind(VDJ[,"V_gene2"], VDJ[,"J_gene2"]), 1, paste, collapse = ":")
  SHM1[which(SHM1=="-")] = 0
  SHM2[which(SHM2=="-")] = 0
  SHM1 = as.numeric(SHM1)
  SHM2 = as.numeric(SHM2)
  names(SHM1) = rownames(VDJ)
  names(SHM2) = rownames(VDJ)
  
  isotypes = sort(unique(isotype))
  IGHKLs = sort(unique(IGHKL))
  V_gene1s = sort(unique(V_gene1))
  V_gene2s = sort(unique(V_gene2))
  VJ_gene1s = sort(unique(VJ_gene1))
  VJ_gene2s = sort(unique(VJ_gene2))
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  
  group_ids = NULL
  for(c in c(1:length(cell_types_broad))){
    w = intersect(w_clone, names(cell_type_broad)[which(cell_type_broad == cell_types_broad[c])])
    group_ids = c(group_ids, list(w))
  }
  for(c in c(1:length(cell_types))){
    w = intersect(w_clone, names(cell_type)[which(cell_type == cell_types[c])])
    group_ids = c(group_ids, list(w))
  }
  group_ids = c(group_ids, list(w_clone))
  names(group_ids) = c(apply(cbind(cell_types_broad,"broad"),1, paste, collapse = ":"), 
                       apply(cbind(cell_types,"detailed"),1, paste, collapse = ":") , "All" )
  
  list_isotype_percentage = NULL
  list_IGHKL_percentage = NULL
  list_meanSHM1 = NULL ## by isotype
  list_meanSHM2 = NULL ## by isotype
  list_meanSHM12 = NULL ## by isotype
  list_IGHV_gene_percentage = NULL
  list_IGKLV_gene_percentage = NULL
  list_IGHVJ_gene_percentage = NULL
  list_IGKLVJ_gene_percentage = NULL
  
  for(c in c(1:length(group_ids))){
    m_isotype_percentage = matrix(data = 0,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_IGHKL_percentage = matrix(data = 0,nrow = length(samples), ncol = length(IGHKLs), dimnames = c(list(samples), list(IGHKLs)))
    m_meanSHM1 = matrix(data = NA,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_meanSHM2 = matrix(data = NA,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_meanSHM12 = matrix(data = NA,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_IGHV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene1s), dimnames = c(list(samples), list(V_gene1s)))
    m_IGKLV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene2s), dimnames = c(list(samples), list(V_gene2s)))
    m_IGHVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene1s), dimnames = c(list(samples), list(VJ_gene1s)))
    m_IGKLVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene2s), dimnames = c(list(samples), list(VJ_gene2s)))
    
    w  = group_ids[[c]]
    t = table(sample[w],isotype[w])
    m_isotype_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],IGHKL[w])
    m_IGHKL_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],V_gene1[w])
    m_IGHV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],V_gene2[w])
    m_IGKLV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene1[w])
    m_IGHVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene2[w])
    m_IGKLVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    df = data.frame(w)
    df$SHM1 =SHM1[w]
    df$SHM2 =SHM2[w]
    df$SHM12 =SHM1[w]+SHM2[w]
    df$sample =sample[w]
    df$isotype =isotype[w]
    
    dd<-with(df, aggregate(SHM1, by=list(sample=sample, isotype=isotype), mean))
    dn<-with(df, aggregate(SHM1, by=list(sample=sample, isotype=isotype), length))
    dd = dd[which(dn$x>=5),]
    ddd<-reshape(dd, idvar='sample', timevar='isotype', direction='wide')
    colnames(ddd) = gsub("x.", "", colnames(ddd))
    rownames(ddd) = ddd$sample
    ddd = ddd[,intersect(colnames(ddd), isotypes)]
    m_meanSHM1[rownames(ddd), colnames(ddd)] = as.matrix(ddd)
    
    dd<-with(df, aggregate(SHM2, by=list(sample=sample, isotype=isotype), mean))
    dn<-with(df, aggregate(SHM2, by=list(sample=sample, isotype=isotype), length))
    dd = dd[which(dn$x>=5),]
    ddd<-reshape(dd, idvar='sample', timevar='isotype', direction='wide')
    colnames(ddd) = gsub("x.", "", colnames(ddd))
    rownames(ddd) = ddd$sample
    ddd = ddd[,intersect(colnames(ddd), isotypes)]
    m_meanSHM2[rownames(ddd), colnames(ddd)] = as.matrix(ddd)
    
    dd<-with(df, aggregate(SHM12, by=list(sample=sample, isotype=isotype), mean))
    dn<-with(df, aggregate(SHM12, by=list(sample=sample, isotype=isotype), length))
    dd = dd[which(dn$x>=5),]
    ddd<-reshape(dd, idvar='sample', timevar='isotype', direction='wide')
    colnames(ddd) = gsub("x.", "", colnames(ddd))
    rownames(ddd) = ddd$sample
    ddd = ddd[,intersect(colnames(ddd), isotypes)]
    m_meanSHM12[rownames(ddd), colnames(ddd)] = as.matrix(ddd)
    
    list_isotype_percentage = c(list_isotype_percentage, list( m_isotype_percentage))
    list_IGHKL_percentage = c(list_IGHKL_percentage, list( m_IGHKL_percentage))
    list_meanSHM1 = c(list_meanSHM1, list( m_meanSHM1))
    list_meanSHM2 = c(list_meanSHM2, list( m_meanSHM2))
    list_meanSHM12 = c(list_meanSHM12, list( m_meanSHM12))
    list_IGHV_gene_percentage = c(list_IGHV_gene_percentage, list( m_IGHV_gene_percentage))
    list_IGKLV_gene_percentage = c(list_IGKLV_gene_percentage, list( m_IGKLV_gene_percentage))
    list_IGHVJ_gene_percentage = c(list_IGHVJ_gene_percentage, list( m_IGHVJ_gene_percentage))
    list_IGKLVJ_gene_percentage = c(list_IGKLVJ_gene_percentage, list( m_IGKLVJ_gene_percentage))
  }
  
  names(list_isotype_percentage) = names(group_ids) 
  names(list_IGHKL_percentage) = names(group_ids) 
  names(list_meanSHM1) = names(group_ids) 
  names(list_meanSHM2) = names(group_ids) 
  names(list_meanSHM12) = names(group_ids) 
  names(list_IGHV_gene_percentage) = names(group_ids) 
  names(list_IGKLV_gene_percentage) = names(group_ids) 
  names(list_IGHVJ_gene_percentage ) = names(group_ids) 
  names(list_IGKLVJ_gene_percentage) = names(group_ids) 
  
  list_VDJ_freatures_aggregated  = c(list(list_isotype_percentage),
                                     list(list_IGHKL_percentage),
                                     list(list_meanSHM1),
                                     list(list_meanSHM2),
                                     list(list_meanSHM12),
                                     list(list_IGHV_gene_percentage),
                                     list(list_IGKLV_gene_percentage),
                                     list(list_IGHVJ_gene_percentage),
                                     list(list_IGKLVJ_gene_percentage)
  )
  names(list_VDJ_freatures_aggregated)  = c("isotype usage","IGK/L usage","mean SHM IGH","mean SHM IGK/L","mean SHM IGH+K/L",
                                            "IGHV gene usage", "IGK/L V gene usage","IGHVJ gene usage","IGK/K VJ gene usage")
  
  saveRDS(file = concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")), list_VDJ_freatures_aggregated)
  writeLines(concat(c("\nOutput as RDS object with the following features calculated per cell type:\n",
                      paste(c(" ", names(list_VDJ_freatures_aggregated)), collapse = "\n\t"),
                      "\n\nOutput location:", concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")))))
}
BCR_isotype_SHM_VJ(VDJ_list_BCR, output_directory)

Plot_clone_sizes_per_cell_type_BCR<-function(VDJ_list_BCR, output_directory,input_directory_groups){
  VDJ_list = VDJ_list_BCR
  type = "BCR"
  type1 = "B_cells"
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  source = strsplit(sample,"-",fixed = T)
  for(i in c(1:length(source))){source[i]=source[[i]][2]}
  source = unlist(source)
  table(source)

  sources = sort(unique(source))
  
  for(s in c(1:length(sources))){
    w_source = names(source)[which(source==sources[s])]
    w = intersect(w_source, w_clone)
    t = table(cell_type_broad[w], clone[w])
    list_clone_sizes = NULL
    cell_type_group = NULL
    for(i in c(1:length(t[,1]))){
      x = t[i,]
      x = x[which(x!=0)]
      list_clone_sizes = c(list_clone_sizes, x)
      cell_type_group = c(cell_type_group, rep(rownames(t)[i], length(x)))
    }
    df = data.frame(c(1:length(list_clone_sizes)))
    df$clone_sizes = list_clone_sizes
    df$cell_type = cell_type_group
    
    library(ggplot2)
    fileout1=concat(c(output_directory, "Clone_size_distribution_BCR_", sources[s],".pdf"))
    w=6
    pdf(file=fileout1, height=w*0.85, width=w*1.1)
    par(mfrow= c(1,1), mar = c(1,1,1,1))
    if(s ==1){
      ggplot(df, aes(x = clone_sizes)) +
        geom_histogram(fill = "darkgreen", colour = "darkgreen", binwidth = 1, pad=T) +
        facet_grid(cell_type ~ .) +
        scale_y_continuous(trans='log10')}
    if(s ==2){
      ggplot(df, aes(x = clone_sizes)) +
        geom_histogram(fill = "darkgreen", colour = "darkgreen", pad=T) +
        facet_grid(cell_type ~ .) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') 
    }
    dev.off()
    print(summary(df$clone_sizes[which(cell_type=="B cell naive")]))
    x = t["B cell naive",]
    x = x[which(x!=0)]
    
  }
  
  
}
Plot_clone_sizes_per_cell_type_BCR(VDJ_list_BCR, output_directory,input_directory_groups)

### TCR

## identify CD4 and CD8 T cell subtypes: the user may need to edit this depending on annotations
VDJ_list = VDJ_list_TCR
cell_type_broad = VDJ_list[[2]]
cell_type = VDJ_list[[3]]
cell_type_broads = sort(unique(cell_type_broad))
cell_types = sort(unique(cell_type))
CD4 = sort(unique(c(cell_type_broads[grep("CD4",cell_type_broads)], cell_types[grep("CD4",cell_types)])))
CD8 = sort(unique(c(cell_type_broads[grep("CD8",cell_type_broads)], cell_types[grep("CD8",cell_types)], 
                    cell_type_broads[grep("MAIT",cell_type_broads)], cell_types[grep("MAIT",cell_types)])))
CD48 = c(list(CD4),list(CD8))
names(CD48)= c("CD4","CD8")

TCR_subsampled_inter_subset_clonality<-function(VDJ_list_TCR, output_directory, batch, CD48, threshold_sample_depth = 5){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    wx = names(sample)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[wx])
    w_clone = intersect(names(which(clone!='-')),wx)
    
    ### get total cell numbers for total and cell subsets
    count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
    totals = rowSums(count_per_sample)
    count_per_sample = cbind(count_per_sample, totals)
    
    print("Calculating duplicate rates")
    m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
    table(sample[w_clone], cell_type_broad[w_clone])
    repeats = 1000
    thresholds1 = threshold_sample_depth
    for(c in c(1:length(cell_types_broad))){
      print (concat(c(c, " of ",length(cell_types_broad)," done" )))
      for(s in c(1:length(samples))){
        w2 = intersect(which(sample == samples[s]), which(cell_type_broad== cell_types_broad[c]))
        id_sub = intersect((names(sample)[w2]), w_clone)
        clon = clone[id_sub]
        clon = clon[which(clon!='-')]
        if(length(clon)>= thresholds1){
          n_dupl = NULL
          if(max(table(clon))==1){n_dupl=0
          }else{
            for(r in c(1:repeats)){
              rand = sample(clon, thresholds1)
              t = table(rand)
              n_dupl = c(n_dupl, length(which(t>1)))
            }}
          m_n_duplicate_clones[samples[s], cell_types_broad[c]] = mean(n_dupl)
        }}
    }
    all = c(cell_types )
    non_naive = setdiff(cell_types,"naive B cells")
    
    groups_cell_type = c(list(non_naive), list(all))
    groups_cell_type_names = c("non-naive","all")
    
    m_n_duplicate_clones_sub = matrix(data = -1,nrow = length(samples), ncol = length(groups_cell_type_names), dimnames = c(list(samples), list(groups_cell_type_names)))
    
    repeats = 1000
    thresholds1 = threshold_sample_depth
    for(c in c(1:length(groups_cell_type_names))){
      print (concat(c(c, " of ",length(groups_cell_type_names)," done" )))
      for(s in c(1:length(samples))){
        w2 = intersect(which(sample == samples[s]), which(cell_type %in% groups_cell_type[[c]]))
        id_sub = intersect((names(sample)[w2]), w_clone)
        clon = clone[id_sub]
        clon = clon[which(clon!='-')]
        if(length(clon)>= thresholds1){
          n_dupl = NULL
          if(max(table(clon))==1){n_dupl=0
          }else{
            for(r in c(1:repeats)){
              rand = sample(clon, thresholds1)
              t = table(rand)
              n_dupl = c(n_dupl, length(which(t>1)))
            }}
          m_n_duplicate_clones_sub[samples[s], groups_cell_type_names[c]] = mean(n_dupl)
        }}
    }
    m_n_duplicate_clones1 = cbind(m_n_duplicate_clones, m_n_duplicate_clones_sub)
    
    out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_intra.txt"))
    write.table(m_n_duplicate_clones1, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  }
}
TCR_subsampled_inter_subset_clonality(VDJ_list_TCR, output_directory,batch,CD48, threshold_sample_depth = 5)

TCR_duplicate_counts_subsampled_inter<-function(VDJ_list_TCR, output_directory, batch,CD48){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    wx = names(sample)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[wx])
    w_clone = intersect(names(which(clone!='-')),wx)
    
    table(cell_type[w_clone])
    table(cell_type[wx])
    
    ### get total cell numbers for total and cell subsets
    count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
    totals = rowSums(count_per_sample)
    count_per_sample = cbind(count_per_sample, totals)
    
    print("Calculating inter duplicate rates")
    m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
    repeats = 1000
    thresholds1 = min(c(50,quantile(totals, 0.05)))
    for(s in c(1:length(samples))){
      print (concat(c(s, " of ",length(samples)," done" )))
      w1 = intersect(names(sample)[which(sample == samples[s])], w_clone)
      repeated_measures = matrix(data = 0,nrow = length(c(1:repeats)), ncol = length(cell_types_broad), dimnames = c(list(c(1:repeats)), list(cell_types_broad)))
      if(length(w1)>= thresholds1){
        for(r in c(1:repeats)){
          rand = sample(w1, thresholds1)
          t = table(clone[rand], cell_type_broad[rand])
          expanded = which(rowSums(t)>1)
          if(length(expanded)==1){
            repeated_measures[r,colnames(t)] = repeated_measures[r,colnames(t)]+t[expanded,]
          }else{
            if(length(expanded)>1){
              repeated_measures[r,colnames(t)] = repeated_measures[r,colnames(t)]+colSums(t[expanded,])
            }}}
        if(sum(repeated_measures)>=5){
          m_n_duplicate_clones[samples[s],]= colSums(repeated_measures)*100/sum(repeated_measures)
        }}
    }
    
    out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_inter.txt"))
    write.table(m_n_duplicate_clones, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  }
}
TCR_duplicate_counts_subsampled_inter(VDJ_list_TCR, output_directory,batch,CD48)

### TCR v gene usages
TCR_VJ<-function(VDJ_list_TCR, output_directory,CD48){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[3]]
  cell_type = VDJ_list[[2]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  V_gene1 = VDJ[,"V_gene1"]
  V_gene2 = VDJ[,"V_gene2"]
  VJ_gene1 = apply(cbind(VDJ[,"V_gene1"], VDJ[,"J_gene1"]), 1, paste, collapse = ":")
  VJ_gene2 = apply(cbind(VDJ[,"V_gene2"], VDJ[,"J_gene2"]), 1, paste, collapse = ":")
  
  V_gene1s = sort(unique(V_gene1))
  V_gene2s = sort(unique(V_gene2))
  VJ_gene1s = sort(unique(VJ_gene1))
  VJ_gene2s = sort(unique(VJ_gene2))
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  
  group_ids = NULL
  for(c in c(1:length(cell_types_broad))){
    w = intersect(w_clone, names(cell_type_broad)[which(cell_type_broad == cell_types_broad[c])])
    group_ids = c(group_ids, list(w))
  }
  for(c in c(1:length(cell_types))){
    w = intersect(w_clone, names(cell_type)[which(cell_type == cell_types[c])])
    group_ids = c(group_ids, list(w))
  }
  group_ids = c(group_ids, list(w_clone),CD48)
  
  names(group_ids) = c(apply(cbind(cell_types_broad,"broad"),1, paste, collapse = ":"), 
                       apply(cbind(cell_types,"detailed"),1, paste, collapse = ":") ,names(CD48) )
  
  list_IGHV_gene_percentage = NULL
  list_IGKLV_gene_percentage = NULL
  list_IGHVJ_gene_percentage = NULL
  list_IGKLVJ_gene_percentage = NULL
  for(c in c(1:length(group_ids))){
    m_IGHV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene1s), dimnames = c(list(samples), list(V_gene1s)))
    m_IGKLV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene2s), dimnames = c(list(samples), list(V_gene2s)))
    m_IGHVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene1s), dimnames = c(list(samples), list(VJ_gene1s)))
    m_IGKLVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene2s), dimnames = c(list(samples), list(VJ_gene2s)))
    
    t = table(sample[w],V_gene1[w])
    m_IGHV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],V_gene2[w])
    m_IGKLV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene1[w])
    m_IGHVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene2[w])
    m_IGKLVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    list_IGHV_gene_percentage = c(list_IGHV_gene_percentage, list( m_IGHV_gene_percentage))
    list_IGKLV_gene_percentage = c(list_IGKLV_gene_percentage, list( m_IGKLV_gene_percentage))
    list_IGHVJ_gene_percentage = c(list_IGHVJ_gene_percentage, list( m_IGHVJ_gene_percentage))
    list_IGKLVJ_gene_percentage = c(list_IGKLVJ_gene_percentage, list( m_IGKLVJ_gene_percentage))
  }
  
  names(list_IGHV_gene_percentage) = names(group_ids) 
  names(list_IGKLV_gene_percentage) = names(group_ids) 
  names(list_IGHVJ_gene_percentage ) = names(group_ids) 
  names(list_IGKLVJ_gene_percentage) = names(group_ids) 
  
  list_VDJ_freatures_aggregated  = c(list(list_IGHV_gene_percentage),
                                     list(list_IGKLV_gene_percentage),
                                     list(list_IGHVJ_gene_percentage),
                                     list(list_IGKLVJ_gene_percentage)
  )
  names(list_VDJ_freatures_aggregated)  = c("TRAV gene usage", "TRBV gene usage","TRAVJ gene usage","TRBVJ gene usage")
  
  saveRDS(file = concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")), list_VDJ_freatures_aggregated)
  writeLines(concat(c("\nOutput as RDS object with the following features calculated per cell type:\n",
                      paste(c(" ", names(list_VDJ_freatures_aggregated)), collapse = "\n\t"),
                      "\n\nOutput location:", concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")))))
  
  
  }
TCR_VJ(VDJ_list_TCR, output_directory)

Plot_clone_sizes_per_cell_type_TCR<-function(VDJ_list_TCR, output_directory,input_directory_groups){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  VDJ = VDJ_list[[1]]
  cell_type  = VDJ_list[[2]]
  cell_type_broad = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  source = strsplit(sample,"-",fixed = T)
  for(i in c(1:length(source))){source[i]=source[[i]][2]}
  source = unlist(source)
  table(source)
  cell_types_broad
  
  sources = sort(unique(source))
  
  for(s in c(1:length(sources))){
    w_source = names(source)[which(source==sources[s])]
    w = intersect(w_source, w_clone)
    t = table(cell_type_broad[w], clone[w])
    list_clone_sizes = NULL
    cell_type_group = NULL
    for(i in c(1:length(t[,1]))){
      x = t[i,]
      x = x[which(x!=0)]
      list_clone_sizes = c(list_clone_sizes, x)
      cell_type_group = c(cell_type_group, rep(rownames(t)[i], length(x)))
    }
    df = data.frame(c(1:length(list_clone_sizes)))
    df$clone_sizes = list_clone_sizes
    df$cell_type = gsub("T cell","",cell_type_group)
    
    library(ggplot2)
    fileout1=concat(c(output_directory, "Clone_size_distribution_TCR_", sources[s],".pdf"))
    w=6
    pdf(file=fileout1, height=w*1.7, width=w*1.1)
    par(mfrow= c(1,1), mar = c(1,1,1,1))
    if(s ==1){
      ggplot(df, aes(x = clone_sizes)) +
        geom_histogram(fill = "darkgreen", colour = "darkgreen", binwidth = 1, pad=T) +
        facet_grid(cell_type ~ .) +
        scale_y_continuous(trans='log10')}
    if(s ==2){
      ggplot(df, aes(x = clone_sizes)) +
        geom_histogram(fill = "darkgreen", colour = "darkgreen", pad=T) +
        facet_grid(cell_type ~ .) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') 
    }
    dev.off()
    
  }
  
  
}
Plot_clone_sizes_per_cell_type_TCR(VDJ_list_TCR, output_directory,input_directory_groups)

print("You can now plot between patient groups using these outputs")
  



