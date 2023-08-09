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
groups_PCA = readRDS(file = concat(c(input_directory_groups, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))

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


################# clonal overlap between subsets

### BCR overlap within samples between subsets
BCR_clonal_overlap_per_site<-function(VDJ_list_BCR, output_directory,groups_PCA){
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
  
  m_cell_types = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
  m_cell_types_proportions = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
  for(s in c(1:length(samples))){
    w = intersect(w_clone , names(which(sample == samples[s])))
    t = table(cell_type[w])
    m_cell_types[samples[s], names(t)]= t
    m_cell_types_proportions[samples[s], names(t)]= t*100/sum(t)
  }
  totals = rowSums(m_cell_types)
  
  ############# get overlap between compartments
  c1 = NULL
  c2 = NULL
  for(i1 in c(1:length(cell_types))){
    for(i2 in c(i1:length(cell_types))){
      if(i1<i2){
        c1 = c(c1, cell_types[i1])
        c2 = c(c2, cell_types[i2])}}}
  sharing_names = apply(cbind(c1,c2), 1, paste,collapse = " - ")
  
  print("Running without subsampling")
  m_overlap_raw = matrix(data = 0,nrow = length(samples), ncol = length(sharing_names), dimnames = c(list(samples), list(sharing_names)))
  
  for(s in c(1:length(samples))){
    print (s)
    w_sample = intersect(w_clone , names(which(sample == samples[s])))
    for(i in c(1:length(c1))){
      w1 = intersect(w_sample, names(which(cell_type==c1[i])))
      w2 = intersect(w_sample, names(which(cell_type==c2[i])))
      shared_clones = intersect(clone[w1],clone[w2])
      m_overlap_raw[samples[s], sharing_names[i]] = length(shared_clones)
    }}
  
  # normalise 
  list_overlapping_per_sample = NULL
  for(s in c(1:length(samples))){
    m_overlapping = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
    m_overlapping_total = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
    wa = intersect(names(sample)[which(sample== samples[s])], w_clone)
    t = table(clone[wa], cell_type[wa])
    w = which(apply(t, 1, function(x){length(which(x!=0))})>1)
    if(length(w)>0){
      # print(r)
      t=rbind(t[w,])
      for(i in c(1:length(w))){
        x = names(t[i,][which(t[i,]!=0)])
        for(i1 in c(1:(length(x)-1))){
          for(i2 in c((i1+1):length(x))){
            m_overlapping_total[x[i1],x[i2]] = m_overlapping_total[x[i1],x[i2]]+1
          }}
      }
    }
    
    scale = NULL
    if(sum(m_overlapping_total)==0){scale = 1
    m_overlapping_total = m_overlapping_total#/sum(m_overlapping_total)
    }else{m_overlapping_total = m_overlapping_total/sum(m_overlapping_total)}
   
    scale = 1
    list_overlapping_per_sample = c(list_overlapping_per_sample, list(m_overlapping_total*mean(scale)))
    print (s)
  }
  names(list_overlapping_per_sample) = samples
    
  groups = NULL
  links = sharing_names
  
  overlapping_per_link = matrix(data = 0,nrow = length(samples), ncol = length(links), dimnames = c(list(samples), list(links)))
  for(s in c(1:length(samples))){
    mat = list_overlapping_per_sample[[s]]
    for(c1 in c(1:length(cell_types))){
      for(c2 in c(c1:length(cell_types))){
        if(c1<c2){
          if(is.na(mat[cell_types[c1],cell_types[c2]])==F){
            link = concat(c(cell_types[c1], " - ", cell_types[c2]))
            overlapping_per_link[samples[s],link] = mat[cell_types[c1],cell_types[c2]]
          }}}}}
  
  m_overlap_sampled1 = overlapping_per_link
  saveRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")), m_overlap_sampled1)
  
  
  #################### plot by PCA group
  m_overlap_sampled1  = readRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")))
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  mat = m_overlap_sampled1
  nz = apply(mat,2,function(x){length(which(x!=0))})
  mat = mat[,which(nz>=1)]
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  w = which(group_PCA_list %in% rownames(mat))
  mat_stat = mat[group_PCA_list[w],]
  
  factor = factor(factor_PCA[w])
  mat_stat1 = mat_stat^0.25
  fit = manova(formula = mat_stat1 ~ factor)
  p1 = summary.aov(fit)
  nam = gsub(" Response ","",names(p1))
  p_value = NULL
  means = NULL
  
  i1 = 0
  for(i in p1){
    i1 = i1+1
    p_value = c(p_value, i$'Pr(>F)'[1]) 
    if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
    }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
  }
  p_value[which(is.na(p_value))] = 2
  print(min(p_value))
  length(which(p_value<0.05))
  
  colnames(means) = paste("mean.group.", c(1:length(means[1,])))
  combined_p_value = cbind(p_value ,means)
  name = colnames(mat_stat)
  rownames(combined_p_value) = nam
  p.site = rep(concat(c("biopsy")), length(nam))
  analysis = "overlap"
  p.analysis = rep(analysis, length(nam))
  x = cbind(p.site, p.analysis, combined_p_value)
  summary_stats = x
  
  
  groups = NULL
  for(i in c(1:length(mat_stat[1,]))){
    g1 = NULL
    for(s in c(1:length(groups_PCA))){
      g1 = c(g1, list(mat_stat[ gsub("_","-",groups_PCA[[s]]),i]))}
    groups = c(groups, list(g1))
  }
  names(groups) = colnames(mat_stat)
  
  ######################################
  analysis = "Clonal_overlap_B_cell_subsets"
  
  outfile = concat(c(output_directory,"Stats_",analysis,"_biopsy_", batch,".txt"))
  
  write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  summary_stats_biopsy = summary_stats
  
  fileout1=concat(c(output_directory,"",analysis,"_biopsy_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.8, width=w*3.7)
  par(mfrow= c(1,1), mar = c(20,4,2.5,1))
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  p_values = p_value
  factors1 = paste("group",c(1:length(groups_PCA)))
  factors = names(groups)
  main = concat(c("clonal overlap:biopsy"))
  max = max(c(unlist(groups), unlist(groups))*1.3)
  min = 0
  b = (max-min)*0.035
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = min(c(max,100))
  range = max-min
  if(range>50){scale = c(0:100)*20}
  if(range<=50){scale = c(0:100)*10}
  if(range<=30){scale = c(0:100)*5}
  if(range <15){scale = c(0:100)*2.5}
  if(range <5){scale = c(0:100)*1}
  if(range <4){scale = c(0:100)*0.5}
  if(range <1.5){scale = c(0:1000)*0.2}
  if(range <0.5){scale = c(0:100)*0.1}
  if(range <0.1){scale = c(0:100)*0.01}
  if(range <0.01){scale = c(0:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  plot(c(1.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.25/max(shift)
  
  for(i in c(1:l)){
    for(i1 in c(1:length(groups[[i]]))){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
    }
  }
  
  y = max(unlist(groups))
  for(i in c(1:length(groups))){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i]<signif_threshold){
      pval1 = "*"
      y = max(unlist(groups[[i]]))+b
      # if(p_values[i] <signif_threshold/10){pval1 = "**"}
      # if(p_values[i] <signif_threshold/100){pval1 = "***"}
      # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
      text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
    }
  }
  
  
  dev.off()
  
  #################### plot by PCA group blood
  
  mat = m_overlap_sampled1
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  group_PCA_list = gsub("biopsy","blood",group_PCA_list)
  w = which(group_PCA_list %in% rownames(mat))
  mat_stat = mat[group_PCA_list[w],]
  nz = apply(mat_stat ,2,function(x){length(which(x!=0))})
  mat_stat  = mat_stat [,which(nz>=1)]
  
  factor = factor(factor_PCA[w])
  mat_stat1 = mat_stat^0.25
  fit = manova(formula = mat_stat1 ~ factor)
  p1 = summary.aov(fit)
  nam = gsub(" Response ","",names(p1))
  p_value = NULL
  means = NULL
  
  i1 = 0
  for(i in p1){
    i1 = i1+1
    p_value = c(p_value, i$'Pr(>F)'[1]) 
    if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
    }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
  }
  p_value[which(is.na(p_value))] = 2
  print(min(p_value))
  
  colnames(means) = paste("mean.group.", c(1:length(means[1,])))
  combined_p_value = cbind(p_value ,means)
  name = colnames(mat_stat)
  rownames(combined_p_value) = nam
  p.site = rep(concat(c("biopsy")), length(nam))
  p.analysis = rep(analysis, length(nam))
  x = cbind(p.site, p.analysis, combined_p_value)
  summary_stats = x
  
  analysis = "Clonal_overlap_B_cell_subsets"
  
  outfile = concat(c(output_directory,"Stats_",analysis,"_blood_", batch,".txt"))
  
  write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  groups = NULL
  for(i in c(1:length(mat_stat[1,]))){
    g1 = NULL
    for(s in c(1:length(groups_PCA))){
      g1 = c(g1, list(mat_stat[gsub("_biopsy","-blood",groups_PCA[[s]]),i]))}
    groups = c(groups, list(g1))
  }
  names(groups) = colnames(mat_stat)
  
  ######################################
  analysis = "Clonal_overlap_B_cell_subsets"
  fileout1=concat(c(output_directory,"",analysis,"_blood_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.8, width=w*3.7)
  par(mfrow= c(1,1), mar = c(20,4,2.5,1))
  p_values = p_value
  factors1 = paste("group",c(1:length(groups_PCA)))
  factors = names(groups)
  main = concat(c("clonal overlap:blood"))
  max = max(c(unlist(groups), unlist(groups))*1.2)
  min = 0
  b = (max-min)*0.035
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = min(c(max,100))
  range = max-min
  if(range>50){scale = c(0:100)*20}
  if(range<=50){scale = c(0:100)*10}
  if(range<=30){scale = c(0:100)*5}
  if(range <15){scale = c(0:100)*2.5}
  if(range <5){scale = c(0:100)*1}
  if(range <4){scale = c(0:100)*0.5}
  if(range <1.5){scale = c(0:1000)*0.2}
  if(range <0.5){scale = c(0:100)*0.1}
  if(range <0.1){scale = c(0:100)*0.01}
  if(range <0.01){scale = c(0:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.25/max(shift)
  
  for(i in c(1:l)){
    for(i1 in c(1:length(groups[[i]]))){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
    }
  }
  
  y = max(unlist(groups))
  for(i in c(1:length(groups))){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i]<signif_threshold){
      pval1 = "*"
      y = max(unlist(groups[[i]]))+b
      # if(p_values[i] <signif_threshold/10){pval1 = "**"}
      # if(p_values[i] <signif_threshold/100){pval1 = "***"}
      # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
      text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
    }
  }
  
  
  dev.off()
  
  ###################### plot network plot for clonal sharing
  pca_groups = paste("Group", c(1:length(groups_PCA)))
  mean_cell_type_proportions = matrix(data = 0,nrow = length(pca_groups), ncol = length(cell_types), dimnames = c(list(pca_groups), list(cell_types)))
  
  mean_cell_overlap = matrix(data = 0,nrow = length(pca_groups), ncol = length(colnames(m_overlap_sampled1)), dimnames = c(list(pca_groups), list(colnames(m_overlap_sampled1))))
  
  for(i in c(1:length(pca_groups))){
    mean_cell_type_proportions [i,]= apply(m_cell_types_proportions[gsub("_","-",groups_PCA[[i]]), ], 2, median)
    mean_cell_overlap [i,] = apply(m_overlap_sampled1[gsub("_","-",groups_PCA[[i]]), ], 2, mean)
  }
  
  
  ### 
  analysis = "Clonal_overlap_B_cell_subsets"
  file = concat(c(output_directory,"Stats_",analysis,"_biopsy_", batch,".txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p_value = as.numeric(p[,"p_value"])
  transition = rownames(p)
  names(p_value) = transition
  mean1 = as.numeric(p[,"mean.group..1"])
  mean2 = as.numeric(p[,"mean.group..2"])
  names(mean1) = transition
  names(mean2) = transition
  pca_group_means = c(list(mean1), list(mean2))
  transition_matrix = cbind(mean1, mean2)
  transition_split=strsplit(transition, " - ", fixed = T)
  transition1 = NULL
  transition2 = NULL
  for(i in c(1:length(transition_split))){
    transition1 = c(transition1, transition_split[[i]][[1]])
    transition2 = c(transition2, transition_split[[i]][[2]])
  }
  all_cell_types = cell_types
  
  ##################
  library(igraph)
  
  signif = which(p_value<0.05)
  
  fileout1=concat(c(output_directory,"",analysis,"_network_biopsy_", batch,".pdf"))
  w=2.7
  pdf(file=fileout1, height=w*1, width=w*3)
  par(mfrow= c(1,3), mar = c(1,1,3,1))
  
  c = 1
  mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
  main = concat(c("immuno-group ",c))
  g <- graph.empty( n=0, directed=FALSE)
  names_plot = all_cell_types
  names_plot = gsub("B cell  ","",names_plot)
  names_plot = gsub("B cell ","",names_plot)
  names_plot = gsub("activated pre-","activated\npre-",names_plot)
  names_plot = gsub("memory activated","memory\nactivated",names_plot)
  
   g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
  names <- V(g)$name
  ids <- 1:length(names)
  names(ids) <- names
  
  w = which(mean_edge_strength!=0)
  w = c(1:length(mean_edge_strength))
  edge_strength = mean_edge_strength[w]
  p_value_sub = p_value[w]
  from <- transition1[w]
  to <- transition2[w]
  w = intersect(which(from %in% names),which(to %in% names))
  p_value_sub = p_value_sub[w]
  edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
  g <- add.edges(g, t(edges), weight= edge_strength[w])
  signif_edges = names(which(p_value_sub<0.05))
  t = transition_matrix[signif_edges,]
  
  sizes = mean_cell_type_proportions[c, all_cell_types]
  sizes_scaled = sizes^0.3
  sizes_scaled = sizes_scaled*5
  
  V(g)$size<-sizes_scaled
  V(g)$label.cex<-0.5
  V(g)$name = names_plot
  V(g)$color = "grey"
    layout1 =layout_in_circle(g)
    
    edge_strength_plot = edge_strength^0.5
    edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
    
    col = rgb(0.5,0.5,0.5,alpha = 0.5)
    cols = rep(col, length(edge_strength_plot))
    names(cols) = names(p_value_sub)
    cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("blue",alpha = 0.65)
    cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, max))))] = add.alpha("red",alpha = 0.65)
    
    # col = "black"
    plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.5),ylim = range(layout1[,2]*1.5))
    
    ## Apply labels manually
    #Specify x and y coordinates of labels, adjust outward as desired
    x = layout1[,1]*1.5
    y = layout1[,2]*1.5
    
    #create vector of angles for text based on number of nodes 
    # (flipping the orientation of the words half way around so none appear 
    # upside down)
    angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                   90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
    
    #Apply the text labels with a loop with angle as srt
    for (i in 1:length(x)) {
      text(x=x[i], y=y[i], labels=V(g)$name[i], adj=NULL, 
           pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
    }
    dev.off()
    
    
    ####### for blood
    ### 
    analysis = "Clonal_overlap_B_cell_subsets"
    file = concat(c(output_directory,"Stats_",analysis,"_blood_", batch,".txt"))
    p <- as.matrix(read.csv(file, head=T, sep="\t"))
    p_value = as.numeric(p[,"p_value"])
    transition = rownames(p)
    names(p_value) = transition
    mean1 = as.numeric(p[,"mean.group..1"])
    mean2 = as.numeric(p[,"mean.group..2"])
    names(mean1) = transition
    names(mean2) = transition
    pca_group_means = c(list(mean1), list(mean2))
    transition_matrix = cbind(mean1, mean2)
    transition_split=strsplit(transition, " - ", fixed = T)
    transition1 = NULL
    transition2 = NULL
    for(i in c(1:length(transition_split))){
      transition1 = c(transition1, transition_split[[i]][[1]])
      transition2 = c(transition2, transition_split[[i]][[2]])
    }
    all_cell_types = cell_types
    
    ##################
    library(igraph)
    
    signif = which(p_value<0.05)
    
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/BCR_TCR/",analysis,"_network_blood_", batch,".pdf"))
    w=2.7
    pdf(file=fileout1, height=w*1, width=w*3)
    par(mfrow= c(1,3), mar = c(1,1,3,1))
    
    c = 1
    mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
    main = concat(c("immuno-group ",c))
    g <- graph.empty( n=0, directed=FALSE)
    names_plot = all_cell_types
    names_plot = gsub("B cell  ","",names_plot)
    names_plot = gsub("B cell ","",names_plot)
    names_plot = gsub("activated pre-","activated\npre-",names_plot)
    names_plot = gsub("memory activated","memory\nactivated",names_plot)
    g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
    names <- V(g)$name
    ids <- 1:length(names)
    names(ids) <- names
    
    w = which(mean_edge_strength!=0)
    w = c(1:length(mean_edge_strength))
    edge_strength = mean_edge_strength[w]
    p_value_sub = p_value[w]
    from <- transition1[w]
    to <- transition2[w]
    w = intersect(which(from %in% names),which(to %in% names))
    p_value_sub = p_value_sub[w]
    edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
    g <- add.edges(g, t(edges), weight= edge_strength[w])
    signif_edges = names(which(p_value_sub<0.05))
    t = transition_matrix[signif_edges,]
    
    sizes = mean_cell_type_proportions[c, all_cell_types]
    sizes_scaled = sizes^0.3
    sizes_scaled = sizes_scaled*5
    
    V(g)$size<-sizes_scaled
    V(g)$label.cex<-0.5
    V(g)$name = names_plot
    V(g)$color = "grey"
      layout1 =layout_in_circle(g)
      edge_strength_plot = edge_strength^0.5
      edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
      
      col = rgb(0.5,0.5,0.5,alpha = 0.5)
      cols = rep(col, length(edge_strength_plot))
      names(cols) = names(p_value_sub)
      if(length(signif_edges)<1){
        cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("blue",alpha = 0.65)
        cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, max))))] = add.alpha("red",alpha = 0.65)
      }
      
      plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.5),ylim = range(layout1[,2]*1.5))
      
      ## Apply labels manually
      #Specify x and y coordinates of labels, adjust outward as desired
      x = layout1[,1]*1.5
      y = layout1[,2]*1.5
      
      #create vector of angles for text based on number of nodes 
      # (flipping the orientation of the words half way around so none appear 
      # upside down)
      angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                     90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
      
      #Apply the text labels with a loop with angle as srt
      for (i in 1:length(x)) {
        text(x=x[i], y=y[i], labels=V(g)$name[i], adj=NULL, 
             pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
      }
      dev.off()
}
BCR_clonal_overlap_per_site(VDJ_list_BCR, output_directory,groups_PCA)

### TCR overlap within samples between subsets
TCR_clonal_overlap_per_site<-function(VDJ_list_TCR, output_directory,CD48, groups_PCA){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type = VDJ_list[[2]]
  cell_type_broad = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone1 = VDJ[,"clone1"]
  clone2 = VDJ[,"clone2"]
  clone = clone1
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type1 = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    w_cd48 = names(cell_type)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[w_cd48])
    w_clone = intersect(names(which(clone!='-')),w_cd48)
    chain = names(CD48)[cd48]
    cell_types1 = cell_types_use

    ## cell counts per cell type
    w_clone = intersect(names(which(clone1!='-')),names(which(clone2!='-')))
    w_clone = intersect(w_clone,w_cd48)
    cell_types = sort(unique(cell_type[w_clone]))
    m_cell_types = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
    m_cell_types_proportions = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
    for(s in c(1:length(samples))){
      w = intersect(w_clone , names(which(sample == samples[s])))
      t = table(cell_type[w])
      m_cell_types[samples[s], names(t)]= t
      m_cell_types_proportions[samples[s], names(t)]= t*100/sum(t)
    }
    totals = rowSums(m_cell_types)
    
    ############# get overlap between compartments
    c1 = NULL
    c2 = NULL
    for(i1 in c(1:length(cell_types1))){
      for(i2 in c(i1:length(cell_types1))){
        if(i1<i2){
          c1 = c(c1, cell_types1[i1])
          c2 = c(c2, cell_types1[i2])}}}
    sharing_names = apply(cbind(c1,c2), 1, paste,collapse = " - ")
    
    print("Running without subsampling")
    m_overlap_raw = matrix(data = 0,nrow = length(samples), ncol = length(sharing_names), dimnames = c(list(samples), list(sharing_names)))
    
    for(s in c(1:length(samples))){
      print(s)
      w_sample = intersect(w_clone , names(which(sample == samples[s])))
      t = table(clone[w_sample], cell_type[w_sample])
      a = apply(t, 1, function(x){length(which(x!=0))})
      t[which(a>1),]
      for(i in c(1:length(c1))){
        w1 = intersect(w_sample, names(which(cell_type==c1[i])))
        w2 = intersect(w_sample, names(which(cell_type==c2[i])))
        shared_clones = intersect(clone[w1],clone[w2])
        m_overlap_raw[samples[s], sharing_names[i]] = length(shared_clones)
      }
    }
    
    list_overlapping_per_sample = NULL
    for(s in c(1:length(samples))){
      m_overlapping = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
      m_overlapping_total = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
      wa = intersect(names(sample)[which(sample== samples[s])], w_clone)
      t = table(clone[wa], cell_type[wa])
      w = which(apply(t, 1, function(x){length(which(x!=0))})>1)
      if(length(w)>0){
        # print(r)
        t=rbind(t[w,])
        for(i in c(1:length(w))){
          x = names(t[i,][which(t[i,]!=0)])
          for(i1 in c(1:(length(x)-1))){
            for(i2 in c((i1+1):length(x))){
              m_overlapping_total[x[i1],x[i2]] = m_overlapping_total[x[i1],x[i2]]+1
            }}
        }
      }
      
      scale = NULL
      if(sum(m_overlapping_total)==0){scale = 1
      m_overlapping_total = m_overlapping_total#/sum(m_overlapping_total)
      }else{m_overlapping_total = m_overlapping_total/sum(m_overlapping_total)}
      
      scale = 1
      list_overlapping_per_sample = c(list_overlapping_per_sample, list(m_overlapping_total*mean(scale)))
      print (s)
    }
    names(list_overlapping_per_sample) = samples
    
    groups = NULL
    links = sharing_names
    
    overlapping_per_link = matrix(data = 0,nrow = length(samples), ncol = length(links), dimnames = c(list(samples), list(links)))
    for(s in c(1:length(samples))){
      mat = list_overlapping_per_sample[[s]]
      for(c1 in c(1:length(cell_types))){
        for(c2 in c(c1:length(cell_types))){
          if(c1<c2){
            if(is.na(mat[cell_types[c1],cell_types[c2]])==F){
              link = concat(c(cell_types[c1], " - ", cell_types[c2]))
              overlapping_per_link[samples[s],link] = mat[cell_types[c1],cell_types[c2]]
            }}}}}
    
    m_overlap_sampled1 = overlapping_per_link
    saveRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")), m_overlap_sampled1)
    print(concat(c(type1, " done")))
  }
  
  print("Plots and statistics running")
  
  #################### plot by PCA group
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type1 = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    w_cd48 = names(sample)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[wx])
    w_clone = intersect(names(which(clone!='-')),w_cd48)
    chain = names(CD48)[cd48]
    cell_types1 = cell_types_use
    
    
    m_overlap_sampled1  =readRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")))
    
    Means_factor = function(factor, x){
      m = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        m = c(m, mean(x1))}
      return(m)}
    
    mat = m_overlap_sampled1
    nz = apply(mat,2,function(x){length(which(x!=0))})
    mat = mat[,which(nz>=1)]
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]]))
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    w = which(group_PCA_list %in% rownames(mat))
    mat_stat = mat[group_PCA_list[w],]
    
    factor = factor(factor_PCA[w])
    mat_stat1 = mat_stat^1
    fit = manova(formula = mat_stat1 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
      }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    names(p_value) =nam
    print(head(sort(p_value), 10))
    analysis= concat(c("Clonal_overlap_T_cell_subsets_",type1))
    colnames(means) = paste("mean.group.", c(1:length(means[1,])))
    combined_p_value = cbind(p_value ,means)
    name = colnames(mat_stat)
    rownames(combined_p_value) = nam
    p.site = rep(concat(c("biopsy")), length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_stats = x
    
    
    groups = NULL
    for(i in c(1:length(mat_stat[1,]))){
      g1 = NULL
      for(s in c(1:length(groups_PCA))){
        g1 = c(g1, list(mat_stat[gsub("_","-",groups_PCA[[s]]),i]))}
      groups = c(groups, list(g1))
    }
    names(groups) = colnames(mat_stat)
    
    ######################################
    analysis = "Clonal_overlap_T_cell_subsets"
    
    outfile = concat(c(output_directory,"Stats_",analysis,"_biopsy_", type1,"_", batch,".txt"))
    
    write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    summary_stats_biopsy = summary_stats
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    fileout1=concat(c(output_directory,"",analysis,"_biopsy_", type1,"_", batch,".pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1.8, width=w*8.7)
    par(mfrow= c(1,1), mar = c(20,4,2.5,1))
    p_values = p_value
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = names(groups)
    main = concat(c("clonal overlap:biopsy"))
    max = max(c(unlist(groups), unlist(groups))*1.2)
    min = 0
    b = (max-min)*0.035
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 200#length(groups)
    max_scale = min(c(max,100))
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*10}
    if(range<=30){scale = c(0:100)*5}
    if(range <15){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(3, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.18
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:length(groups[[i]]))){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
      }
    }
    
    y = max(unlist(groups))
    for(i in c(1:length(groups))){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        y = max(unlist(groups[[i]]))+b
        # if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
        text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
      }
    }
    
    dev.off()
    
    #################### plot by PCA group blood
    
    mat = m_overlap_sampled1
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    group_PCA_list = gsub("_biopsy","-blood",group_PCA_list)
    w = which(group_PCA_list %in% rownames(mat))
    mat_stat = mat[group_PCA_list[w],]
    nz = apply(mat_stat ,2,function(x){length(which(x!=0))})
    mat_stat  = mat_stat [,which(nz>=2)]
    
    factor = factor(factor_PCA[w])
    mat_stat1 = mat_stat^1
    fit = manova(formula = mat_stat1 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
      }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    print(min(p_value))
    
    colnames(means) = paste("mean.group.", c(1:length(means[1,])))
    combined_p_value = cbind(p_value ,means)
    name = colnames(mat_stat)
    rownames(combined_p_value) = nam
    p.site = rep(concat(c("biopsy")), length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_stats = x
    
    analysis = "Clonal_overlap_T_cell_subsets"
    
    outfile = concat(c(output_directory,"Stats_",analysis,"_blood_", type1,"_", batch,".txt"))
    
    write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    
    groups = NULL
    for(i in c(1:length(mat_stat[1,]))){
      g1 = NULL
      for(s in c(1:length(groups_PCA))){
        g1 = c(g1, list(mat_stat[gsub("_biopsy","-blood",groups_PCA[[s]]),i]))}
      groups = c(groups, list(g1))
    }
    names(groups) = colnames(mat_stat)
    
    ######################################
    analysis = "Clonal_overlap_T_cell_subsets"
    fileout1=concat(c(output_directory,"",analysis,"_blood_", type1,"_", batch,".pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1.8, width=w*8.7)
    par(mfrow= c(1,1), mar = c(20,4,2.5,1))
    p_values = p_value
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = names(groups)
    main = concat(c("clonal overlap:blood"))
    max = max(c(unlist(groups), unlist(groups))*1.2)
    min = 0
    b = (max-min)*0.035
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 80#length(groups)
    max_scale = min(c(max,100))
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*10}
    if(range<=30){scale = c(0:100)*5}
    if(range <15){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(3, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.18
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:length(groups[[i]]))){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
      }
    }
    
    y = max(unlist(groups))
    for(i in c(1:length(groups))){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        y = max(unlist(groups[[i]]))+b
        # if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
        text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.7)
      }
    }
    
    dev.off()
    
    ###################### plot network plot for clonal sharing
    pca_groups = paste("Group", c(1:length(groups_PCA)))
    mean_cell_type_proportions = matrix(data = 0,nrow = length(pca_groups), ncol = length(cell_types1), dimnames = c(list(pca_groups), list(cell_types1)))
    for(i in c(1:length(pca_groups))){
      mean_cell_type_proportions [i,colnames(m_cell_types_proportions)]= apply(m_cell_types_proportions[gsub("_","-",groups_PCA[[i]]), ], 2, median)
    }
    
    
    ### 
    analysis = "Clonal_overlap_T_cell_subsets"
    file = concat(c(output_directory,"Stats_",analysis,"_biopsy_", type1,"_", batch,".txt"))
    p <- as.matrix(read.csv(file, head=T, sep="\t"))
    p_value = as.numeric(p[,"p_value"])
    transition = rownames(p)
    names(p_value) = transition
    mean1 = as.numeric(p[,"mean.group..1"])
    mean2 = as.numeric(p[,"mean.group..2"])
    names(mean1) = transition
    names(mean2) = transition
    pca_group_means = c(list(mean1), list(mean2))
    transition_matrix = cbind(mean1, mean2)
    transition_split=strsplit(transition, " - ", fixed = T)
    transition1 = NULL
    transition2 = NULL
    for(i in c(1:length(transition_split))){
      transition1 = c(transition1, transition_split[[i]][[1]])
      transition2 = c(transition2, transition_split[[i]][[2]])
    }
    all_cell_types = unique(c(transition1,transition2))
    
    ##################
    library(igraph)
    
    signif = which(p_value<0.05)
    
    c = 1
    mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
    main = concat(c("immuno-group ",c))
    g <- graph.empty( n=0, directed=FALSE)
    if(cd48==1){
      all_cell_types_broad = gsub("Activated ", "", all_cell_types)
      order = order(all_cell_types_broad)
      all_cell_types = all_cell_types[order]}
    # order= c(grep("CD4", all_cell_types),grep("Treg", all_cell_types),grep("Tfh", all_cell_types)
    # order = intersect(cell_classifications [which(cell_classifications[,1]=="CD4"), 2], all_cell_types)
    # all_cell_types = c(order , setdiff(all_cell_types, order))
    # all_cell_types = setdiff(all_cell_types, "MAIT")
    names_plot = all_cell_types
    names_plot = gsub("_"," ",names_plot)
    names_plot = gsub(" EffectorMem","\nEffectorMem",names_plot)
    names_plot = gsub("Exhausted ","Exhausted\n",names_plot)
    names_plot = gsub("chemokine ","chemokine\n",names_plot)
    
    g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
    names <- V(g)$name
    ids <- 1:length(names)
    names(ids) <- names
    
    w = which(mean_edge_strength!=0)
    w = c(1:length(mean_edge_strength))
    edge_strength = mean_edge_strength[w]
    p_value_sub = p_value[w]
    from <- transition1[w]
    to <- transition2[w]
    w = intersect(which(from %in% names),which(to %in% names))
    p_value_sub = p_value_sub[w]
    edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
    g <- add.edges(g, t(edges), weight= edge_strength[w])
    signif_edges = names(which(p_value_sub<0.05))
    t = transition_matrix[signif_edges,]
    
    sizes = mean_cell_type_proportions[c, all_cell_types]
    sizes_scaled = sizes^0.3
    sizes_scaled = sizes_scaled*5
    
    V(g)$size<-sizes_scaled
    V(g)$label.cex<-0.5
    V(g)$name = names_plot
    V(g)$color = "grey"
      layout1 =layout_in_circle(g)
      edge_strength_plot = edge_strength[w]^0.5
      edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
      
      col = rgb(0.5,0.5,0.5,alpha = 0.5)
      cols = rep(col, length(edge_strength_plot))
      names(cols) = names(p_value_sub)
      c = 1
      cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, min))))] = add.alpha("blue",alpha = 0.65)
      cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, max))))] = add.alpha("red",alpha = 0.65)
      #cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("purple",alpha = 0.65)
      
      fileout1=concat(c(output_directory,"",analysis,"_network_biopsy_", type1,"_", batch,".pdf"))
      w=2.7
      pdf(file=fileout1, height=w*1, width=w*3)
      par(mfrow= c(1,3), mar = c(2,1,2,1))
      # col = "black"
      plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.7),ylim = range(layout1[,2]*1.7))
      
      ## Apply labels manually
      #Specify x and y coordinates of labels, adjust outward as desired
      x = layout1[,1]*1.65
      y = layout1[,2]*1.65
      
      #create vector of angles for text based on number of nodes 
      # (flipping the orientation of the words half way around so none appear 
      # upside down)
      angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                     90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
      
      #Apply the text labels with a loop with angle as srt
      lab = gsub("T cell ","",V(g)$name)
      lab = gsub("Activated ","Activated\n", lab)
      lab = gsub("NK cell ","", lab)
      lab = gsub("CD4","CD4\n", lab)
      lab = gsub("CD8","CD8\n", lab)
      
      for (i in 1:length(x)) {
        text(x=x[i], y=y[i], labels=lab[i], adj=NULL, 
             pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
      }
      dev.off()
      
      
      ####### for blood
      ### 
      analysis = "Clonal_overlap_T_cell_subsets"
      file = concat(c(output_directory,"Stats_",analysis,"_blood_", type1,"_", batch,".txt"))
      p <- as.matrix(read.csv(file, head=T, sep="\t"))
      p_value = as.numeric(p[,"p_value"])
      transition = rownames(p)
      names(p_value) = transition
      mean1 = as.numeric(p[,"mean.group..1"])
      mean2 = as.numeric(p[,"mean.group..2"])
      names(mean1) = transition
      names(mean2) = transition
      pca_group_means = c(list(mean1), list(mean2))
      transition_matrix = cbind(mean1, mean2)
      transition_split=strsplit(transition, " - ", fixed = T)
      transition1 = NULL
      transition2 = NULL
      for(i in c(1:length(transition_split))){
        transition1 = c(transition1, transition_split[[i]][[1]])
        transition2 = c(transition2, transition_split[[i]][[2]])
      }
      all_cell_types = unique(c(transition1,transition2))
      
      ##################
      library(igraph)
      
      signif = which(p_value<0.05)
      if(cd48==1){
        all_cell_types_broad = gsub("Activated ", "", all_cell_types)
        order = order(all_cell_types_broad)
        all_cell_types = all_cell_types[order]}
      
      
      c = 1
      mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
      main = concat(c("immuno-group ",c))
      g <- graph.empty( n=0, directed=FALSE)
      names_plot = all_cell_types
      names_plot = gsub("_"," ",names_plot)
      names_plot = gsub(" EffectorMem","\nEffectorMem",names_plot)
      names_plot = gsub("Exhausted ","Exhausted\n",names_plot)
      names_plot = gsub("chemokine ","chemokine\n",names_plot)
      g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
      names <- V(g)$name
      ids <- 1:length(names)
      names(ids) <- names
      
      w = which(mean_edge_strength!=0)
      w = c(1:length(mean_edge_strength))
      edge_strength = mean_edge_strength[w]
      p_value_sub = p_value[w]
      from <- transition1[w]
      to <- transition2[w]
      w = intersect(which(from %in% names),which(to %in% names))
      p_value_sub = p_value_sub[w]
      edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
      g <- add.edges(g, t(edges), weight= edge_strength[w])
      signif_edges = names(which(p_value_sub<0.05))
      t = transition_matrix[signif_edges,]
      
      sizes = mean_cell_type_proportions[c, all_cell_types]
      sizes_scaled = sizes^0.3
      sizes_scaled = sizes_scaled*5
      
      V(g)$size<-sizes_scaled
      V(g)$label.cex<-0.5
      V(g)$name = names_plot
      V(g)$color = "grey"
        layout1 =layout_in_circle(g)
        edge_strength_plot = edge_strength[w]^0.5
        edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
        
        col = rgb(0.5,0.5,0.5,alpha = 0.5)
        cols = rep(col, length(edge_strength_plot))
        names(cols) = names(p_value_sub)
        c = 1
        cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, min))))] = add.alpha("blue",alpha = 0.65)
        cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, max))))] = add.alpha("red",alpha = 0.65)
        #cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("purple",alpha = 0.65)
        
        fileout1=concat(c(output_directory,"",analysis,"_network_blood_", type1,"_", batch,".pdf"))
        w=2.7
        pdf(file=fileout1, height=w*1, width=w*3)
        par(mfrow= c(1,3), mar = c(2,1,2,1))
        # col = "black"
        plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.7),ylim = range(layout1[,2]*1.7))
        
        ## Apply labels manually
        #Specify x and y coordinates of labels, adjust outward as desired
        x = layout1[,1]*1.65
        y = layout1[,2]*1.65
        
        #create vector of angles for text based on number of nodes 
        # (flipping the orientation of the words half way around so none appear 
        # upside down)
        angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                       90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
        
        #Apply the text labels with a loop with angle as srt
        lab = gsub("T cell ","",V(g)$name)
        lab = gsub("Activated ","Activated\n", lab)
        lab = gsub("NK cell ","", lab)
        
        for (i in 1:length(x)) {
          text(x=x[i], y=y[i], labels=lab[i], adj=NULL, 
               pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
        }
        dev.off()
  }
  
}
TCR_clonal_overlap_per_site(VDJ_list_TCR, output_directory,CD48, groups_PCA)


################# clonal overlap between blood and tumour
BCR_immunosurveillance<-function(VDJ_list_BCR, output_directory, batch, groups_PCA){
  VDJ_list = VDJ_list_BCR
  type = "BCR"
  type1 = "B_cells"
  analysis = "Clonal_overlap_between_sites"
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[6]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  w_clone = names(which(clone!='-'))
  clone = VDJ[,"clone2"]
  clone1 = VDJ[,"clone1"]
  clone2 = VDJ[,"clone2"]
  cdr3_1 = VDJ[,"cdr3_aa1"]
  cdr3_2 = VDJ[,"cdr3_aa2"]
  names(cdr3_1) = rownames(VDJ)
  names(cdr3_2) = rownames(VDJ)
  names(clone) = rownames(VDJ)
  names(clone1) = rownames(VDJ)
  names(clone2) = rownames(VDJ)
  constant_region1 = VDJ[,"constant_region1"]
  names(constant_region1) = rownames(VDJ)
  
  V_gene_10X1 = VDJ[,"V_gene_10X1"]
  V_mm1 = VDJ[,"V_mm1"]
  names(constant_region1) = rownames(VDJ)
  names(V_gene_10X1) = rownames(VDJ)
  names(V_mm1) = rownames(VDJ)
  constant_region1s = sort(unique(constant_region1))
  V_gene_10X1s = sort(unique(V_gene_10X1))
  cell_ids = rownames(VDJ)
  id = cell_ids
  #### get subsampling threshold
  min = 1
  
  #### 
  m_isotypes_shared = matrix(data = -1,nrow = length(samples), ncol = length(constant_region1s), dimnames = c(list(samples), list(constant_region1s)))
  m_isotypes_private = matrix(data = -1,nrow = length(samples), ncol = length(constant_region1s), dimnames = c(list(samples), list(constant_region1s)))
  m_v_genes_shared = matrix(data = -1,nrow = length(samples), ncol = length(V_gene_10X1s), dimnames = c(list(samples), list(V_gene_10X1s)))
  m_v_genes_private = matrix(data = -1,nrow = length(samples), ncol = length(V_gene_10X1s), dimnames = c(list(samples), list(V_gene_10X1s)))	
  m_cell_type_shared = matrix(data = -1,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
  m_cell_type_private = matrix(data = -1,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))	
  m_CDR3_length_chain1_p = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types), dimnames = c(list(Patients), list(cell_types)))	
  m_CDR3_length_chain1_sh = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types), dimnames = c(list(Patients), list(cell_types)))	
  m_CDR3_length_chain2_p = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types), dimnames = c(list(Patients), list(cell_types)))	
  m_CDR3_length_chain2_sh = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types), dimnames = c(list(Patients), list(cell_types)))	
  proportion_shared_cells = rep(-1, length(samples))
  names(proportion_shared_cells)=samples
  SHM_IGHV4_34_overlapping = NULL
  nam1 = NULL
  shared_private_clone_sizes = NULL
  
  ###
  cell_ids = intersect(names(clone1), id)
  w_clone1 = intersect(intersect(names(clone1)[which(clone1!='-')], names(clone1)[which(clone1!='0')]), id)
  w_clone2 = intersect(intersect(names(clone2)[which(clone2!='-')], names(clone2)[which(clone2!='0')]), id)
  shared_ids = rep("UNKNOWN",length(cell_ids))
  shared_ids[which(Sample.Type=="biopsy")] = "UNKNOWN TIL"
  shared_ids[which(Sample.Type=="blood")] = "UNKNOWN blood"
  ### consider shared heavy and/or light chain
  names(shared_ids) = cell_ids
  print("Getting overlap")
  for(s in c(1:length(Patients))){
    print(s)
    w1 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))]
    w2 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))]
    w1 = intersect(w1, cell_ids[which(cell_ids %in% w_clone1)])
    w2 = intersect(w2, cell_ids[which(cell_ids %in% w_clone1)])
    shared_clones = intersect(clone1[w1] , clone1[w2])
    cl1a = intersect(names(clone1 [which(clone1 %in% shared_clones)]), w1)
    cl2a = intersect(names(clone1 [which(clone1 %in% shared_clones)]), w2)
    
    w1 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))]
    w2 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))]
    w1 = intersect(w1, cell_ids[which(cell_ids %in% w_clone2)])
    w2 = intersect(w2, cell_ids[which(cell_ids %in% w_clone2)])
    shared_clones = intersect(clone2[w1] , clone2[w2])
    cl1b = intersect(names(clone2 [which(clone2 %in% shared_clones)]), w1)
    cl2b = intersect(names(clone2 [which(clone2 %in% shared_clones)]), w2)
    cls = c(list(c(cl1a, cl1b)),list(c(cl2a, cl2b))) ## shared clones. light chain
    
    tx = table(shared_ids)
    proportion_shared_cells[concat(c(Patients[s],"-blood"))] = tx["immunosurveilling blood"]*100/sum(tx[c("immunosurveilling blood","private blood")])
    proportion_shared_cells[concat(c(Patients[s],"-biopsy"))] = tx["immunosurveilling TIL"]*100/sum(tx[c("immunosurveilling TIL","private TIL")])
    
    shared_ids[cl1a] = "immunosurveilling TIL"
    shared_ids[cl2a] = "immunosurveilling blood"
    shared_ids[cl1b] = "immunosurveilling TIL"
    shared_ids[cl2b] = "immunosurveilling blood"
    
    w1 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))]
    w2 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))]
    pl1 = setdiff(w1, unlist(cls))
    pl2 = setdiff(w2, unlist(cls))
    pls = c(list(pl1),list(pl2)) ## private clones
    shared_ids[pl1] = "private TIL"
    shared_ids[pl2] = "private blood"
    
    t_biopsy = length(w1)
    t_blood  = length(w2)
    shared_clones_biopsy = as.numeric((table(clone [cls[[1]]])-1)*100/t_biopsy)
    shared_clones_blood = as.numeric ((table(clone [cls[[2]]])-1)*100/t_blood)
    private_clones_biopsy = as.numeric ((table(clone [pls[[1]]])-1)*100/t_biopsy )
    private_clones_blood = as.numeric ((table(clone [pls[[2]]])-1)*100/t_blood )
    
    # private clones
    cdr31 = nchar( cdr3_1[pls[[1]]])
    cdr32 = nchar( cdr3_2[pls[[1]]])
    cdr31 = cdr31[which(is.na(cdr31)==F)]
    cdr32 = cdr32[which(is.na(cdr32)==F)]
    t1a = table(cell_type[names(cdr31)], cdr31)
    t1b = table(cell_type[names(cdr32)], cdr32)
    
    cdr31 = nchar( cdr3_1[cls[[1]]])
    cdr32 = nchar( cdr3_2[cls[[1]]])
    cdr31 = cdr31[which(is.na(cdr31)==F)]
    cdr32 = cdr32[which(is.na(cdr32)==F)]
    t2a = table(cell_type[names(cdr31)], cdr31)
    t2b = table(cell_type[names(cdr32)], cdr32)
    
    ct_use = names(which(rowSums(t2a)>=5))
    apply(t1a, 1, function(x){ sum(x*as.numeric(colnames(t1a))/sum(x))})
    
    m_CDR3_length_chain1_p[Patients[s], ct_use] = apply(t1a, 1, function(x){ sum(x*as.numeric(colnames(t1a))/sum(x))})[ct_use]
    m_CDR3_length_chain1_sh[Patients[s], ct_use] = apply(t2a, 1, function(x){ sum(x*as.numeric(colnames(t2a))/sum(x))})[ct_use]
    m_CDR3_length_chain2_p[Patients[s], ct_use] = apply(t1b, 1, function(x){ sum(x*as.numeric(colnames(t1b))/sum(x))})[ct_use]
    m_CDR3_length_chain2_sh[Patients[s], ct_use] = apply(t2b, 1, function(x){ sum(x*as.numeric(colnames(t2b))/sum(x))})[ct_use]
    
    g = c(list(private_clones_blood),list(shared_clones_blood),list(private_clones_biopsy),list(shared_clones_biopsy))
    names(g) = c("private blood","immunosurv. blood","private tumour","immunosurv. tumour")
    # boxplot(g)
    # t.test(g[[3]],y=g[[4]])
    shared_private_clone_sizes = c(shared_private_clone_sizes, list(g))
    
    
    for(i in c(1:length(cls))){
      if(length(cls[[i]])>min){
        sample1 = concat(c(Patients[s],"-",Sample.Types[i]))
        t = table(constant_region1[cls[[i]]])
        m_isotypes_shared[sample1, ]=0
        m_isotypes_shared[sample1, names(t)] = t*100/sum(t)
        t = table(V_gene_10X1[cls[[i]]])
        m_v_genes_shared[sample1, ]=0
        m_v_genes_shared[sample1, names(t)] = t*100/sum(t)
        t = table(cell_type[cls[[i]]])
        m_cell_type_shared[sample1, ]=0
        m_cell_type_shared[sample1, names(t)] = t*100/sum(t)
        cls_IGHV434 = names(which(V_gene_10X1[cls[[i]]] == "IGHV4-34"))
        SHM_IGHV4_34_overlapping  = c(SHM_IGHV4_34_overlapping , list(V_mm1 [cls_IGHV434]))
        nam1 = c(nam1,sample1)
      }
      if(length(pls[[i]])>min){
        sample1 = concat(c(Patients[s],"-",Sample.Types[i]))
        t = table(constant_region1[pls[[i]]])
        m_isotypes_private[sample1, ]=0
        m_isotypes_private[sample1, names(t)] = t*100/sum(t)
        t = table(V_gene_10X1[pls[[i]]])
        m_v_genes_private[sample1, ]=0
        m_v_genes_private[sample1, names(t)] = t*100/sum(t)
        t = table(cell_type[pls[[i]]])
        m_cell_type_private[sample1, ]=0
        m_cell_type_private[sample1, names(t)] = t*100/sum(t)
        cls_IGHV434 = names(which(V_gene_10X1[cls[[i]]] == "IGHV4-34"))
        SHM_IGHV4_34_overlapping  = c(SHM_IGHV4_34_overlapping , list(V_mm1 [cls_IGHV434]))
        nam1 = c(nam1,sample1)
      }
    }
  }
  names(SHM_IGHV4_34_overlapping) = nam1
  names(shared_private_clone_sizes) = Patients
  
  analysis = "cell_type_private"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(cbind(m_cell_type_private), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "cell_type_shared"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(cbind(m_cell_type_shared), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "isotype_private"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(cbind(m_isotypes_private), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "isotype_shared"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(cbind(m_isotypes_shared), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "V_genes_private"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(cbind(m_v_genes_private), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "V_genes_shared"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(cbind(m_v_genes_shared), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  
  ########################### plot clone sizes between private and shared in blood and biopsy
  analysis="Clone_sizes_private_immunosurvielling_BCR"
  fileout1=concat(c(output_directory,"Overlap_",analysis,"_boxplots_", batch,".pdf"))
  w=2.35
  pdf(file=fileout1, height=w*1.2*3, width=w*0.6*6)
  par(mfrow= c(3,6), mar = c(11,5,4,0.1))
  
  array1 = NULL
  array2 = NULL
  array3 = NULL
  array4 = NULL
  for(s in c(1:length(Patients))){
    groups = shared_private_clone_sizes[[s]]
    array1 = c(array1, mean(groups[[1]]))
    array2 = c(array2, mean(groups[[2]]))
    array3 = c(array3, mean(groups[[3]]))
    array4 = c(array4, mean(groups[[4]]))
    factors = names(groups)
    main = Patients[s]
    max = max(c(unlist(groups), unlist(groups))*1.25)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = length(factors)
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.3
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.3/max(shift)
    
    library(RColorBrewer)
    cols = add.alpha (brewer.pal(8, "Set1")[c(1,1,2,2)], alpha = 0.5)
    cols[c(2,4)] = add.alpha (cols[c(2,4)], alpha = 0.5)
    #cols = add.alpha (brewer.pal(8, "Set1")[c(1,2,1,2)], alpha = 0.5)
    cols1 = add.alpha (c(cols,alpha = 0.5))
    #cols2 =  brewer.pal(8, "Set1")[c(1,2,1,2)]
    
    for(i in c(1:l)){
      points1=as.numeric(groups[[i]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i,width,cols[i],1, cols1[i])
      points(rep(i, length(points1)),points1, pch =21, col=cols[i], cex = 0.7)
    }
    if(min(c(length(groups[[1]]),length(groups[[2]]),length(groups[[3]]),length(groups[[4]])))>3){
      pval1 = t.test(groups[[1]], groups[[2]])$p.value
      pval2 = t.test(groups[[3]], groups[[4]])$p.value
      pvals = c(pval1, pval2)
      start = c(1,3)
      for(i in c(1:length(start))){	
        b = max*0.035
        signif_threshold = 0.05
        pval1 = signif(pvals[i], digits = 2)
        y = max(unlist(groups[c(start[i], start[i]+1)]))
        y = y+1*b
        segments(start[i],y+b, start[i]+1,y+b,lwd = 3, col = "darkgrey")
        text(mean(c(start[i], start[i]+1)), y+2.5*b, labels = pval1, cex = cex)
      }
    }
  }
  
  mean_shared_private_clone_sizes = cbind(array1, array2, array3, array4)
  colnames(mean_shared_private_clone_sizes) = names(shared_private_clone_sizes[[1]])
  rownames(mean_shared_private_clone_sizes) = Patients
  
  
  groups_PCA1 = c(groups_PCA, list(unlist(groups_PCA)))
  names(groups_PCA1) = c("ME", "AE", "All")
  for(s in c(1:length(groups_PCA1))){
    mat = mean_shared_private_clone_sizes[gsub("_biopsy","",groups_PCA1[[s]]),]
    groups = NULL
    for(i in c(1:length(mat[1,]))){
      x = mat[,i]
      x = x[which(is.na(x)==F)]
      groups = c(groups, list(x))
    }
    factors = colnames(mat)
    main = names(groups_PCA1)[s]
    if(s==4){main = "all"}
    max = max(mean_shared_private_clone_sizes[which(is.na(mean_shared_private_clone_sizes)==F)])*1.25
    min = 0
    b = (max-min)*0.034
    ylab = "mean clone size (%)"
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = length(factors)
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.3
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.3/max(shift)
    
    library(RColorBrewer)
    cols = add.alpha (brewer.pal(8, "Set1")[c(1,1,2,2)], alpha = 0.5)
    cols[c(2,4)] = add.alpha (cols[c(2,4)], alpha = 0.5)
    #cols = add.alpha (brewer.pal(8, "Set1")[c(1,2,1,2)], alpha = 0.5)
    cols1 = add.alpha (c(cols,alpha = 0.5))
    #cols2 =  brewer.pal(8, "Set1")[c(1,2,1,2)]
    
    for(i in c(1:l)){
      points1=as.numeric(groups[[i]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i,width,cols[i],1, cols1[i])
      points(rep(i, length(points1)),points1, pch =21, col=cols[i], cex = 0.7)
    }
    if(min(c(length(groups[[1]]),length(groups[[2]]),length(groups[[3]]),length(groups[[4]])))>=3){
      pval1 = wilcox.test(groups[[1]],groups[[2]],paired = F)$p.value
      pval2 = wilcox.test(groups[[3]], groups[[4]],paired = F)$p.value
      pvals = c(pval1, pval2)
      start = c(1,3)
      for(i in c(1:length(start))){	
        b = max*0.035
        signif_threshold = 0.05
        #pval1 = signif(pvals[i], digits = 2)
        pval1 = "NS"
        if(pvals[i]<0.05){pval1 = "*"}
        if(pvals[i]<0.005){pval1 = "**"}
        if(pvals[i]<0.0005){pval1 = "***"}
        y = max(unlist(groups))
        y = y+0.5*b
        segments(start[i],y+b, start[i]+1,y+b,lwd = 3, col = "darkgrey")
        text(mean(c(start[i], start[i]+1)), y+4*b, labels = pval1, cex = cex)
      }
    }
  }
  
  dev.off()
  
  ########################### plot CDR3 lengths between private and shared in blood and biopsy
  analysis=concat(c("Clone_sizes_private_immunosurvielling_", type1))
  fileout1=concat(c(output_directory,"Overlap_",analysis,"_boxplots_", batch,".pdf"))
  w=2.35
  pdf(file=fileout1, height=w*1.2*3, width=w*0.6*6)
  par(mfrow= c(3,5), mar = c(11,5,4,0.1))
  
  chain1_cdr3 = c(list(m_CDR3_length_chain1_sh), list(m_CDR3_length_chain1_p))
  chain2_cdr3 = c(list(m_CDR3_length_chain2_sh), list(m_CDR3_length_chain2_p))
  chain12_cdr3 = c(list(chain1_cdr3), list(chain2_cdr3))
  names(chain12_cdr3) = c("IGH","IGL")
  
  for(s in c(1:length(chain12_cdr3))){
    shared = chain12_cdr3[[s]][[1]]
    private = chain12_cdr3[[s]][[2]]
    cell_type_use = colnames(shared)[intersect(which(apply(shared, 2, function(x){length(which(x!=-1))})>=4), which(apply(private, 2, function(x){length(which(x!=-1))})>=4))]
    
    groups = NULL
    for(c in c(1:length(cell_type_use))){
      x = shared[,cell_type_use[c]]
      y = private[,cell_type_use[c]]
      x = x[which(x!=-1)]
      y = y[which(y!=-1)]
      g1 = c(list(x), list(y))
      names(g1) = c("shared","private")
      groups = c(groups, list(g1))
    }
    names(groups) = cell_type_use
    
    factors = names(groups)
    main = concat(c(names(chain12_cdr3)[s]," CDR3 lengths"))
    max = 22
    min = 10
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = length(factors)
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.15
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.2/max(shift)
    
    library(RColorBrewer)
    cols = add.alpha (brewer.pal(8, "Set1")[c(1,2)], alpha = 0.95)
    #cols = add.alpha (brewer.pal(8, "Set1")[c(1,2,1,2)], alpha = 0.5)
    cols1 = add.alpha (cols,alpha = 0.5)
    #cols2 =  brewer.pal(8, "Set1")[c(1,2,1,2)]
    b = max*0.035
    signif_threshold = 0.05
    pval1 = signif(pvals[i], digits = 2)
    y = max(unlist(groups))
    y = y+1*b
    for(i in c(1:l)){
      pval1 = t.test(groups[[i]][[1]], groups[[i]][[2]], paired = T)$p.value
      l2 = length(groups[[i]][[1]])
      for(i1 in c(1:l2)){
        segments(i-shift[1], groups[[i]][[1]][i1], i-shift[2], groups[[i]][[2]][i1], col = "grey")
      }
      pval2 = "NS"
      if(pval1<0.05){pval2 ="*"}
      segments(i-shift[1],y+b, i-shift[2],y+b,lwd = 3, col = "darkgrey")
      text(i, y+2.5*b, labels = pval2, cex = cex+0.1)
    }
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols1[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21,bg = cols1[i1], col=cols[i1], cex = 0.7)
      }}
  }
  
  
  dev.off()
  
  ############ overlap
  count_per_sample = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
  prop_per_sample = count_per_sample
  w_clone = names(which(clone!='-'))
  for(i in c(1:length(samples))){
    w = intersect(which(sample == samples[i]), which(id %in% w_clone))
    t = table(cell_type[w])
    count_per_sample[i,names(t)] = t
    prop_per_sample[i,names(t)] = t*100/sum(t)
  }
  totals = rowSums(count_per_sample)
  
  threshold = 50
  length(which(totals<threshold)) ### which to not include
  repeats = 1000
  
  overlap_types = c("overlap_cells", "overlap_clones")
  m_overlap_numbers = matrix(data = -1,nrow = length(samples), ncol = length(overlap_types), dimnames = c(list(samples), list(overlap_types)))	
  overlap_clone_counts = rep(-1,length(Patients))
  names(overlap_clone_counts) = Patients
  
  for(s in c(1:length(Patients))){
    w1 = intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))
    w2 = intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))
    w1 = intersect(w1, which(id %in% w_clone))
    w2 = intersect(w2, which(id %in% w_clone))
    shared_clones = intersect(clone[id[w1]] , clone[id[w2]])
    cl1 = clone[id[w1]]
    cl2 = clone[id[w2]]
    cls = c(list(cl1),list(cl2)) ## shared clones
    if(min(c(length(cl1), length(cl2)))> threshold){
      overlap1 = NULL
      overlap2 = NULL
      overlap_clones = NULL
      for(r in c(1:repeats)){
        rand1 = sample(cl1, threshold)
        rand2 = sample(cl2, threshold)
        o = intersect(rand1, rand2)
        overlap_clones = c(overlap_clones, length(o))
        overlap1 = c(overlap1, length(which(rand1 %in% o))*100/length(rand1))
        overlap2 = c(overlap2, length(which(rand2 %in% o))*100/length(rand2))
      }
      m_overlap_numbers[concat(c(Patients[s],"-",Sample.Types[1])), ] = c(mean(overlap1), mean(overlap_clones))
      m_overlap_numbers[concat(c(Patients[s],"-",Sample.Types[2])), ] = c(mean(overlap2), mean(overlap_clones))
      overlap_clone_counts[Patients[s]] = mean(overlap_clones)
    }
  }
  
  type = c(rep("shared", length(samples)), rep("private", length(samples)))
  m_isotypes = rbind(m_isotypes_shared, m_isotypes_private)
  m_isotypes = cbind(type, m_isotypes)
  
  m_v_genes = rbind(m_v_genes_shared, m_v_genes_private)
  m_v_genes = cbind(type, m_v_genes)
  
  m_cell_type = rbind(m_cell_type_shared, m_cell_type_private)
  m_cell_type = cbind(type, m_cell_type)
  
  type = "BCR"
  analysis = "VDJ_overlap_isotype"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(m_isotypes, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "VDJ_overlap_V_genes"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(m_v_genes, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "VDJ_overlap_cell_type"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(m_cell_type, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "VDJ_overlap_rel_proportions"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(m_overlap_numbers, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "VDJ_overlap_rel_proportions_n_clones"
  out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(cbind(overlap_clone_counts), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
}
BCR_immunosurveillance(VDJ_list_BCR, output_directory, batch, groups_PCA)

TCR_immunosurveillance<-function(VDJ_list_TCR, output_directory, batch, groups_PCA){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  analysis = "Clonal_overlap_between_sites"
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[6]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  w_clone = names(which(clone!='-'))
  clone = VDJ[,"clone2"]
  clone1 = VDJ[,"clone1"]
  clone2 = VDJ[,"clone2"]
  cdr3_1 = VDJ[,"cdr3_aa1"]
  cdr3_2 = VDJ[,"cdr3_aa2"]
  names(cdr3_1) = rownames(VDJ)
  names(cdr3_2) = rownames(VDJ)
  names(clone) = rownames(VDJ)
  names(clone1) = rownames(VDJ)
  names(clone2) = rownames(VDJ)
  constant_region1 = VDJ[,"constant_region1"]
  names(constant_region1) = rownames(VDJ)
  
  V_gene_10X1 = VDJ[,"V_gene_10X1"]
  V_mm1 = VDJ[,"V_mm1"]
  names(constant_region1) = rownames(VDJ)
  names(V_gene_10X1) = rownames(VDJ)
  names(V_mm1) = rownames(VDJ)
  constant_region1s = sort(unique(constant_region1))
  V_gene_10X1s = sort(unique(V_gene_10X1))
  cell_ids = rownames(VDJ)
  id = cell_ids
  
  #### get subsampling threshold
  min = 2
  
  #### 
  w_clone = intersect(names(which(clone1!='-')),names(which(clone2!='-')))
  w_clone = names(which(clone2!='-'))
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type1 = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    w_cd48 = names(cell_type)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[w_cd48])
    w_clone = intersect(names(which(clone2!='-')),w_cd48)
    chain = names(CD48)[cd48]
    cell_types1 = cell_types_use
    
    m_v_genes_shared = matrix(data = -1,nrow = length(samples), ncol = length(V_gene_10X1s), dimnames = c(list(samples), list(V_gene_10X1s)))
    m_v_genes_private = matrix(data = -1,nrow = length(samples), ncol = length(V_gene_10X1s), dimnames = c(list(samples), list(V_gene_10X1s)))	
    m_cell_type_shared = matrix(data = -1,nrow = length(samples), ncol = length(cell_types1), dimnames = c(list(samples), list(cell_types1)))
    m_cell_type_private = matrix(data = -1,nrow = length(samples), ncol = length(cell_types1), dimnames = c(list(samples), list(cell_types1)))	
    m_CDR3_length_chain1_p = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types1), dimnames = c(list(Patients), list(cell_types1)))	
    m_CDR3_length_chain1_sh = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types1), dimnames = c(list(Patients), list(cell_types1)))	
    m_CDR3_length_chain2_p = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types1), dimnames = c(list(Patients), list(cell_types1)))	
    m_CDR3_length_chain2_sh = matrix(data = -1,nrow = length(Patients), ncol = length(cell_types1), dimnames = c(list(Patients), list(cell_types1)))	
    shared_private_clone_sizes = NULL
    
    print("Getting overlap")
    for(s in c(1:length(Patients))){
      print(s)
      w1 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))]
      w2 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))]
      w1 = intersect(w1, w_cd48)
      w2 = intersect(w2, w_cd48)
      w1 = intersect(w1, w_clone)
      w2 = intersect(w2, w_clone)
      shared_clones = intersect(clone[w1] , clone[w2])
      cl1 = intersect(names(clone [which(clone %in% shared_clones)]), w1)
      cl2 = intersect(names(clone [which(clone %in% shared_clones)]), w2)
      cls = c(list(cl1),list(cl2)) ## shared clones
      
      pl1 = intersect(names(clone [which(clone %in% shared_clones==F)]), w1)
      pl2 = intersect(names(clone [which(clone %in% shared_clones==F)]), w2)
      pls = c(list(pl1),list(pl2)) ## private clones
      
      t_biopsy = length(w1)
      t_blood  = length(w2)
      shared_clones_biopsy = as.numeric((table(clone[cl1])-1)*100/t_biopsy)
      shared_clones_blood = as.numeric ((table(clone [cl2])-1)*100/t_blood)
      private_clones_biopsy = as.numeric ((table(clone [pl1])-1)*100/t_biopsy )
      private_clones_blood = as.numeric ((table(clone [pl2])-1)*100/t_blood )
      
      # private clones
      cdr31 = nchar( cdr3_1[pls[[1]]])
      cdr32 = nchar( cdr3_2[pls[[1]]])
      cdr31 = cdr31[which(is.na(cdr31)==F)]
      cdr32 = cdr32[which(is.na(cdr32)==F)]
      t1a = table(cell_type[names(cdr31)], cdr31)
      t1b = table(cell_type[names(cdr32)], cdr32)
      
      cdr31 = nchar( cdr3_1[cls[[1]]])
      cdr32 = nchar( cdr3_2[cls[[1]]])
      cdr31 = cdr31[which(is.na(cdr31)==F)]
      cdr32 = cdr32[which(is.na(cdr32)==F)]
      t2a = table(cell_type[names(cdr31)], cdr31)
      t2b = table(cell_type[names(cdr32)], cdr32)
      
      ct_use = names(which(rowSums(t2a)>=5))
      apply(t1a, 1, function(x){ sum(x*as.numeric(colnames(t1a))/sum(x))})
      
      m_CDR3_length_chain1_p[Patients[s], ct_use] = apply(t1a, 1, function(x){ sum(x*as.numeric(colnames(t1a))/sum(x))})[ct_use]
      m_CDR3_length_chain1_sh[Patients[s], ct_use] = apply(t2a, 1, function(x){ sum(x*as.numeric(colnames(t2a))/sum(x))})[ct_use]
      m_CDR3_length_chain2_p[Patients[s], ct_use] = apply(t1b, 1, function(x){ sum(x*as.numeric(colnames(t1b))/sum(x))})[ct_use]
      m_CDR3_length_chain2_sh[Patients[s], ct_use] = apply(t2b, 1, function(x){ sum(x*as.numeric(colnames(t2b))/sum(x))})[ct_use]
      
      g = c(list(private_clones_blood),list(shared_clones_blood),list(private_clones_biopsy),list(shared_clones_biopsy))
      names(g) = c("private blood","immunosurv. blood","private tumour","immunosurv. tumour")
      # boxplot(g)
      # t.test(g[[3]],y=g[[4]])
      shared_private_clone_sizes = c(shared_private_clone_sizes, list(g))
      
      for(i in c(1:length(cls))){
        if(length(cls[[i]])>min){
          sample1 = concat(c(Patients[s],"-",Sample.Types[i]))
          t = table(V_gene_10X1[cls[[i]]])
          m_v_genes_shared[sample1, ]=0
          m_v_genes_shared[sample1, names(t)] = t*100/sum(t)
          t = table(cell_type[cls[[i]]])
          m_cell_type_shared[sample1, ]=0
          m_cell_type_shared[sample1, names(t)] = t*100/sum(t)
        }
        if(length(pls[[i]])>min){
          sample1 = concat(c(Patients[s],"-",Sample.Types[i]))
          t = table(V_gene_10X1[pls[[i]]])
          m_v_genes_private[sample1, ]=0
          m_v_genes_private[sample1, names(t)] = t*100/sum(t)
          t = table(cell_type[pls[[i]]])
          m_cell_type_private[sample1, ]=0
          m_cell_type_private[sample1, names(t)] = t*100/sum(t)
        }
      }
    }
    names(shared_private_clone_sizes) = Patients
    
    type1 = concat(c(type, "_", chain))
    analysis = "cell_type_private"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(cbind(m_cell_type_private), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    analysis = "cell_type_shared"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(cbind(m_cell_type_shared), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    analysis = "V_genes_private"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(cbind(m_v_genes_private), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    analysis = "V_genes_shared"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(cbind(m_v_genes_shared), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    
    ########################### plot clone sizes between private and shared in blood and biopsy
    analysis=concat(c("Clone_sizes_private_immunosurvielling_", type1))
    fileout1=concat(c(output_directory,"Overlap_",analysis,"_boxplots_", batch,".pdf"))
    w=2.35
    pdf(file=fileout1, height=w*1.2*3, width=w*0.6*6)
    par(mfrow= c(3,6), mar = c(11,5,4,0.1))
    
    array1 = NULL
    array2 = NULL
    array3 = NULL
    array4 = NULL
    for(s in c(1:length(Patients))){
      groups = shared_private_clone_sizes[[s]]
      array1 = c(array1, mean(groups[[1]]))
      array2 = c(array2, mean(groups[[2]]))
      array3 = c(array3, mean(groups[[3]]))
      array4 = c(array4, mean(groups[[4]]))
      factors = names(groups)
      main = Patients[s]
      max = max(c(unlist(groups), unlist(groups))*1.25)
      min = 0
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = length(factors)
      max_scale = max
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*5}
      if(range <10){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.3
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.3/max(shift)
      
      library(RColorBrewer)
      cols = add.alpha (brewer.pal(8, "Set1")[c(1,1,2,2)], alpha = 0.5)
      cols[c(2,4)] = add.alpha (cols[c(2,4)], alpha = 0.5)
      #cols = add.alpha (brewer.pal(8, "Set1")[c(1,2,1,2)], alpha = 0.5)
      cols1 = add.alpha (c(cols,alpha = 0.5))
      #cols2 =  brewer.pal(8, "Set1")[c(1,2,1,2)]
      
      for(i in c(1:l)){
        points1=as.numeric(groups[[i]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i,width,cols[i],1, cols1[i])
        points(rep(i, length(points1)),points1, pch =21, col=cols[i], cex = 0.7)
      }
      if(min(c(length(groups[[1]]),length(groups[[2]]),length(groups[[3]]),length(groups[[4]])))>3){
        pval1 = t.test(groups[[1]], groups[[2]])$p.value
        pval2 = t.test(groups[[3]], groups[[4]])$p.value
        pvals = c(pval1, pval2)
        start = c(1,3)
        for(i in c(1:length(start))){	
          b = max*0.035
          signif_threshold = 0.05
          pval1 = signif(pvals[i], digits = 2)
          y = max(unlist(groups[c(start[i], start[i]+1)]))
          y = y+1*b
          segments(start[i],y+b, start[i]+1,y+b,lwd = 3, col = "darkgrey")
          text(mean(c(start[i], start[i]+1)), y+2.5*b, labels = pval1, cex = cex)
        }
      }
    }
    
    mean_shared_private_clone_sizes = cbind(array1, array2, array3, array4)
    colnames(mean_shared_private_clone_sizes) = names(shared_private_clone_sizes[[1]])
    rownames(mean_shared_private_clone_sizes) = Patients
    
    groups_PCA1 = c(groups_PCA, list(unlist(groups_PCA)))
    names(groups_PCA1) = c("ME", "AE", "All")
    for(s in c(1:length(groups_PCA1))){
      mat = mean_shared_private_clone_sizes[gsub("_biopsy","",groups_PCA1[[s]]),]
      groups = NULL
      for(i in c(1:length(mat[1,]))){
        x = mat[,i]
        x = x[which(is.na(x)==F)]
        groups = c(groups, list(x))
      }
      factors = colnames(mat)
      main = names(groups_PCA1)[s]
      if(s==4){main = "all"}
      max = max(mean_shared_private_clone_sizes[which(is.na(mean_shared_private_clone_sizes)==F)])*1.25
      min = 0
      b = (max-min)*0.034
      ylab = "mean clone size (%)"
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = length(factors)
      max_scale = max
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*5}
      if(range <10){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.3
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.3/max(shift)
      
      library(RColorBrewer)
      cols = add.alpha (brewer.pal(8, "Set1")[c(1,1,2,2)], alpha = 0.5)
      cols[c(2,4)] = add.alpha (cols[c(2,4)], alpha = 0.5)
      #cols = add.alpha (brewer.pal(8, "Set1")[c(1,2,1,2)], alpha = 0.5)
      cols1 = add.alpha (c(cols,alpha = 0.5))
      #cols2 =  brewer.pal(8, "Set1")[c(1,2,1,2)]
      
      for(i in c(1:l)){
        points1=as.numeric(groups[[i]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i,width,cols[i],1, cols1[i])
        points(rep(i, length(points1)),points1, pch =21, col=cols[i], cex = 0.7)
      }
      if(min(c(length(groups[[1]]),length(groups[[2]]),length(groups[[3]]),length(groups[[4]])))>=3){
        pval1 = wilcox.test(groups[[1]],groups[[2]],paired = F)$p.value
        pval2 = wilcox.test(groups[[3]], groups[[4]],paired = F)$p.value
        pvals = c(pval1, pval2)
        start = c(1,3)
        for(i in c(1:length(start))){	
          b = max*0.035
          signif_threshold = 0.05
          #pval1 = signif(pvals[i], digits = 2)
          pval1 = "NS"
          if(pvals[i]<0.05){pval1 = "*"}
          if(pvals[i]<0.005){pval1 = "**"}
          if(pvals[i]<0.0005){pval1 = "***"}
          y = max(unlist(groups))
          y = y+0.5*b
          segments(start[i],y+b, start[i]+1,y+b,lwd = 3, col = "darkgrey")
          text(mean(c(start[i], start[i]+1)), y+4*b, labels = pval1, cex = cex)
        }
      }
    }
    
    dev.off()
    
    ########################### plot CDR3 lengths between private and shared in blood and biopsy
    analysis=concat(c("Clone_CDR3_lengths_private_immunosurvielling_", type1))
    fileout1=concat(c(output_directory,"Overlap_",analysis,"_boxplots_", batch,".pdf"))
    w=2.35
    pdf(file=fileout1, height=w*1.2*3, width=w*0.6*6)
    par(mfrow= c(3,5), mar = c(11,5,4,0.1))
    
    chain1_cdr3 = c(list(m_CDR3_length_chain1_sh), list(m_CDR3_length_chain1_p))
    chain2_cdr3 = c(list(m_CDR3_length_chain2_sh), list(m_CDR3_length_chain2_p))
    chain12_cdr3 = c(list(chain1_cdr3), list(chain2_cdr3))
    names(chain12_cdr3) = c("TRA","TRB")
    
    for(s in c(1:length(chain12_cdr3))){
      shared = chain12_cdr3[[s]][[1]]
      private = chain12_cdr3[[s]][[2]]
      cell_type_use = colnames(shared)[intersect(which(apply(shared, 2, function(x){length(which(x!=-1))})>=9), which(apply(private, 2, function(x){length(which(x!=-1))})>=9))]
      
      groups = NULL
      for(c in c(1:length(cell_type_use))){
        x = shared[,cell_type_use[c]]
        y = private[,cell_type_use[c]]
        x = x[which(x!=-1)]
        y = y[which(y!=-1)]
        g1 = c(list(x), list(y))
        names(g1) = c("shared","private")
        groups = c(groups, list(g1))
      }
      names(groups) = cell_type_use
      
      factors = names(groups)
      main = concat(c(names(chain12_cdr3)[s]," CDR3 lengths"))
      max = 17
      min = 8
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = length(factors)
      max_scale = max
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*5}
      if(range <10){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.15
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.2/max(shift)
      
      library(RColorBrewer)
      cols = add.alpha (brewer.pal(8, "Set1")[c(1,2)], alpha = 0.95)
      #cols = add.alpha (brewer.pal(8, "Set1")[c(1,2,1,2)], alpha = 0.5)
      cols1 = add.alpha (cols,alpha = 0.5)
      #cols2 =  brewer.pal(8, "Set1")[c(1,2,1,2)]
      b = max*0.025
      signif_threshold = 0.05
      pval1 = signif(pvals[i], digits = 2)
      y = max(unlist(groups))
      y = y+1*b
      for(i in c(1:l)){
        pval1 = t.test(groups[[i]][[1]], groups[[i]][[2]], paired = T)$p.value
        l2 = length(groups[[i]][[1]])
        for(i1 in c(1:l2)){
          segments(i-shift[1], groups[[i]][[1]][i1], i-shift[2], groups[[i]][[2]][i1], col = "grey")
        }
        pval2 = "NS"
        if(pval1<0.05){pval2 ="*"}
        segments(i-shift[1],y+b, i-shift[2],y+b,lwd = 3, col = "darkgrey")
        text(i, y+2.5*b, labels = pval2, cex = cex+0.1)
      }
      
      for(i in c(1:l)){
        for(i1 in c(1:l1)){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols1[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21,bg = cols1[i1], col=cols[i1], cex = 0.7)
        }}
    }
    
    
    dev.off()
    
    
    ############ overlap
    count_per_sample = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
    prop_per_sample = count_per_sample
    w_clone = names(which(clone!='-'))
    for(i in c(1:length(samples))){
      w = intersect(which(sample == samples[i]), which(id %in% w_clone))
      t = table(cell_type[w])
      count_per_sample[i,names(t)] = t
      prop_per_sample[i,names(t)] = t*100/sum(t)
    }
    totals = rowSums(count_per_sample)
    
    threshold = 50
    length(which(totals<threshold)) ### which to not include
    repeats = 1000
    
    overlap_types = c("overlap_cells", "overlap_clones")
    m_overlap_numbers = matrix(data = -1,nrow = length(samples), ncol = length(overlap_types), dimnames = c(list(samples), list(overlap_types)))	
    overlap_clone_counts = rep(-1,length(Patients))
    names(overlap_clone_counts) = Patients
    
    for(s in c(1:length(Patients))){
      w1 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))]
      w2 = id[intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))]
      w1 = intersect(w1, w_cd48)
      w2 = intersect(w2, w_cd48)
      w1 = intersect(w1, w_clone)
      w2 = intersect(w2, w_clone)
      shared_clones = intersect(clone[w1] , clone[w2])
      cl1 = clone[w1]
      cl2 = clone[w2]
      cls = c(list(cl1),list(cl2)) ## shared clones
      if(min(c(length(cl1), length(cl2)))> threshold){
        overlap1 = NULL
        overlap2 = NULL
        overlap_clones = NULL
        for(r in c(1:repeats)){
          rand1 = sample(cl1, threshold)
          rand2 = sample(cl2, threshold)
          o = intersect(rand1, rand2)
          overlap_clones = c(overlap_clones, length(o))
          overlap1 = c(overlap1, length(which(rand1 %in% o))*100/length(rand1))
          overlap2 = c(overlap2, length(which(rand2 %in% o))*100/length(rand2))
        }
        m_overlap_numbers[concat(c(Patients[s],"-",Sample.Types[1])), ] = c(mean(overlap1), mean(overlap_clones))
        m_overlap_numbers[concat(c(Patients[s],"-",Sample.Types[2])), ] = c(mean(overlap2), mean(overlap_clones))
        overlap_clone_counts[Patients[s]] = mean(overlap_clones)
      }
    }
    
    type = c(rep("shared", length(samples)), rep("private", length(samples)))
    
    m_v_genes = rbind(m_v_genes_shared, m_v_genes_private)
    m_v_genes = cbind(type, m_v_genes)
    
    m_cell_type = rbind(m_cell_type_shared, m_cell_type_private)
    m_cell_type = cbind(type, m_cell_type)
    
    
    type = "TCR"
    
    analysis = "VDJ_overlap_V_genes"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(m_v_genes, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    analysis = "VDJ_overlap_cell_type"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(m_cell_type, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    analysis = "VDJ_overlap_rel_proportions"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(m_overlap_numbers, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    analysis = "VDJ_overlap_rel_proportions_n_clones"
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(cbind(overlap_clone_counts), file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  }
  
}
TCR_immunosurveillance(VDJ_list_TCR, output_directory, batch, groups_PCA)


Get_proportions_of_immunosurveilling_cell_types_B_cells<-function(VDJ_list_BCR, output_directory, batch, groups_PCA){
  VDJ_list = VDJ_list_BCR
  type = "BCR"
  type1 = "B_cells"
  analysis = "Clonal_overlap_between_sites"
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[6]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  w_clone = names(which(clone!='-'))
  clone1 = VDJ[,"clone1"]
  clone = VDJ[,"clone1"]
  clone2 = apply(cbind("L:",VDJ[,"clone2"]),1,paste, collapse = "")
  clone2[which(clone2=="L:-")] = "-"
  
  cdr3_1 = VDJ[,"cdr3_aa1"]
  cdr3_2 = VDJ[,"cdr3_aa2"]
  names(cdr3_1) = rownames(VDJ)
  names(cdr3_2) = rownames(VDJ)
  names(clone) = rownames(VDJ)
  names(clone1) = rownames(VDJ)
  names(clone2) = rownames(VDJ)
  constant_region1 = VDJ[,"constant_region1"]
  names(constant_region1) = rownames(VDJ)
  
  V_gene_10X1 = VDJ[,"V_gene_10X1"]
  V_mm1 = VDJ[,"V_mm1"]
  names(constant_region1) = rownames(VDJ)
  names(V_gene_10X1) = rownames(VDJ)
  names(V_mm1) = rownames(VDJ)
  constant_region1s = sort(unique(constant_region1))
  V_gene_10X1s = sort(unique(V_gene_10X1))
  cell_ids = rownames(VDJ)
  id = cell_ids
  
  # consider both heavy and light chains
  clone_use = clone2
  w = which(clone_use=="-")
  clone_use[w] = clone1[w]
  
  length(which(clone_use!="-"))
  
  #### identify shared/private clones check by heavy chain
  blood_ids = id[which(Sample.Type=="blood")]
  tumour_ids = id[which(Sample.Type=="biopsy")]
  inter = intersect(names(clone_use), blood_ids)
  blood_clones = unique(clone_use[inter])
  inter = intersect(names(clone_use), tumour_ids)
  tumour_clones = unique(clone_use[inter])
  shared = intersect(blood_clones, tumour_clones)
  shared = shared[which(shared!='-')]
  shared_tumour_ids = intersect(names(clone_use[which(clone_use %in% shared)]), tumour_ids)
  
  #### get tumour cell type proportions
  shared_tumour_ids = unique(c(shared_tumour_ids, shared_tumour_ids))
  inter = intersect(names(clone), tumour_ids) ## only include cells with VDJ
  private_tumour_ids = setdiff(inter, shared_tumour_ids)
  
  # check
  length(shared_tumour_ids) + length(private_tumour_ids)
  length(intersect(tumour_ids, names(clone)))
  
  t1 =  t(table( cell_type[shared_tumour_ids], sample[shared_tumour_ids]))
  t2 =  t(table( cell_type[private_tumour_ids], sample[private_tumour_ids]))
  t1b = t2*0
  t1b[rownames(t1), colnames(t1)] = t1
  t1 = t1b
  rowSums(t1)
  # normalise
  for(i in c(1:length(rownames(t1)))){
    if(sum(t1[i,])>=3){t1[i,] = t1[i,]*100/sum(t1[i,])
    }else{t1[i,] = -1}
    if(sum(t2[i,])>=3){t2[i,] = t2[i,]*100/sum(t2[i,])
    }else{t2[i,] = -1}
  }
  
  analysisa = "cell_type_prop"
  ct = "_private"
  file_table1 = concat(c(output_directory,"Overlap_", analysisa, ct,"_information_",type1,"_PDAC150Ka.txt"))
  write.table(t2, file = file_table1, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ct = "_shared"
  file_table2 = concat(c(output_directory,"Overlap_",  analysisa, ct,"_information_",type1,"_PDAC150Ka.txt"))
  write.table(t1, file = file_table2, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  # V gene usages between shared and non-shared? plot by PCA? 
  inter = intersect(c(names(clone),names(clone2)), tumour_ids)
  t0 =  t(table( V_gene_10X1[inter], sample[inter]))*0
  t1 =  t(table( V_gene_10X1[shared_tumour_ids], sample[shared_tumour_ids]))
  t2 =  t(table( V_gene_10X1[private_tumour_ids], sample[private_tumour_ids]))
  t1b = t0
  t1b[rownames(t1), colnames(t1)] = t1
  t1 = t1b
  t2b = t0
  t2b[rownames(t2), colnames(t2)] = t2
  t2 = t2b
  
  rowSums(t1)
  # normalise
  for(i in c(1:length(rownames(t1)))){
    if(sum(t1[i,])>=3){t1[i,] = t1[i,]*100/sum(t1[i,])
    }else{t1[i,] = -1}
    if(sum(t2[i,])>=3){t2[i,] = t2[i,]*100/sum(t2[i,])
    }else{t2[i,] = -1}
  }
  
  analysisa = "V_gene_prop"
  ct = "_private"
  file_table1 = concat(c(output_directory,"Overlap_",  analysisa, ct,"_information_",type1,"_PDAC150Ka.txt"))
  write.table(t2, file = file_table1, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ct = "_shared"
  file_table2 = concat(c(output_directory,"Overlap_",  analysisa, ct,"_information_",type1,"_PDAC150Ka.txt"))
  write.table(t1, file = file_table2, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ######## get which cell type are connected between blood and tumour
  isotypes = sort(unique(constant_region1))
  #isotypes = isotypes[which(isotypes!="-")]
  ct1 = NULL
  ct2 = NULL
  for(i in c(1:length(cell_types))){
    for(j in c(1:length(cell_types))){
      ct1 = c(ct1, cell_types[i])
      ct2 = c(ct2, cell_types[j])
    }
  }
  ct12 = apply(cbind(ct1,ct2), 1, paste, collapse= ' - ')
  w_clone = names(which(clone_use!='-'))
  
  iso1 = NULL
  iso2 = NULL
  for(i in c(1:length(isotypes))){
    for(j in c(1:length(isotypes))){
      iso1 = c(iso1, isotypes[i])
      iso2 = c(iso2, isotypes[j])
    }
  }
  iso12 = apply(cbind(iso1,iso2), 1, paste, collapse= ' : ')
  
  
  
  m_pat_shared_freq = t(matrix(data = 0,nrow = length(Patients), ncol = length(ct12), dimnames = c(list(Patients), list(ct12))))
  m_pat_shared_freq_isotype = t(matrix(data = 0,nrow = length(Patients), ncol = length(iso12), dimnames = c(list(Patients), list(iso12))))
  for(s in c(1:length(Patients))){
    w1 = intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))
    w2 = intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))
    w1 = intersect(w1, which(id %in% w_clone))
    w2 = intersect(w2, which(id %in% w_clone))
    clone1 = clone_use[id[w1]]
    clone2 = clone_use[id[w2]]
    clone1 = clone1[which(names(clone1) %in% names(cell_type))]
    clone2 = clone2[which(names(clone2) %in% names(cell_type))]
    ####### absolute counts
    shared_clones = intersect(clone1, clone2)
    if(length(shared_clones)>0){
      for(c in c(1:length(shared_clones))){
        x1 = sort(unique(cell_type [names(clone1[which(clone1 %in% shared_clones[c])])]))
        x2 = sort(unique(cell_type [names(clone2[which(clone2 %in% shared_clones[c])])]))
        for(c1 in c(1:length(x1))){
          for(c2 in c(1:length(x2))){
            nam = concat(c(x1[c1]," - ",x2[c2]))
            m_pat_shared_freq[nam,Patients[s]] = m_pat_shared_freq[nam,Patients[s]] +1
          }}
        
        x1 = sort(unique(constant_region1 [names(clone1[which(clone1 %in% shared_clones[c])])]))
        x2 = sort(unique(constant_region1 [names(clone2[which(clone2 %in% shared_clones[c])])]))
        for(c1 in c(1:length(x1))){
          for(c2 in c(1:length(x2))){
            nam = concat(c(x1[c1]," : ",x2[c2]))
            m_pat_shared_freq_isotype[nam,Patients[s]] = m_pat_shared_freq_isotype[nam,Patients[s]] +1
          }}
      }
    }
    m_pat_shared_freq[,Patients[s]] = m_pat_shared_freq[,Patients[s]]*100/length(w2)
  }
  m_pat_shared_freq_isotype = m_pat_shared_freq_isotype[grep("-",rownames(m_pat_shared_freq_isotype),invert = T),]
  
  
  analysis = "BCR_cell_type_biopsy_blood_overlap"
  out_file_table = concat(c(output_directory,"Overlap_", analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(m_pat_shared_freq, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "BCR_isotype_biopsy_blood_overlap"
  out_file_table = concat(c(output_directory,"Overlap_", analysis,"_information_",type,"_PDAC150Ka.txt"))
  write.table(m_pat_shared_freq_isotype, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  analysis = "BCR_cell_type_biopsy_blood_overlap"
  file = concat(c(output_directory,"Overlap_", analysis,"_information_",type,"_PDAC150Ka.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  m_pat_shared_freq = p
  repeats = 10000
  
  m_pat_shared_freq1 = m_pat_shared_freq[which(rowSums(m_pat_shared_freq)>0),]
  for(i in c(1:length(m_pat_shared_freq1[1,]))){m_pat_shared_freq1[,i] = m_pat_shared_freq1[,i]*100/sum(m_pat_shared_freq1[,i])}
  
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  mat_stat = t(m_pat_shared_freq1)
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("-biopsy","",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  mat1 = mat_stat[gsub("_biopsy","",group_PCA_list),]
  wx = which(apply(mat_stat, 2, function(x){length(which(x!=0))})>=3)
  wx = which(is.na(rowSums(mat_stat))==F)
  
  #for(i in c(1:length(mat1[,1]))){mat1[i,] = mat1[i,]+runif(length(mat1[i,]))*1e-10}
  
  factor = factor(factor_PCA)
  fit = manova(formula = mat1[wx,]^1 ~ factor[wx])
  p1 = summary.aov(fit)
  nam = gsub(" Response ","",names(p1))
  p_value = NULL
  means = NULL
  
  i1 = 0
  for(i in p1){
    i1 = i1+1
    p_value = c(p_value, i$'Pr(>F)'[1]) 
    if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
    }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
  }
  p_value[which(is.na(p_value))] = 2
  print(sort(p_value))
  
  colnames(means) = paste("mean.group.", c(1:length(means[1,])))
  combined_p_value = cbind(p_value ,means)
  rownames(combined_p_value) = nam
  p.site = rep("biopsy", length(nam))
  p.analysis = rep(analysis, length(nam))
  x = cbind(p.site, p.analysis, combined_p_value)
  summary_tables = x
  
  table(apply(mat1, 2, function(x){length(which(x!=0))}))
  
  fileout1=concat(c(output_directory,"Overlap_", analysis,"_PCA_groups_boxplots_", batch,".pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1.5*1, width=w*1*2.5)
  par(mfrow= c(1,1), mar = c(16,5,3,3))
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)
  
  mat1 = mat_stat
  wx = which(apply(mat_stat, 2, function(x){length(which(x!=0))})>=3)
  mat1 = mat1[,wx]
  groups = NULL
  p_values = combined_p_value[,1]
  for(i in c(1:length(mat1[1,]))){
    g1 = NULL
    data = NULL
    factor = NULL
    for(g in c(1:length(groups_PCA))){
      x = mat1 [gsub("_biopsy","",groups_PCA[[g]]), i]
      x = x[which(x!=-1)]
      data = c(data,x )
      factor = c(factor, rep(g, length(x)))
      g1 = c(g1, list(x))}
    groups = c(groups, list(g1))
  }
  factors1 = paste("group",c(1:length(groups_PCA)))
  factors = gsub("_"," ",colnames(mat1))
  factors = gsub("activated","act.", factors)
  factors = gsub("memory","mem.", factors)
  factors = gsub("Ag presenting","Ag pres.", factors)
  factors = gsub("SHM diversifying","SHM div.", factors)
  
  main = concat(c("Relative clone sharing blood biopsy"))
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = 0
  b = (max-min)*0.034
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = max
  range = max-min
  if(range>50){scale = c(0:100)*20}
  if(range<=50){scale = c(0:100)*5}
  if(range <10){scale = c(0:100)*2.5}
  if(range <5){scale = c(0:100)*1}
  if(range <4){scale = c(0:100)*0.5}
  if(range <1.5){scale = c(0:1000)*0.2}
  if(range <0.5){scale = c(0:100)*0.1}
  if(range <0.1){scale = c(0:100)*0.01}
  if(range <0.01){scale = c(0:100)*0.001}
  if(range <0.0002){scale = c(0:100)*0.00002}
  cex = 0.9
  Fun<-function(x){x}
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 0.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.25/max(shift)
  
  for(i in c(1:l)){
    for(i1 in c(1:l1)){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.7)
    }}
  
  for(i1 in c(1:length(p_values))){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i1]<signif_threshold){
      i = which(colnames(mat1)==names(p_values)[i1])
      pval1 = "*"
      if(p_values[i] <signif_threshold/10){pval1 = "**"}
      if(p_values[i] <signif_threshold/100){pval1 = "***"}
      y = max(unlist(groups[[i]]))
      y = y+1*b
      # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
      text(i, y+2*b, labels = pval1, cex = 1.3)
    }
  }
  
  dev.off()
  
  ################################################################
  #overall clonal sharing: 
  m_pat_shared_freq1 = m_pat_shared_freq
  for(i in c(1:length(m_pat_shared_freq1[1,]))){m_pat_shared_freq1[,i] = m_pat_shared_freq1[,i]*100/sum(m_pat_shared_freq1[,i])}
  wx = which(is.na(colSums(m_pat_shared_freq1))==F)
  m_pat_shared_freq1 = m_pat_shared_freq1[,wx]
  mat_stat = t(m_pat_shared_freq1)
  width = 1
  height = -1*length(cell_types)
  
  fileout1=concat(c(output_directory,"Overlap_",analysis,"_network_", batch,".pdf"))
  w=3.35
  pdf(file=fileout1, height=w*0.75, width=w*1.5)
  par(mfrow= c(1,1), mar = c(2,10,2,10))
  
  plot(c(0, width),c(-0.5, height-0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "clone sharing between blood and pancreatic tumour", axes = F)
  
  mtext(side = 2, text = cell_types, line = 0.1,cex= cex-0.1, at = -1*c(1:length(cell_types)), las = 1, font = 1)
  mtext(side = 4, text = cell_types, line = 0.1,cex= cex-0.1, at = -1*c(1:length(cell_types)), las = 1, font = 1)
  
  mtext(side = 2, text = gsub("  "," ",Sample.Types[1],fixed = T), line = 0.1,cex= cex-0.1, at = 0, las = 1, font = 2)
  mtext(side = 4, text = gsub("  "," ",Sample.Types[2],fixed = T), line = 0.1,cex= cex-0.1, at = 0, las = 1, font = 2)
  
  edge_strength = apply(mat_stat,2,function(x){mean(x[which(is.na(x)==F)])})
  edge_strength = edge_strength
  edge_strength = edge_strength*1#/max(edge_strength)
  col = add.alpha("darkblue", alpha = 0.5)
  p_val_sub = p_values
  names(p_val_sub)==colnames(mat1)
  cols = rep(add.alpha("black", alpha = 0.8), length(edge_strength))
  names(cols) = names(p_val_sub)
  
  split = strsplit(colnames(mat_stat)," - ",fixed = T)
  f1 = NULL
  f2 = NULL
  
  for(i in c(1:length(split))){
    from = which(cell_types== split[[i]][1])
    to = which(cell_types== split[[i]][2])
    
    if(length(c(from,to))==2){
      f1 = c(f1, from)
      f2 = c(f2, to)
      #segments(0,-1*from, x1 = width, y1 = -1*to, lwd = 1, col = "black")}
      segments(0,-1*from, x1 = width, y1 = -1*to, lwd = edge_strength[i], col = cols[i])
    }
  }
  
  dev.off()
  
  ######### isotype overlap
  analysis = "BCR_isotype_biopsy_blood_overlap"
  file = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  m_pat_shared_freq = p
  repeats = 10000
  
  m_pat_shared_freq1 = m_pat_shared_freq[which(rowSums(m_pat_shared_freq)>0),]
  for(i in c(1:length(m_pat_shared_freq1[1,]))){m_pat_shared_freq1[,i] = m_pat_shared_freq1[,i]*100/sum(m_pat_shared_freq1[,i])}
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  mat_stat = t(m_pat_shared_freq1)
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("-biopsy","",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  mat_stat = mat_stat[gsub("_biopsy","",group_PCA_list),]
  
  factor = factor(factor_PCA)
  fit = manova(formula = mat_stat ~ factor)
  p1 = summary.aov(fit)
  nam = gsub(" Response ","",names(p1))
  p_value = NULL
  means = NULL
  
  i1 = 0
  for(i in p1){
    i1 = i1+1
    p_value = c(p_value, i$'Pr(>F)'[1]) 
    if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
    }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
  }
  p_value[which(is.na(p_value))] = 2
  print(sort(p_value))
  
  colnames(means) = paste("mean.group.", c(1:length(means[1,])))
  combined_p_value = cbind(p_value ,means)
  rownames(combined_p_value) = nam
  p.site = rep("biopsy", length(nam))
  p.analysis = rep(analysis, length(nam))
  x = cbind(p.site, p.analysis, combined_p_value)
  summary_tables = x
  
  table(apply(mat1, 2, function(x){length(which(x!=0))}))
  
  fileout1=concat(c(output_directory,"Overlap_",analysis,"_PCA_groups_boxplots_", batch,".pdf"))
  w=3.25
  pdf(file=fileout1, height=w*1.4*1, width=w*1*1.5)
  par(mfrow= c(1,1), mar = c(16,5,3,3))
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)
  
  mat1 = mat_stat
  groups = NULL
  p_values = combined_p_value[,1]
  for(i in c(1:length(mat1[1,]))){
    g1 = NULL
    data = NULL
    factor = NULL
    for(g in c(1:length(groups_PCA))){
      x = mat1 [gsub("_biopsy","",groups_PCA[[g]]), i]
      x = x[which(x!=-1)]
      data = c(data,x )
      factor = c(factor, rep(g, length(x)))
      g1 = c(g1, list(x))}
    groups = c(groups, list(g1))
  }
  factors1 = paste("group",c(1:length(groups_PCA)))
  factors = gsub("_"," ",colnames(mat1))
  factors = gsub("activated","act.", factors)
  factors = gsub("memory","mem.", factors)
  factors = gsub("Ag presenting","Ag pres.", factors)
  factors = gsub("SHM diversifying","SHM div.", factors)
  
  main = concat(c("Relative clone sharing blood biopsy"))
  max = max(c(unlist(groups), unlist(groups))*1.15)
  min = 0
  b = (max-min)*0.034
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = max
  range = max-min
  if(range>50){scale = c(0:100)*20}
  if(range<=50){scale = c(0:100)*5}
  if(range <10){scale = c(0:100)*2.5}
  if(range <5){scale = c(0:100)*1}
  if(range <4){scale = c(0:100)*0.5}
  if(range <1.5){scale = c(0:1000)*0.2}
  if(range <0.5){scale = c(0:100)*0.1}
  if(range <0.1){scale = c(0:100)*0.01}
  if(range <0.01){scale = c(0:100)*0.001}
  if(range <0.0002){scale = c(0:100)*0.00002}
  cex = 0.9
  Fun<-function(x){x}
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 0.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.25/max(shift)
  
  for(i in c(1:l)){
    for(i1 in c(1:l1)){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.7)
    }}
  
  for(i1 in c(1:length(p_values))){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i1]<signif_threshold){
      i = which(colnames(mat1)==names(p_values)[i1])
      pval1 = "*"
      if(p_values[i] <signif_threshold/10){pval1 = "**"}
      if(p_values[i] <signif_threshold/100){pval1 = "***"}
      y = max(unlist(groups[[i]]))
      y = y+1*b
      # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
      text(i, y+2*b, labels = pval1, cex = 1.3)
    }
  }
  
  dev.off()
  
  
  
}
Get_proportions_of_immunosurveilling_cell_types_B_cells(VDJ_list_BCR, output_directory, batch, groups_PCA)

Get_proportions_of_immunosurveilling_cell_types_T_cells<-function(VDJ_list_TCR, output_directory, batch, groups_PCA){
  types = c("TCR_CD4","TCR_CD8")
  names(types) = c("CD4","CD8")
  
  type = "TCR"
  file = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Data and references/VDJ_information_",type,"_PDAC150Ka.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  table(as.numeric(p[,"n_umis2"]))
  w = which(p[,"n_umis2"]=="1") ## exclude nUMIs ==1
  p[w,] = '-'
  
  clone1 = p[,"clone1"]
  names(clone1) = rownames(p)
  length(which(clone1!='-'))
  clone2 = p[,"clone2"]
  names(clone2) = rownames(p)
  length(which(clone2!='-'))
  V_gene_10X1 = p[,"V_gene2"]
  names(V_gene_10X1) = rownames(p)
  V_gene_10X1s = sort(unique(V_gene_10X1))
  clone = apply(cbind(clone1, clone2), 1, paste, collapse = "-")
  clone = clone2
  
  file = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Overall UMAPs and annotations/Cell_annotation_ALL_PDAC150Ka.txt"
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  p1 <- p1[grep("T cell", p1[,"cell_type"]), ]
  id1 = p1[,"barcode"]
  id2 = p1[,"barcode"]
  sample = p1[,"sample"]
  origin_sample = sort(unique(sample))
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
  names(sample)= id
  Patients = sort(unique(Patient))
  Sample.Types = rev(sort(unique(Sample.Type)))
  cell_types = sort(unique(cell_type))
  
  cd4cd8 = cell_type
  names(cd4cd8) = id
  cd4cd8[grep("CD4", cd4cd8)]= "CD4"
  cd4cd8[grep("CD8", cd4cd8)]= "CD8"
  cd4cd8[grep("NK like", cd4cd8)]= "CD8"
  cd4cd8[grep("MAIT", cd4cd8)]= "CD8"
  cd4cd8[which(cd4cd8 %in% c("CD4","CD8")==F)] = NA
  table(cd4cd8[names(clone)])
  cd4cd8s = sort(unique(cd4cd8))
  
  for(cd in c(1:length(types))){
    type1 = types[cd]
    w_cd = which(cd4cd8==names(types)[cd])
    
    clone_use = clone
    w_clone = which(clone!='-')
    #### identify shared/private clones
    blood_ids = id[intersect(which(Sample.Type=="blood"),w_cd)]
    tumour_ids = id[intersect(which(Sample.Type=="biopsy"),w_cd)]
    inter = intersect(names(clone), blood_ids)
    blood_clones = unique(clone[inter])
    inter = intersect(names(clone), tumour_ids)
    tumour_clones = unique(clone[inter])
    shared = intersect(blood_clones, tumour_clones)
    private = setdiff(tumour_clones,blood_clones)
    private_blood = setdiff(blood_clones,tumour_clones)
    shared = shared[which(shared!='-')]
    #### get tumour cell type proportions
    shared_tumour_ids = intersect(names(clone[which(clone %in% shared)]), tumour_ids)
    private_tumour_ids = intersect(names(clone[which(clone %in% private)]), tumour_ids)
    
    shared_blood_ids = intersect(names(clone[which(clone %in% shared)]), blood_ids)
    private_blood_ids = intersect(names(clone[which(clone %in% private_blood)]), blood_ids)
    
    # check
    length(shared_tumour_ids) + length(private_tumour_ids)
    length(intersect(tumour_ids, names(clone)))
    
    t1 =  t(table( cell_type[shared_tumour_ids], sample[shared_tumour_ids]))
    t2 =  t(table( cell_type[private_tumour_ids], sample[private_tumour_ids]))
    
    t3 =  t(table( cell_type[shared_blood_ids], sample[shared_blood_ids]))
    t4 =  t(table( cell_type[private_blood_ids], sample[private_blood_ids]))
    
    t1dt1pt2 = t1*100/(t1+t2)
    total = t1+t2
    
    t1dt1pt2_blood = t3*100/(t3+t4)
    total_blood = t3+t4
    
    t1b = t2*0
    t1b[rownames(t1), colnames(t1)] = t1
    t1 = t1b
    rowSums(t1)
    # normalise
    for(i in c(1:length(rownames(t1)))){
      if(sum(t1[i,])>=3){t1[i,] = t1[i,]*100/sum(t1[i,])
      }else{t1[i,] = -1}
      if(sum(t2[i,])>=3){t2[i,] = t2[i,]*100/sum(t2[i,])
      }else{t2[i,] = -1}
    }
    
    list_proportion_shared= c(list(t1dt1pt2),list(total),list(t1dt1pt2_blood), list(total_blood ))
    names(list_proportion_shared) = c("prop tumour","total tumour", "prop blood","total blood")
    saveRDS(file = concat(c(output_directory,"Overlap_Proportion_per_cell_type_shared_information_",type1,"_PDAC150Ka.txt")),list_proportion_shared)
    
    analysisa = "cell_type_prop"
    ct = "_private"
    file_table1 = concat(c(output_directory,"Overlap_",analysisa, ct,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(t2, file = file_table1, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    ct = "_shared"
    file_table2 = concat(c(output_directory,"Overlap_", analysisa, ct,"_information_",type1,"_PDAC150Ka.txt"))
    write.table(t1, file = file_table2, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    ### clone sizes
    clone_sizes_s =table(clone[shared_tumour_ids])
    clone_sizes_p =table(clone[private_tumour_ids])
    groups = c(list(clone_sizes_s), list(clone_sizes_p))
    boxplot(groups)
    t.test(groups[[1]]^0.5, y = groups[[2]]^0.5)
    
    
    ######## get which cell type are connected between blood and tumour
    cell_types = sort(unique(cell_type[w_cd]))
    ct1 = NULL
    ct2 = NULL
    for(i in c(1:length(cell_types))){
      for(j in c(1:length(cell_types))){
        ct1 = c(ct1, cell_types[i])
        ct2 = c(ct2, cell_types[j])
      }
    }
    ct12 = apply(cbind(ct1,ct2), 1, paste, collapse= ' - ')
    
    m_pat_shared_freq = t(matrix(data = 0,nrow = length(Patients), ncol = length(ct12), dimnames = c(list(Patients), list(ct12))))
    for(s in c(1:length(Patients))){
      print(s)
      w1 = intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[1]))
      w2 = intersect(which(Patient == Patients[s]), which(Sample.Type== Sample.Types[2]))
      w1 = intersect(w1, which(id %in% names(w_clone)))
      w2 = intersect(w2, which(id %in% names(w_clone)))
      w1 = intersect(w1, which(id %in% names(w_cd)))
      w2 = intersect(w2, which(id %in% names(w_cd)))
      clone1 = clone_use[id[w1]]
      clone2 = clone_use[id[w2]]
      ####### absolute counts
      shared_clones = intersect(clone1, clone2)
      if(length(shared_clones)>0){
        for(c in c(1:length(shared_clones))){
          x1 = sort(unique(cell_type [names(clone1[which(clone1 %in% shared_clones[c])])]))
          x2 = sort(unique(cell_type [names(clone2[which(clone2 %in% shared_clones[c])])]))
          for(c1 in c(1:length(x1))){
            for(c2 in c(1:length(x2))){
              nam = concat(c(x1[c1]," - ",x2[c2]))
              m_pat_shared_freq[nam,Patients[s]] = m_pat_shared_freq[nam,Patients[s]] +1
            }}
        }
      }
      m_pat_shared_freq[,Patients[s]] = m_pat_shared_freq[,Patients[s]]*100/length(w2)
    }
    
    analysis = concat(c(names(types)[cd],"_cell_type_biopsy_blood_overlap"))
    out_file_table = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
    write.table(m_pat_shared_freq, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    
    analysis = concat(c(names(types)[cd],"_cell_type_biopsy_blood_overlap"))
    file = concat(c(output_directory,"Overlap_",analysis,"_information_",type,"_PDAC150Ka.txt"))
    p <- as.matrix(read.csv(file, head=T, sep="\t"))
    m_pat_shared_freq = p
    repeats = 10000
    
    m_pat_shared_freq1 = m_pat_shared_freq[which(rowSums(m_pat_shared_freq)>0),]
    for(i in c(1:length(m_pat_shared_freq1[1,]))){m_pat_shared_freq1[,i] = m_pat_shared_freq1[,i]*100/sum(m_pat_shared_freq1[,i])}
    
    
    Means_factor = function(factor, x){
      m = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        m = c(m, mean(x1))}
      return(m)}
    
    mat_stat = t(m_pat_shared_freq1)
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, gsub("-biopsy","",groups_PCA[[i]]))
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    
    mat1 = mat_stat[gsub("_biopsy","",group_PCA_list),]
    wx = which(apply(mat_stat, 2, function(x){length(which(x!=0))})>=3)
    wx = which(is.na(rowSums(mat_stat))==F)
    
    #for(i in c(1:length(mat1[,1]))){mat1[i,] = mat1[i,]+runif(length(mat1[i,]))*1e-10}
    
    factor = factor(factor_PCA)
    fit = manova(formula = mat1[wx,]^0.7 ~ factor[wx])
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
      }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    print(sort(p_value))
    
    colnames(means) = paste("mean.group.", c(1:length(means[1,])))
    combined_p_value = cbind(p_value ,means)
    rownames(combined_p_value) = nam
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_tables = x
    
    table(apply(mat1, 2, function(x){length(which(x!=0))}))
    
    fileout1=concat(c(output_directory,"Overlap_",analysis,"_PCA_groups_boxplots_", batch,".pdf"))
    w=3.35
    if(cd==1){pdf(file=fileout1, height=w*1.7*1, width=w*1*6.5)}
    if(cd==2){pdf(file=fileout1, height=w*1.7*1, width=w*1*4.5)}
    par(mfrow= c(1,1), mar = c(18,5,3,3))
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)
    
    mat1 = mat_stat
    wx = which(apply(mat_stat, 2, function(x){length(which(x!=0))})>=3)
    mat1 = mat1[,wx]
    groups = NULL
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat1[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        x = mat1 [gsub("_biopsy","",groups_PCA[[g]]), i]
        x = x[which(x!=-1)]
        data = c(data,x )
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
    }
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = gsub("_"," ",colnames(mat1))
    factors = gsub("activated","act.", factors)
    factors = gsub("memory","mem.", factors)
    factors = gsub("Ag presenting","Ag pres.", factors)
    factors = gsub("SHM diversifying","SHM div.", factors)
    
    main = concat(c("Relative clone sharing blood biopsy"))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = length(groups)
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    if(range <0.0002){scale = c(0:100)*0.00002}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 0.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.18
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.7)
      }}
    
    for(i1 in c(1:length(p_values))){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i1]<signif_threshold){
        i = which(colnames(mat1)==names(p_values)[i1])
        pval1 = "*"
        if(p_values[i] <signif_threshold/10){pval1 = "**"}
        if(p_values[i] <signif_threshold/100){pval1 = "***"}
        y = max(unlist(groups[[i]]))
        y = y+1*b
        # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
        text(i, y+2*b, labels = pval1, cex = 1.3)
      }
    }
    
    dev.off()
    
    ################################################################
    #overall clonal sharing: 
    m_pat_shared_freq1 = m_pat_shared_freq
    for(i in c(1:length(m_pat_shared_freq1[1,]))){m_pat_shared_freq1[,i] = m_pat_shared_freq1[,i]*100/sum(m_pat_shared_freq1[,i])}
    wx = which(is.na(colSums(m_pat_shared_freq1))==F)
    m_pat_shared_freq1 = m_pat_shared_freq1[,wx]
    mat_stat = t(m_pat_shared_freq1)
    width = 1
    height = -1*length(cell_types)
    
    fileout1=concat(c(output_directory,"Overlap_",analysis,"_network_", batch,".pdf"))
    w=3.35
    if(cd==1){pdf(file=fileout1, height=w*0.95, width=w*1.5)}
    if(cd==2){pdf(file=fileout1, height=w*0.75, width=w*1.5)}
    par(mfrow= c(1,1), mar = c(2,10,2,10))
    
    
    plot(c(0, width),c(-0.5, height-0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "clone sharing between blood and pancreatic tumour", axes = F)
    
    mtext(side = 2, text = cell_types, line = 0.1,cex= cex-0.1, at = -1*c(1:length(cell_types)), las = 1, font = 1)
    mtext(side = 4, text = cell_types, line = 0.1,cex= cex-0.1, at = -1*c(1:length(cell_types)), las = 1, font = 1)
    
    mtext(side = 2, text = gsub("  "," ",Sample.Types[1],fixed = T), line = 0.1,cex= cex-0.1, at = 0, las = 1, font = 2)
    mtext(side = 4, text = gsub("  "," ",Sample.Types[2],fixed = T), line = 0.1,cex= cex-0.1, at = 0, las = 1, font = 2)
    
    edge_strength = apply(mat_stat,2,function(x){mean(x[which(is.na(x)==F)])})
    edge_strength = edge_strength
    edge_strength = edge_strength*0.9#/max(edge_strength)
    col = add.alpha("darkblue", alpha = 0.5)
    p_val_sub = (edge_strength*0)+2
    names(p_val_sub)==colnames(mat_stat)
    p_val_sub[names(p_values)]=p_values
    cols = rep(add.alpha("black", alpha = 0.8), length(edge_strength))
    names(cols) = names(p_val_sub)
    mean1 = (edge_strength*0)
    mean2 = (edge_strength*0)
    mean1[rownames(summary_tables)] = as.numeric(summary_tables[,"mean.group. 1"])
    mean2[rownames(summary_tables)] = as.numeric(summary_tables[,"mean.group. 2"])
    w = which(p_val_sub <0.05)
    cols[intersect(names(w), names(which(mean1>mean2)))] = add.alpha("red", alpha = 0.8)
    cols[intersect(names(w), names(which(mean1<mean2)))] = add.alpha("blue", alpha = 0.8)
    
    cbind(p_val_sub,mean1,mean2)[w,]
    table(cols)
    
    split = strsplit(colnames(mat_stat)," - ",fixed = T)
    f1 = NULL
    f2 = NULL
    
    for(i in c(1:length(split))){
      from = which(cell_types== split[[i]][1])
      to = which(cell_types== split[[i]][2])
      
      if(length(c(from,to))==2){
        f1 = c(f1, from)
        f2 = c(f2, to)
        #segments(0,-1*from, x1 = width, y1 = -1*to, lwd = 1, col = "black")}
        segments(0,-1*from, x1 = width, y1 = -1*to, lwd = edge_strength[i], col = cols[i])}
    }
    
    dev.off()
    
    ####### by cell type what proportion of cells are immunosurveilling
    
    
    
    list_proportion_shared= readRDS(file = concat(c(output_directory,"Overlap_Proportion_per_cell_type_shared_information_",type1,"_PDAC150Ka.txt")))
    ## tumour
    t1dt1pt2 = list_proportion_shared[[1]]
    total = list_proportion_shared[[2]]
    for(i in c(1:length(t1dt1pt2[,1]))){
      t1dt1pt2[i,which(total[i,]<=10)] = NA
    }
    a = apply(t1dt1pt2, 2, function(x){length(which(is.na(x)))})
    t1dt1pt2 = t1dt1pt2[,which(a<=3)]
    ## blood
    t1dt1pt2_blood = list_proportion_shared[[3]]
    total = list_proportion_shared[[4]]
    for(i in c(1:length(t1dt1pt2[,1]))){
      t1dt1pt2_blood[i,which(total[i,]<=10)] = NA
    }
    a = apply(t1dt1pt2_blood, 2, function(x){length(which(is.na(x)))})
    t1dt1pt2_blood = t1dt1pt2_blood[,which(a<=3)]
    
    props = c(list(t1dt1pt2), list(t1dt1pt2_blood))
    names(props) = c("tumour","blood")
    
    fileout1=concat(c(output_directory,"Overlap_Proportion_per_cell_type_shared_information_",type1,"_PCA_groups_boxplots_", batch,".pdf"))
    w=3.35
    if(cd==1){pdf(file=fileout1, height=w*1.7*1, width=w*1*6.5)}
    if(cd==2){pdf(file=fileout1, height=w*1.7*1, width=w*1*4.5)}
    par(mfrow= c(1,2), mar = c(18,5,3,3))
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)
    
    
    for(ind in c(1:length(props))){
      mat1 = props[[ind]]
      if(ind==1){mat1 = mat1[gsub("_biopsy","-biopsy",group_PCA_list),]}
      if(ind==2){mat1 = mat1[gsub("_biopsy","-blood",group_PCA_list),]}
      name = names(props)[ind]
      factor = factor(factor_PCA)
      fit = manova(formula = mat1 ~ factor)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      
      i1 = 0
      for(i in p1){
        i1 = i1+1
        p_value = c(p_value, i$'Pr(>F)'[1]) 
        if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
        }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      print(sort(p_value))
      
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.site = rep("biopsy", length(nam))
      p.analysis = rep(analysis, length(nam))
      x = cbind(p.site, p.analysis, combined_p_value)
      summary_tables = x
      
      table(apply(mat1, 2, function(x){length(which(x!=0))}))
      
      groups = NULL
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat1[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA))){
          if(ind==1){x = mat1 [gsub("_biopsy","-biopsy",groups_PCA[[g]]), i]}
          if(ind==2){x = mat1 [gsub("_biopsy","-blood",groups_PCA[[g]]), i]}
          x = x[which(x!=-1)]
          data = c(data,x )
          factor = c(factor, rep(g, length(x)))
          g1 = c(g1, list(x))}
        groups = c(groups, list(g1))
      }
      factors1 = paste("group",c(1:length(groups_PCA)))
      factors = gsub("_"," ",colnames(mat1))
      factors = gsub("activated","act.", factors)
      factors = gsub("memory","mem.", factors)
      factors = gsub("Ag presenting","Ag pres.", factors)
      factors = gsub("SHM diversifying","SHM div.", factors)
      
      main = concat(c(name,"\n% clones shared"))
      max = max(c(unlist(groups), unlist(groups))*1.35)
      max = 100
      min = 0
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 20
      max_scale = max
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*5}
      if(range <10){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      if(range <0.0002){scale = c(0:100)*0.00002}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 0.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.18
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.25/max(shift)
      
      for(i in c(1:l)){
        for(i1 in c(1:l1)){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.7)
        }}
      
      for(i1 in c(1:length(p_values))){	
        b = max*0.035
        signif_threshold = 0.05
        if(p_values[i1]<signif_threshold){
          i = which(colnames(mat1)==names(p_values)[i1])
          pval1 = "*"
          if(p_values[i] <signif_threshold/10){pval1 = "**"}
          if(p_values[i] <signif_threshold/100){pval1 = "***"}
          y = max(unlist(groups[[i]]))
          y = y+1*b
          # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
          text(i, y+2*b, labels = pval1, cex = 1.3)
        }
      }
      
    }
    dev.off()
    
    
  }
  
}
Get_proportions_of_immunosurveilling_cell_types_T_cells(VDJ_list_TCR, output_directory, batch, groups_PCA)

Plot_proportions_of_immunosurveilling_cell_types<-function(VDJ_list_TCR, output_directory, batch, groups_PCA){
  ### cell type differences
  analysisa = "cell_type_prop"
  fileout1=concat(c(output_directory,"Immunosurveilling_vv_private_", analysisa,"_boxplots_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.15*3, width=w*1*4.2)
  par(mfrow= c(3,4), mar = c(16,5,3,3))
  
  types = c("B_cells","TCR_CD4","TCR_CD8")
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 1)
  cols =  add.alpha (cols1, alpha = 0.25)
  cols11 = cols1
  summary_tables = NULL
  for(t in types){
    type= t
    ct = "_private"
    file_table1 = concat(c(output_directory,"Overlap_", analysisa, ct,"_information_",type,"_PDAC150Ka.txt"))
    
    ct = "_shared"
    file_table2 = concat(c(output_directory,"Overlap_",  analysisa, ct,"_information_",type,"_PDAC150Ka.txt"))
    
    p1 <- as.matrix(read.csv(file_table1, head=T, sep="\t"))
    p2 <- as.matrix(read.csv(file_table2, head=T, sep="\t"))
    
    col = colnames(p1)
    col = sort(col)
    p1 = p1[,col]
    p2 = p2[,col]
    col = colnames(p1)
    col[which(col=="NK.cell.NK.like.T.cell" )] = "T.cell.NK.like.T.cell" 
    colnames(p1)=col
    colnames(p2)=col
    
    for(pca in c(1:length(groups_PCA))){
      pats = groups_PCA[[pca]]
      list_frequencies = c(list(p2[gsub("_","-",pats),]),list(p1[gsub("_","-",pats),]))
      names_frequencies = c("immunosurveilling","private")
      
      samples_run = rownames(p1)[grep('biopsy',rownames(p1))]
      columns = colnames(p1)
      groups = NULL
      p_values_list = NULL
      means = matrix(data = -1, nrow = length(columns),ncol = length(names_frequencies), dimnames = c(list(columns),list(names_frequencies)))
      for(ind in c(1:length(columns))){
        gr = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(list_frequencies))){
          x = list_frequencies[[g]] [, columns[ind]]
          x = x[which(x!=-1)]
          data = c(data, x)
          factor = c(factor, rep(g, length(x)))
          gr = c(gr, list(x))
          means[ind, g] = mean(x)}
        if(length(gr[[1]])>=3){
          p = wilcox.test(gr[[1]], y = gr[[2]],alternative = "two.sided",paired = F)$p.value
          p_values_list = c(p_values_list, p)
        }
        groups = c(groups, list(gr))
      }		
      x = cbind(columns, rep(type , length(columns)), p_values_list, means)
      #if(length(summary_tables)==0){summary_tables =x
      #}else{summary_tables =rbind(summary_tables, x)}
      factors = columns
      factors = gsub("."," ", factors, fixed = T)
      factors = gsub("_"," ", factors, fixed = T)
      factors = gsub("  "," ", factors, fixed = T)
      main = concat(c(type," group ",pca))
      max = 135#max(c(unlist(groups), unlist(groups))*1.35)
      min = 0
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 15
      max_scale = 90
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*10}
      if(range<=30){scale = c(0:100)*5}
      if(range <10){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.15
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.2/max(shift)
      library(RColorBrewer)
      cols1 =  add.alpha (brewer.pal(6, "Paired")[c(6,2)], alpha = 0.95)
      cols =  add.alpha (cols1, alpha = 0.5)
      #cols1 =  c(cols11[pca], add.alpha (cols11[pca], alpha = 0.25))
      #cols =  add.alpha (cols1, alpha = 0.25)
      
      for(i in c(1:l)){
        for(i1 in c(1:length(groups[[i]]))){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1], cex = 0.7)
        }}
      if(length(p_values_list)>=1){
        for(i in c(1:l)){	
          b = max*0.035
          signif_threshold = 0.05
          if(is.na(p_values_list[i])==F){
            if(p_values_list[i]<signif_threshold){
              pval1 = "*"
              y = max(unlist(groups[[i]]))
              y = y+1*b
              text(i, y+2*b, labels = pval1, cex = 1.3)
            }}
        }
      }
      plot(c(0,1),c(0,1),pch = 21, col = "white",bg = "white", main = "", xlab = "", ylab = "",axes = "F")
      legend("topleft", names_frequencies, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols1, pt.lwd = 2, text.font = 1)
    }
  }
  
  dev.off()
  
  summary_tables1 = summary_tables
  
  ### isotype differences
  analysisa = "isotype"
  analysis_all= "isotype"
  fileout1=concat(c(output_directory,"Immunosurveilling_vv_private_", analysisa,"_boxplots_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.05*3, width=w*1*2.9)
  par(mfrow= c(3,1), mar = c(15,5,3,3))
  
  types = c("BCR")
  summary_tables = NULL
  for(t in types){
    type= t
    ct = "_private"
    file_table1 = concat(c(output_directory,"Overlap_", analysisa, ct,"_information_",type,"_PDAC150Ka.txt"))
    
    ct = "_shared"
    file_table2 = concat(c(output_directory,"Overlap_",analysisa, ct,"_information_",type,"_PDAC150Ka.txt"))
    
    p1 <- as.matrix(read.csv(file_table1, head=T, sep="\t"))
    p2 <- as.matrix(read.csv(file_table2, head=T, sep="\t"))
    p1 = p1[,which(colnames(p1)!="X." )]
    p2 = p2[,which(colnames(p2)!="X." )]
    for(pca in c(1:length(groups_PCA))){
      pats = groups_PCA[[pca]]
      list_frequencies = c(list(p2[gsub("_","-",pats),]),list(p1[gsub("_","-",pats),]))
      names_frequencies = c("immunosurveilling","private")
      
      samples_run = rownames(p1)[grep('biopsy',rownames(p1))]
      columns = colnames(p1)
      
      groups = NULL
      p_values_list = NULL
      means = matrix(data = -1, nrow = length(columns),ncol = length(names_frequencies), dimnames = c(list(columns),list(names_frequencies)))
      for(ind in c(1:length(columns))){
        gr = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(list_frequencies))){
          x = list_frequencies[[g]] [, columns[ind]]
          x = x[which(x!=-1)]
          data = c(data, x)
          factor = c(factor, rep(g, length(x)))
          gr = c(gr, list(x))
          means[ind, g] = mean(x)}
        p = wilcox.test(gr[[1]], y = gr[[2]],alternative = "two.sided",paired = F)$p.value
        p_values_list = c(p_values_list, p)
        groups = c(groups, list(gr))
      }		
      #x = cbind(columns, rep(type , length(columns)), p_values_list, means)
      #if(length(summary_tables)==0){summary_tables =x
      #}else{summary_tables =rbind(summary_tables, x)}
      factors = columns
      main = concat(c(type," group ",pca))
      max = max(c(unlist(groups), unlist(groups))*1.35)
      min = 0
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 40
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*10}
      if(range<=30){scale = c(0:100)*5}
      if(range <10){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(1.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.15
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.2/max(shift)
      library(RColorBrewer)
      cols1 =  c(cols11[pca], add.alpha (cols11[pca], alpha = 0.25))
      cols =  add.alpha (cols1, alpha = 0.25)
      
      for(i in c(1:l)){
        for(i1 in c(1:length(groups[[i]]))){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1], cex = 0.7)
        }}
      for(i in c(1:l)){	
        b = max*0.035
        signif_threshold = 0.05
        if(is.na(p_values_list[i])==F){
          if(p_values_list[i]<signif_threshold){
            pval1 = "*"
            # if(p_values_list[i] <signif_threshold/10){pval1 = "**"}
            # if(p_values_list[i] <signif_threshold/100){pval1 = "***"}
            y = max(unlist(groups[[i]]))
            y = y+1*b
            # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
            text(i, y+2*b, labels = pval1, cex = 1.3)
          }}
      }
    }
    plot(c(0,1),c(0,1),pch = 21, col = "white",bg = "white", main = "", xlab = "", ylab = "",axes = "F")
    legend("topleft", names_frequencies, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols1, pt.lwd = 2, text.font = 1)
  }
  
  
  dev.off()
  
  summary_tables2 = summary_tables
  
  
  
  ### V gene usage differences
  analysisa = "V_genes"
  
  fileout1=concat(c(output_directory,"Immunosurveilling_vv_private_", analysisa,"_boxplots_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.05*3, width=w*1*5.1*2)
  par(mfrow= c(3,2), mar = c(15,5,3,3))
  
  types = c("BCR","TCR_CD4","TCR_CD8")
  summary_tables = NULL
  for(t in types){
    type= t
    ct = "_private"
    file_table1 = concat(c(output_directory,"Overlap_", analysisa, ct,"_information_",type,"_PDAC150Ka.txt"))
    
    ct = "_shared"
    file_table2 = concat(c(output_directory,"Overlap_", analysisa, ct,"_information_",type,"_PDAC150Ka.txt"))
    
    p1 <- as.matrix(read.csv(file_table1, head=T, sep="\t"))
    p2 <- as.matrix(read.csv(file_table2, head=T, sep="\t"))
    p1= p1[,which(colnames(p1) %in% c("X.","None" ) ==F)]
    p2= p2[,which(colnames(p2) %in% c("X.","None" ) ==F)]
    for(pca in c(1:length(groups_PCA))){
      pats = groups_PCA[[pca]]
      list_frequencies = c(list(p2[gsub("_","-",pats),]),list(p1[gsub("_","-",pats),]))
      names_frequencies = c("immunosurveilling","private")
      
      samples_run = rownames(p1)[grep('biopsy',rownames(p1))]
      columns = colnames(p1)
      
      groups = NULL
      p_values_list = NULL
      means = matrix(data = -1, nrow = length(columns),ncol = length(names_frequencies), dimnames = c(list(columns),list(names_frequencies)))
      for(ind in c(1:length(columns))){
        gr = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(list_frequencies))){
          x = list_frequencies[[g]] [, columns[ind]]
          x = x[which(x!=-1)]
          data = c(data, x)
          factor = c(factor, rep(g, length(x)))
          gr = c(gr, list(x))
          means[ind, g] = mean(x)}
        p = wilcox.test(gr[[1]], y = gr[[2]],alternative = "two.sided",paired = F)$p.value
        p_values_list = c(p_values_list, p)
        groups = c(groups, list(gr))
      }		
      x = cbind(columns, rep(type , length(columns)), p_values_list, means)
      #if(length(summary_tables)==0){summary_tables =x
      #}else{summary_tables =rbind(summary_tables, x)}
      factors = columns
      factors = gsub(".","-", factors, fixed = T)
      main = concat(c(type," group ",pca))
      max = 135#max(c(unlist(groups), unlist(groups))*1.35)
      min = 0
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 100
      max_scale = max
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*10}
      if(range<=30){scale = c(0:100)*5}
      if(range <10){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(4, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.15
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.2/max(shift)
      library(RColorBrewer)
      cols1 =  add.alpha (brewer.pal(6, "Paired")[c(6,2)], alpha = 0.95)
      cols =  add.alpha (cols1, alpha = 0.5)
      
      
      for(i in c(1:l)){
        for(i1 in c(1:length(groups[[i]]))){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1], cex = 0.7)
        }}
      for(i in c(1:l)){	
        b = max*0.035
        signif_threshold = 0.05
        if(is.na(p_values_list[i])==F){
          if(p_values_list[i]<signif_threshold){
            pval1 = "*"
            # if(p_values_list[i] <signif_threshold/10){pval1 = "**"}
            # if(p_values_list[i] <signif_threshold/100){pval1 = "***"}
            y = max(unlist(groups[[i]]))
            y = y+1*b
            # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
            text(i, y+2*b, labels = pval1, cex = 1.3)
          }}
      }
    }
  }
  plot(c(0,1),c(0,1),pch = 21, col = "white",bg = "white", main = "", xlab = "", ylab = "",axes = "F")
  legend("topleft", names_frequencies, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols1, pt.lwd = 2, text.font = 1)
  
  dev.off()
  summary_tables3 = summary_tables
  
}
Plot_proportions_of_immunosurveilling_cell_types(VDJ_list_TCR, output_directory, batch, groups_PCA)




