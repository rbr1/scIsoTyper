#### Run on Oxford compute
# srun -p short --cpus-per-task 2 --pty bash
# module purge
# module load HDF5/1.10.5-gompi-2019a
# module load umap-learn/0.3.10-foss-2019a-Python-3.7.2
# module load Seurat/3.1.2-foss-2019a-R-3.6.0
# module load Harmony/1.0.0-foss-2019a-R-3.6.0
# R


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
  out_dir = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/VDJ/"
  out_dir_raw = "/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/VDJ/"
}

# or for Sakina: 
a = 1
if(a==1){
  file="/gpfs3/well/immune-rep/shared/10X_GENOMICS/P210298_SA/P210298/Samples_PDAC_LT.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,1]!=""),]
  sample_id = as.character(p[,"Sample_Name"])
  sample_output_id = as.character(p[,"Sample_Name"])
  GEX.location = as.character(p[,"Location_of_GEX"])
  CITE.location = as.character(p[,"Location_of_CITE"])
  BCR.location = as.character(p[,"Location_of_BCR"])
  BCR.location = gsub("all_contig","filtered_contig",BCR.location)
  TCR.location = as.character(p[,"Location_of_TCR"])
  TCR.location = gsub("all_contig","filtered_contig",TCR.location)
  Overall_sample_group = as.character(p[,"Patient"])
  Site = as.character(p[,"Sample_type"])
  batch = "PDAC_LTS"
  out_dir = "/gpfs3/well/immune-rep/shared/10X_GENOMICS/P210298_SA/P210298/LT_wdir/"
  out_dir_raw = "/gpfs3/well/immune-rep/shared/10X_GENOMICS/P210298_SA/P210298/LT_wdir/"
}


################# check files exist from input folder
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))

Check_files_exist<-function(TCR.location, BCR.location, sample_id){
  TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
  BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
  
  for(c in BCRs){
    if(is.na(BCR.location[c])==F){
      sample = sample_id[c]
      fasta_file = gsub("_annotations.csv",".fasta",BCR.location[c] )
      csv_file = BCR.location[c]
      if(file.exists(fasta_file)==F){print (fasta_file)}
      if(file.exists(fasta_file)==F){print (csv_file)}
    }}
  
  for(c in TCRs){
    if(is.na(TCR.location[c])==F){
      sample = sample_id[c]
      fasta_file = gsub("_annotations.csv",".fasta",TCR.location[c] )
      csv_file = TCR.location[c]
      if(file.exists(fasta_file)==F){print (fasta_file)}
      if(file.exists(fasta_file)==F){print (csv_file)}
    }}
}

Check_files_exist(TCR.location, BCR.location, sample_id) ## no error messages or printing means that everything is found

################# move files so that source dirs can be closed/compressed

Move_files_to_working_directory<-function(TCR.location, BCR.location, sample_id, out_dir_raw){
  for(c in BCRs){
    if(is.na(BCR.location[c])==F){
      sample = sample_id[c]
      fasta_file = gsub("_annotations.csv",".fasta",BCR.location[c] )
      csv_file = BCR.location[c]
      command = concat(c("cp ", csv_file, " ", out_dir_raw, "filtered_contig_annotations_BCR_", sample,".csv"))
      system(command)
      command = concat(c("cp ", fasta_file, " ", out_dir_raw, "filtered_contig_BCR_", sample,".fasta"))
      system(command)
      print(concat(c(c, " ", sample)))
    }}
  for(c in TCRs){
    if(is.na(TCR.location[c])==F){
      sample = sample_id[c]
      fasta_file = gsub("_annotations.csv",".fasta",TCR.location[c] )
      csv_file = TCR.location[c]
      command = concat(c("cp ", csv_file, " ", out_dir_raw, "filtered_contig_annotations_TCR_", sample,".csv"))
      system(command)
      command = concat(c("cp ", fasta_file, " ", out_dir_raw, "filtered_contig_TCR_", sample,".fasta"))
      system(command)
      print(concat(c(c, " ", sample)))}}
}

Move_files_to_working_directory(TCR.location, BCR.location, sample_id, out_dir_raw) ### this will print out the sample IDs that it has moved

################# batch for IMGT

Group_together_fasta_files_for_IMGT<-function(TCR.location, BCR.location, sample_id, out_dir_raw,batch){
  # make a mapping file that shortens the sample IDs for unambiguous assignment of sequences
  out = cbind(sample_id, BCR.location, TCR.location,apply(cbind(c(1:length(sample_id)),substr(sample_id, start = 1, stop = 7)), 1, paste, collapse = "_"))
  colnames(out) = c("sample","BCR location","TCR location","unique_short_ID")
  outfile = concat(c(out_dir_raw, "IMGT_sample_ID_mapping.txt"))
  write.table(out, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  
  library(ape)
  outfile = concat(c(out_dir_raw, "All_filtered_contig_BCR_", batch,"_"))
  index = 1
  total_seqs = NULL
  max_n_sequences_per_file = 1000000-10
  for(c in BCRs){
    if(is.na(BCR.location[c])==F){
      sample = sample_id[c]
      unique_ID = out[match(sample,out[,1]),"unique_short_ID"]
      fasta_file = concat(c(out_dir_raw, "filtered_contig_BCR_", sample,".fasta"))
      seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
      names(seqs) = paste0(names(seqs),"||", unique_ID, sep = "")
      if( length(seqs) + length(total_seqs)> max_n_sequences_per_file){
        write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE)
        index = index+1
        total_seqs = NULL
      }
      total_seqs = c(total_seqs, seqs)
      print (concat(c( "BCR total sequences: ",length(total_seqs), ", number of files read: ",c)))
    }}
  write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE,colsep = "")
  print(concat(c("Location of BCR grouped sequence file: ", concat(c(outfile,index,".fasta")))))
  
  outfile = concat(c(out_dir_raw, "All_filtered_contig_TCR_", batch,"_"))
  index = 1
  total_seqs = NULL
  max_n_sequences_per_file = 1000000-10
  for(c in TCRs){
    if(is.na(TCR.location[c])==F){
      sample = sample_id[c]
      unique_ID = out[match(sample,out[,1]),"unique_short_ID"]
      fasta_file = concat(c(out_dir_raw, "filtered_contig_TCR_", sample,".fasta"))
      seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
      names(seqs) = paste0(names(seqs),"||", unique_ID, sep = "")
      if( length(seqs) + length(total_seqs)> max_n_sequences_per_file){
        write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE)
        index = index+1
        total_seqs = NULL
      }
      total_seqs = c(total_seqs, seqs)
      print (concat(c( "TCR total sequences: ",length(total_seqs), ", number of files read: ",c)))
    }}
  write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE,colsep = "")
  print(concat(c("Location of TCR grouped sequence file: ", concat(c(outfile,index,".fasta")))))
}

Group_together_fasta_files_for_IMGT(TCR.location, BCR.location, sample_id, out_dir_raw,batch) ## this will print out the cumulative number of sequences per sample and the location of the output file


############## upload the above files to IMGT for annotation
############## once processed, please download and unpackage into a directory specified below
############## the code below assumes that you only needed 1 batch run for the BCR and one for the TCR
dir_IMGT_BCR = concat(c(out_dir,"IMGT_BCR1/"))
dir_IMGT_TCR = concat(c(out_dir,"IMGT_TCR1/"))

Unpacking_guidance<-function(dir_IMGT_BCR, dir_IMGT_TCR){
  print(concat(c("For the BCR IMGT output, run:")))
  print(concat(c("tar Jxvf <file.txz> -C ",dir_IMGT_BCR)))
  print(concat(c("For the TCR IMGT output, run:")))
  print(concat(c("tar Jxvf <file.txz> -C ",dir_IMGT_TCR)))
}
Unpacking_guidance(dir_IMGT_BCR, dir_IMGT_TCR)

TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))

############## read through key IMGT files and split by sample

Map_annotations_by_sample<-function(dir_IMGT_TCR, BCRs, TCRs,dir_IMGT_BCR, sample_id){
  file = concat(c(out_dir_raw, "IMGT_sample_ID_mapping.txt"))
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  short_ID = p1[,"unique_short_ID"]
  names(short_ID) = p1[,1]
  long_ID = p1[,1]
  names(long_ID) = short_ID
  
  types = c("BCR","TCR")
  files = c(concat(c(dir_IMGT_BCR,"2_IMGT-gapped-nt-sequences.txt")),concat(c(dir_IMGT_TCR,"2_IMGT-gapped-nt-sequences.txt")))
  for(f in c(1:length(files))){
    type = types[f]
    p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
    p1 = p[sort(c(which(p[,"V.D.J.REGION"]!=''),which(p[,"V.J.REGION"]!=''))),]
    id = as.character(p1[,"Sequence.ID"])
    seq = toupper (as.character(p1[,"V.D.J.REGION"]))
    seq = gsub(".","",seq,fixed = T)
    seq1 = toupper(as.character(p1[,"V.J.REGION"]))
    
    orig_sample = strsplit(id, "||", fixed = T)
    orig_id = orig_sample
    for(i in c(1:length(orig_sample))){
      orig_id[i] = orig_id[[i]][1]
      orig_sample[i] = orig_sample[[i]][2]}
    orig_sample = unlist(orig_sample)
    
    orig_id = unlist(orig_id)
    orig_id = strsplit(orig_id, "__", fixed = T)
    for(i in c(1:length(orig_id))){
      orig_id[i] = orig_id[[i]][1]}
    orig_id = unlist(orig_id)
    
    sample_long_ID = long_ID[orig_sample]
    
    ### split by sample
    used = NULL
    found = NULL
    library(ape)
    for(c in BCRs){
      if(is.na(BCR.location[c])==F){
        sample = sample_id[c]
        fasta_file = concat(c(out_dir_raw, "/filtered_contig_", type,"_", sample,".fasta"))
        seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
        if(type=="BCR"){csv_file = BCR.location[c]}
        if(type=="TCR"){csv_file = TCR.location[c]}
        p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
        seq_names = p2[,"contig_id"]
        w = which(sample_long_ID== sample)
        w1 = which(orig_id %in% seq_names)
        seq_inds = intersect(w,w1)
        found = c(found, length(seq_inds))
        if(length(seq_inds)>0){
          total_seqs = seq[seq_inds]
          names(total_seqs) = orig_id[seq_inds]
          outfile = concat(c(out_dir_raw, "/IMGT_filtered_contig_", type,"_", sample,".fasta"))
          write.dna(total_seqs, outfile, format = "fasta", append = FALSE, nbcol = -1,colw = 10000000)
          used = c(used,seq_inds)
          print (c(sample, length(w1)))
        }}}
    
    print(concat(c(files[f]," run completed")))
  }
  
  types = c("BCR","TCR")
  files = c(concat(c(dir_IMGT_BCR,"8_V-REGION-nt-mutation-statistics.txt")),concat(c(dir_IMGT_TCR,"8_V-REGION-nt-mutation-statistics.txt")))
  files2 = c(concat(c(dir_IMGT_BCR,"6_Junction.txt")),concat(c(dir_IMGT_TCR,"6_Junction.txt")))
  
  for(f in c(1:length(files))){
    type = types[f]
    p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
    p = p[which(p[,"V.DOMAIN.Functionality"]!="No results"),]
    id = as.character(p[,"Sequence.ID"])
    v_mm = p[,"V.REGION.Nb.of.mutations"]
    v_mm =strsplit(v_mm," ",fixed = T)
    for(i in c(1:length(v_mm))){
      v_mm[i] = v_mm[[i]][1]}
    v_mm = unlist(v_mm)
    names(v_mm) = id
    
    p <- as.matrix(read.csv(files2[f], head=T, sep="\t"))
    p1 = p[which(p[,"V.DOMAIN.Functionality"]!="No results"),]
    id = as.character(p1[,"Sequence.ID"])
    junction_aa = p1[,"JUNCTION..AA."]
    junction_nn = p1[,"JUNCTION"]
    v = p1[,"V.GENE.and.allele"]
    j = p1[,"J.GENE.and.allele"]
    v =strsplit(v," ",fixed = T)
    j =strsplit(j," ",fixed = T)
    for(i in c(1:length(v))){
      v[i] = v[[i]][2]
      j[i] = j[[i]][2]}
    v = unlist(v)
    j = unlist(j)
    v_mms = v_mm[id]
    functionality = p1[,"V.DOMAIN.Functionality"]
    x1 = cbind(junction_aa, junction_nn,v,j, v_mms,functionality)
    rownames(x1) = id
    
    
    orig_sample = strsplit(id, "||", fixed = T)
    orig_id = orig_sample
    for(i in c(1:length(orig_sample))){
      orig_id[i] = orig_id[[i]][1]
      orig_sample[i] = orig_sample[[i]][2]}
    orig_sample = unlist(orig_sample)
    
    sample_long_ID = long_ID[orig_sample]
    
    ### split by sample
    used = NULL
    found = NULL
    library(ape)
    for(c in BCRs){
      if(is.na(BCR.location[c])==F){
        sample = sample_id[c]
        fasta_file = concat(c(out_dir_raw, "/filtered_contig_", type,"_", sample,".fasta"))
        seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
        if(type=="BCR"){csv_file = BCR.location[c]}
        if(type=="TCR"){csv_file = TCR.location[c]}
        p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
        seq_names = p2[,"contig_id"]
        w = which(sample_long_ID== sample)
        w1 = which(orig_id %in% seq_names)
        seq_inds = intersect(w,w1)
        x = x1[seq_inds,]
        name = orig_id[seq_inds]
        found = c(found, length(seq_inds))
        outfile = concat(c(out_dir_raw, "IMGT_filtered_annotation_", type,"_", sample,".txt"))
        write.table(x, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
        print (c(sample, length(w1)))
      }}
    print(concat(c(files[f]," run completed")))
    
  }
  
  types = c("BCR","TCR")
  files = c(concat(c(dir_IMGT_BCR,"5_AA-sequences.txt")),concat(c(dir_IMGT_TCR,"5_AA-sequences.txt")))
  for(f in c(1:length(files))){
    type = types[f]
    p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
    p1 = p[sort(c(which(p[,"V.D.J.REGION"]!=''),which(p[,"V.J.REGION"]!=''))),]
    id = as.character(p1[,"Sequence.ID"])
    seq = toupper (as.character(p1[,"V.D.J.REGION"]))
    seq = gsub(".","",seq,fixed = T)
    seq1 = toupper(as.character(p1[,"V.J.REGION"]))
    w = which(seq=='')
    seq[w] = seq1[w]
    w = which(seq=='')
    
    orig_sample = strsplit(id, "||", fixed = T)
    w = setdiff(c(1:length(id)),  grep("||", id, fixed = T))
    orig_id = orig_sample
    for(i in c(1:length(orig_sample))){
      orig_id[i] = orig_id[[i]][1]
      orig_sample[i] = orig_sample[[i]][2]}
    orig_sample = unlist(orig_sample)
    orig_sample1 = strsplit(id, "__", fixed = T)
    orig_sample = unlist(orig_sample)
    
    sample_long_ID = long_ID[orig_sample]
    
    ### split by sample
    found = NULL
    used = NULL
    library(ape)
    for(c in BCRs){
      if(is.na(BCR.location[c])==F){
        sample = sample_id[c]
        fasta_file = concat(c(out_dir_raw, "/filtered_contig_", type,"_", sample,".fasta"))
        seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
        if(type=="BCR"){csv_file = BCR.location[c]}
        if(type=="TCR"){csv_file = TCR.location[c]}
        p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
        seq_names = p2[,"contig_id"]
        w = which(sample_long_ID== sample)
        w1 = which(orig_id %in% seq_names)
        seq_inds = intersect(w,w1)
        total_seqs = seq[seq_inds]
        names(total_seqs) = orig_id[seq_inds]
        found = c(found, length(seq_inds))
        if(length(seq_inds)>0){
          outfile = concat(c(out_dir_raw, "IMGT_filtered_amino_acids_", type,"_", sample,".fasta"))
          write.dna(total_seqs, outfile, format = "fasta", append = FALSE, nbcol = -1,colw = 10000000)
          used = c(used,seq_inds)
          print (c(sample, length(w1)))
        }}}
    print(concat(c(files[f]," run completed")))
  }
  
}

## this will take a while so I advise to get a cup of tea now
Map_annotations_by_sample(dir_IMGT_TCR, BCRs,TCRs, dir_IMGT_BCR, sample_id)

############## gather together BCR IgH/L or TCRA/B chain information
Gather_VDJ_information<-function(sample, fasta_file, csv_file, gene, filtered_fasta_contig_file, 
                                 immune_cell_annotation_file, trimmed_sequences_VDJ, IMGT_filtered_annotation, type){
  if(type=="BCR"){chains = c("IGH","IGL")}
  if(type=="TCR"){chains = c("TRA","TRB")}
  ##
  p0 <- as.matrix(read.csv(csv_file, head=T, sep=","))
  cell_id = p0[,"barcode"]
  chain = p0[,"chain"]
  chain[which(chain %in% c("IGK"))] = "IGL"
  rownames(p0) = cell_id
  contig = p0[,"contig_id"]
  
  ###
  library(ape)
  seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
  inter = intersect(contig, names(seqs))
  ###
  p1 <- as.matrix(read.csv(IMGT_filtered_annotation, head=T, sep="\t"))
  v_IMGT = strsplit(p1[,"v"],"*", fixed = T)
  j_IMGT = strsplit(p1[,"j"],"*", fixed = T)
  cont = strsplit(rownames(p1),"||",fixed = T)
  for(i in c(1:length(v_IMGT))){
    v_IMGT[i] = v_IMGT[[i]][1]
    j_IMGT[i] = j_IMGT[[i]][1]
    cont[i] = cont[[i]][1]}
  v_IMGT = unlist(v_IMGT)
  j_IMGT = unlist(j_IMGT)
  cont = unlist(cont)
  p1[,"v"] = v_IMGT
  p1[,"j"] = j_IMGT
  rownames(p1) = cont
  
  ### put together with matched receptors per cell
  cell_ids = sort(unique(cell_id))
  headers = c("#cell","contig1","chain1","constant_region1","n_umis1","V_gene_10X1","J_gene_10X1","cdr3_aa1", "cdr3_nn1","V_gene1","J_gene1", "V_mm1", "chain_functionality1","mixed_contig_chain1","mixed_contig_n_umis1","contig2","chain2","constant_region2","n_umis2","V_gene2","J_gene2","cdr3_aa2", "cdr3_nn2", "V_gene2","J_gene2","V_mm2", "chain_functionality2","mixed_contig_chain2","mixed_contig_n_umis2")
  mat_cells = matrix(data = "-", nrow = length(cell_ids), ncol = length(headers), dimnames = c(list(cell_ids), list (headers)))
  mat_cells[,"#cell"] = cell_ids
  w_use = which(p0[,"contig_id"] %in% cont ==T)
  w_chain1 = intersect(which(chain ==chains[1]),w_use)
  w_chain2 = intersect(which(chain ==chains[2]),w_use)
  for(c in c(1:length(cell_ids))){
    w = which(cell_id==cell_ids[c])
    ## chain 1: 
    w1 = intersect(w, w_chain1)
    if(length(w1)==1){
      mat_cells[cell_ids[c],c("contig1","chain1","constant_region1","n_umis1","V_gene_10X1",
                              "J_gene_10X1","cdr3_aa1", "cdr3_nn1","V_gene1","J_gene1", "V_mm1", "chain_functionality1")]= 
        c(p0[w1,c("contig_id","chain","c_gene","umis","v_gene","j_gene")], p1[p0[w1,"contig_id"],], p0[w1,'productive'])
    }
    if(length(w1)>1){
      break
      px = p0[w1,]
      px = px[order(as.numeric(px[,'umis']), decreasing = T),]
      mat_cells[cell_ids[c],c("contig1","chain1","constant_region1","n_umis1","V_gene_10X1",
                              "J_gene_10X1","cdr3_aa1", "cdr3_nn1","V_gene1","J_gene1", "V_mm1", "chain_functionality1")]= 
        c(px[1,c("contig_id","chain","c_gene","umis","v_gene","j_gene")], 
          p1[px[1,"contig_id"],],px[1,'productive'])
      mat_cells[cell_ids[c],c("mixed_contig_chain1","mixed_contig_n_umis1")]= c(px[2,c("contig_id","umis")])
    }
    ## chain 2: 
    w1 = intersect(w, w_chain2)
    if(length(w1)==1){
      mat_cells[cell_ids[c],c("contig2","chain2","constant_region2","n_umis2","V_gene2","J_gene2","cdr3_aa2", "cdr3_nn2", "V_gene2","J_gene2","V_mm2", "chain_functionality2")]= 
        c(p0[w1,c("contig_id","chain","c_gene","umis","v_gene","j_gene")], p1[p0[w1,"contig_id"],], p0[w1,'productive'])
    }
    if(length(w1)>1){
      px = p0[w1,]
      px = px[order(as.numeric(px[,'umis']), decreasing = T),]
      mat_cells[cell_ids[c],c("contig2","chain2","constant_region2","n_umis2","V_gene2","J_gene2","cdr3_aa2", "cdr3_nn2", "V_gene2","J_gene2","V_mm2", "chain_functionality2")]= 
        c(px[1,c("contig_id","chain","c_gene","umis","v_gene","j_gene")], 
          p1[px[1,"contig_id"],],px[1,'productive'])
      mat_cells[cell_ids[c],c("mixed_contig_chain2","mixed_contig_n_umis2")]= c(px[2,c("contig_id","umis")])
    }
  }
  write.table(mat_cells, file = immune_cell_annotation_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}

BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
type = "BCR"
for(c in BCRs){
  if(is.na(BCR.location[c])==F){
    sample = sample_id[c]
    fasta_file = concat(c(out_dir_raw, "/IMGT_filtered_contig_", type,"_", sample,".fasta"))
    csv_file = BCR.location[c]
    gene = "IG"
    filtered_fasta_contig_file = concat(c(out_dir_raw, "Sequences_",sample,"_",gene))
    immune_cell_annotation_file = concat(c(out_dir_raw, "Cell_annotation_",sample,"_",gene,".txt"))
    trimmed_sequences_VDJ = concat(c(out_dir_raw, "IMGT_filtered_amino_acids_",sample,"_",gene))
    IMGT_filtered_annotation = concat(c(out_dir_raw,"IMGT_filtered_annotation_",type,"_",sample,".txt"))
    Gather_VDJ_information(sample, fasta_file, csv_file, gene, filtered_fasta_contig_file, immune_cell_annotation_file, trimmed_sequences_VDJ, IMGT_filtered_annotation, type)
    print(sample)
  }
}

TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
type = "TCR"
for(c in TCRs){
  if(is.na(TCR.location[c])==F){
    sample = sample_id[c]
    fasta_file = concat(c(out_dir_raw, "/IMGT_filtered_contig_", type,"_", sample,".fasta"))
    csv_file = TCR.location[c]
    gene = "TCR"
    filtered_fasta_contig_file = concat(c(out_dir_raw, "Sequences_",sample,"_",gene))
    immune_cell_annotation_file = concat(c(out_dir_raw, "Cell_annotation_",sample,"_",gene,".txt"))
    trimmed_sequences_VDJ = concat(c(out_dir_raw, "Trimmed_sequences_",sample,"_",gene))
    IMGT_filtered_annotation = concat(c(out_dir_raw,"IMGT_filtered_annotation_",type,"_",sample,".txt"))
    Gather_VDJ_information(sample, fasta_file, csv_file, gene, filtered_fasta_contig_file, immune_cell_annotation_file, trimmed_sequences_VDJ, IMGT_filtered_annotation, type)
    print(sample)
  }
}

############## clonality analysis TCRs/BCRs per group
Overall_sample_groups = sort(unique(Overall_sample_group))

Group_sequence_files<-function(combined_sequence_file, cell_info_files, 
                               sample_fastas,id,samples_group,amino_acid_fastas,output_file_prefix,cluster_file){
  library(ape)
  library(igraph)
  if(length(grep("IGH",id))!=0){gene = "IGH"}
  if(length(grep("TRA",id))!=0){gene = "TRA"}
  if(length(grep("IGL",id))!=0){gene = "IGL"}
  if(length(grep("TRB",id))!=0){gene = "TRB"}
  
  all_seqs = NULL
  for(f in c(1:length(cell_info_files))){
    p <- as.matrix(read.csv(cell_info_files[f], head=T, sep="\t"))
    seq_id = p[,"X.cell" ]
    if(gene %in% c("IGH","TRA")){contig = p[,"contig1"]}
    if(gene %in% c("IGL","TRB")){contig = p[,"contig2"]}
    CDR3a = gsub("*","-",p[,"cdr3_aa1"],fixed = T)
    CDR3b = gsub("*","-",p[,"cdr3_aa2"],fixed = T)
    
    seqs <- read.dna(amino_acid_fastas[f], format = "fasta", as.character = T)
    inter = intersect(contig, names(seqs))
    seqs = seqs[inter]
    names(seqs) = paste0(names(seqs),"||", samples_group[f], sep = "")
    all_seqs = c(all_seqs, seqs)
    print(f)
  }
  ## combine identical sequences
  unique_seqs = unique(all_seqs)
  unique_ID = apply(cbind("uniq_",c(1:length(unique_seqs))), 1, paste, collapse = "")
  names(unique_seqs) = unique_ID
  unique_ID_match = match(all_seqs, unique_seqs)
  names_unique_map = cbind(names(all_seqs), unique_ID[unique_ID_match])
  write.table(names_unique_map, file = concat(c(output_file_prefix, "_names_unique_map.txt")), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ### check
  #which(table(names_unique_map[,2])>400)
  #all_seqs[names_unique_map[which(names_unique_map[,2]=="uniq_8"),1]]
  
  # write out
  write.dna(unique_seqs, combined_sequence_file, format = "fasta", append = FALSE)
  
  ## run CD-hit
  cd_hit_directory = "/well/immune-rep/shared/CODE/cd-hit-v4.5.7-2011-12-16/"
  command= concat(c(cd_hit_directory,"cd-hit -i ",combined_sequence_file," -o ",tmp_file," -c 0.98 -d 180 -T 10  -M 0 -AL 40 "))
  print (command)
  system(command)
  
  ## read cd-hit output
  output = concat(c(tmp_file, ".clstr"))
  p <- as.matrix(read.csv(output, head=F, sep="\t"))
  seq_id = NULL
  cluster_match = NULL
  for(i in c(1:length(p[,1]))){
    if(length(grep(">",p[i,1]))==1){cluster = gsub(">Cluster ","",p[i,1])
    }else{
      seq_id = c(seq_id, p[i,2])
      cluster_match = c(cluster_match, cluster)
    }
  }
  seq_id = strsplit(seq_id, ">", fixed = T)
  for(i in c(1:length(seq_id))){seq_id[i] = seq_id[[i]][2]}
  seq_id = strsplit(unlist(seq_id), "...", fixed = T)
  for(i in c(1:length(seq_id))){seq_id[i] = seq_id[[i]][1]}
  seq_id = unlist(seq_id)
  clustered_sequences = cbind(seq_id, cluster_match)
  
  a = names_unique_map[which(names_unique_map[,2] %in% clustered_sequences[which(cluster_match==2185),1]),1]
  
  
  ### check distances between sequences in enlarged clusters
  t = table(cluster_match)
  enlarged_clusters = names(which(t>1))
  tmp_file_aln = concat(c(tmp_file, ".aln"))
  tmp_file_aln_out = concat(c(tmp_file, ".aln_out"))
  all_aas = sort(unique(unlist(all_seqs)))
  getKmers = function(sequence, size=3) {
    kmers = c()
    for (x in 1:(nchar(sequence) - size + 1)) {
      kmers = c(kmers, substr(sequence, x, x+size-1))
    }
    return(table(kmers))
  }
  
  edges1 = NULL
  edges2 = NULL
  edge_dist = NULL
  done = c(1:1000)*100
  if(length(enlarged_clusters)>0){
    for(c in c(1:length(enlarged_clusters))){
      if(c %in% done){print(concat(c(c, " of ",length(enlarged_clusters)," clusters done")))}
      w = clustered_sequences[which(cluster_match==enlarged_clusters[c]),1]
      seqs = unique_seqs[w]
      lens =  lapply(seqs, function(x){length(x)})
      kmercounts = lapply(seqs, function(x){getKmers(concat(x), size=3)})
      all_kmers = sort(unique(unlist(lapply(kmercounts, function(x){names(x)}))))
      t = table(all_kmers)*0
      kmercounts_mat = lapply(kmercounts, function(x){t1 = t
      t1[names(x)] = x
      return(t1)
      })
      kmercounts_mat <- do.call(rbind,lapply(kmercounts_mat,matrix,ncol=length(kmercounts_mat[[1]]),byrow=TRUE))
      rownames(kmercounts_mat) = names(seqs)
      colnames(kmercounts_mat) = names(t)
      d = as.matrix(dist(kmercounts_mat, method = "maximum"))
      if(max(d[1,2])<=1){
        for(i1 in c(1:length(seqs))){
          for(i2 in c(i1:length(seqs))){
            if(i1<i2){
              edges1 = c(edges1, names(seqs)[i1])
              edges2 = c(edges2, names(seqs)[i2])
              edge_dist = c(edge_dist, d[i1,i2])
            }}
        }
      }else{
        d[which(d>1)] = 100
        g  <- graph.adjacency(d, weighted=TRUE)
        g_mst <- mst(g)
        edges = as_edgelist(g_mst, names = TRUE)
        edge_strength = E(g_mst)$weight
        w = which(edge_strength<=1)
        if(length(w)>=1){
          edges = cbind(edges[w,1], edges[w,2])
          edges1 = c(edges1, edges[,1])
          edges2 = c(edges2, edges[,2])
          edge_dist = c(edge_dist, edge_strength[w])
        }
      }
    }
  }
  #### make overall network and output edgelist and cluster IDs
  edge_out = cbind(edges1, edges2, edge_dist)
  
  g <- graph.empty(n=0, directed=T)
  unique_seq_ids = names(unique_seqs)
  g <- igraph::add.vertices(g, length(unique_seq_ids), name= unique_seq_ids) 
  names <- V(g)$name
  ids <- 1:length(names)
  names(ids) <- names
  edges <- matrix(c(ids[edges1], ids[edges2]), nc=2)
  g <- add.edges(g, t(edges), weight= edge_dist)
  saveRDS(file = concat(c(output_file_prefix, "_graph.RDS")), g)
  write.table(edge_out, file = concat(c(output_file_prefix, "_edgefile.txt")), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ### connected components
  cluster = components(g)$ membership
  cluster_all = cluster[names_unique_map[,2]]
  cluster_output = cbind(cluster_all, names_unique_map)
  cluster_output = cluster_output[order(as.numeric(cluster_output[,1])),]
  colnames(cluster_output) = c("cluster", "barcode","unique sequence ID")
  write.table(cluster_output, file = cluster_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  print("Done")
}
### BCR
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
for(o in c(1:length(Overall_sample_groups))){
  samples_group = sample_id [intersect(BCRs ,which(Overall_sample_group ==Overall_sample_groups[o]))]
  print(o)
  if(length(samples_group)>0){
    p=1
    #check which have BCRs/TCRs
    grouped_sample_id = Overall_sample_groups[o]
    chains = c("IGH","IGL")
    for(ch in c(1:length(chains))){
      print(ch)
      chain = chains[ch]
      if(p<= length(samples_group)){
        samples_ids = paste0(samples_group[p], collapse = ",")
        sample_fastas = paste0(concat(c(out_dir_raw,"Trimmed_sequences_")), samples_group[p],"_IG_",chain,".fasta" )
        amino_acid_fastas = paste0(concat(c(out_dir_raw, "IMGT_filtered_amino_acids_BCR_")), samples_group[p],".fasta" )
        #info = file.info(sample_fastas)
        #w = which(info$size != 0)
        # sample_fastas = sample_fastas[w]
        #amino_acid_fastas = amino_acid_fastas[w]
        
        cell_info_files = paste0(concat(c(out_dir_raw,"Cell_annotation_")), samples_group[p],"_IG.txt" )
        id = concat(c(grouped_sample_id ,"_",chain))
        
        combined_sequence_file = concat(c( out_dir_raw,"Sequences_combined_",Overall_sample_groups[o],"_",chain,".fasta"))
        tmp_file = concat(c( out_dir_raw,"Tmp_",Overall_sample_groups[o],"_",chain))
        output_file_prefix = concat(c( out_dir_raw,"Edges_",Overall_sample_groups[o],"_",chain))
        cluster_file = concat(c( out_dir_raw,"Cluster_identities_",Overall_sample_groups[o],"_",chain,".txt"))
        
        Group_sequence_files(combined_sequence_file, cell_info_files, sample_fastas,id,samples_group,amino_acid_fastas,output_file_prefix,cluster_file)
        p = p+1
      }
    }
  }
}

### TCR
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
for(o in c(1:length(Overall_sample_groups))){
  #o=5
  samples_group = sample_id [intersect(TCRs ,which(Overall_sample_group ==Overall_sample_groups[o]))]
  print(o)
  if(length(samples_group)>0){
    p=1
    #check which have BCRs/TCRs
    grouped_sample_id = Overall_sample_groups[o]
    chains = c("TRA","TRB")
    for(ch in c(1:length(chains))){
      print(ch)
      chain = chains[ch]
      if(p<= length(samples_group)){
        samples_ids = paste0(samples_group[p], collapse = ",")
        sample_fastas = paste0(concat(c(out_dir_raw,"Trimmed_sequences_")), samples_group[p],"_TCR_",chain,".fasta" )
        amino_acid_fastas = paste0(concat(c(out_dir_raw, "IMGT_filtered_amino_acids_TCR_")), samples_group[p],".fasta" )
        #info = file.info(sample_fastas)
        #w = which(info$size != 0)
        # sample_fastas = sample_fastas[w]
        #amino_acid_fastas = amino_acid_fastas[w]
        
        cell_info_files = paste0(concat(c(out_dir_raw,"Cell_annotation_")), samples_group[p],"_TCR.txt" )
        id = concat(c(grouped_sample_id ,"_",chain))
        
        combined_sequence_file = concat(c( out_dir_raw,"Sequences_combined_",Overall_sample_groups[o],"_",chain,".fasta"))
        tmp_file = concat(c( out_dir_raw,"Tmp_",Overall_sample_groups[o],"_",chain))
        output_file_prefix = concat(c( out_dir_raw,"Edges_",Overall_sample_groups[o],"_",chain))
        cluster_file = concat(c( out_dir_raw,"Cluster_identities_",Overall_sample_groups[o],"_",chain,".txt"))
        
        Group_sequence_files(combined_sequence_file, cell_info_files, sample_fastas,id,samples_group,amino_acid_fastas,output_file_prefix,cluster_file)
        p = p+1
      }
    }
  }
}
################# Get VDJ annotation information object

BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))

labels.vector =c() ### get all cell IDs
for(c in c(1:length(sample_id))){
  if(c %in% BCRs){
    sample = sample_id[c]
    cell_info_BCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_IG.txt" )
    p1 <- as.matrix(read.csv(cell_info_BCR, head=TRUE, sep="\t"))
    cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
    labels.vector = c(labels.vector, cells)}
  if(c %in% TCRs){
    sample = sample_id[c]
    cell_info_TCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_TCR.txt" )
    p1 <- as.matrix(read.csv(cell_info_TCR, head=TRUE, sep="\t"))
    cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
    labels.vector = c(labels.vector, cells)
  }
  print (c)
}
labels.vector = sort(unique(labels.vector))
sample.vector = strsplit(labels.vector, "||", fixed = T)
for(i in c(1:length(sample.vector))){sample.vector[i]=sample.vector[[i]][2]}
sample.vector = unlist(sample.vector)

headers = c( "X.cell","contig1","chain1","constant_region1","n_umis1","V_gene_10X1",
             "J_gene_10X1","cdr3_aa1","cdr3_nn1","V_gene1","J_gene1","V_mm1" , "chain_functionality1","mixed_contig_chain1","mixed_contig_n_umis1","contig2", "chain2", "constant_region2",
             "n_umis2","V_gene2","J_gene2", "cdr3_aa2","cdr3_nn2","V_gene2.1", "J_gene2.1","V_mm2","chain_functionality2","mixed_contig_chain2","mixed_contig_n_umis2")
m_VDJ_BCR = matrix(data = "-", nrow = length(labels.vector), ncol = length(headers), dimnames = c(list(labels.vector), list(headers)))
m_VDJ_TCR = matrix(data = "-", nrow = length(labels.vector), ncol = length(headers), dimnames = c(list(labels.vector), list(headers)))


for(c in c(1:length(sample_id))){
  sample = sample_id[c]
  cell_info_BCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_IG.txt" )
  p1 <- as.matrix(read.csv(cell_info_BCR, head=TRUE, sep="\t"))
  cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
  w = which(cells %in% labels.vector)
  m_VDJ_BCR[cells[w], ] = p1[w,headers]
  print(concat(c(c," ",sample_id[c])))
  print(table(m_VDJ_BCR[,"chain1"]))
}

for(c in c(1:length(sample_id))){
  sample = sample_id[c]
  cell_info_TCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_TCR.txt" )
  p1 <- as.matrix(read.csv(cell_info_TCR, head=TRUE, sep="\t"))
  cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
  w = which(cells %in% labels.vector)
  m_VDJ_TCR[cells[w], ] = p1[w,headers]
  print(concat(c(c," ",sample_id[c])))
  print(table(m_VDJ_TCR[,"chain1"]))
}

clone1 = rep("-", length(labels.vector))
clone2 = rep("-", length(labels.vector))
m_VDJ_BCR = cbind(m_VDJ_BCR, clone1, clone2)
m_VDJ_TCR = cbind(m_VDJ_TCR, clone1, clone2)

Overall_sample_groups = sort(unique(Overall_sample_group))
run = NULL
contig1 = apply(cbind(m_VDJ_TCR[, "contig1"], sample.vector), 1, paste, collapse="||")
contig2 = apply(cbind(m_VDJ_TCR[, "contig2"], sample.vector), 1, paste, collapse="||")
for(o in c(1:length(Overall_sample_groups))){
  cluster_file1 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_TRA.txt" )
  cluster_file2 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_TRB.txt" )
  cluster_files = c(cluster_file1, cluster_file2)
  types = c("TRA","TRB")
  info = file.info(cluster_files)
  w = which(info$size != 0)
  if(length(w)>0){
    for(c in c(1:length(cluster_files))){
      p1 <- as.matrix(read.csv(cluster_files[c], head=F, sep="\t"))
      p1=p1[-1,]
      cells = p1[,2]#gsub("-1",concat(c("||", sample_output_id[c])), p1[,3])
      clone = p1[,1]
      clone = paste (as.numeric(clone), concat(c("||",Overall_sample_groups[o])))
      clone = gsub(" ", "", clone)
      cells = strsplit(cells,"||", fixed = T)
      cell=NULL
      sample_source = NULL
      for(i in c(1:length(cells))){
        cell = c(cell, cells[[i]][1])
        sample_source = c(sample_source, cells[[i]][2])
      }
      m = sample_id [match(sample_source, sample_id)]
      cells = apply(cbind(cell, m), 1, paste, collapse="||")
      names(clone) = cells
      if(types[c] %in% c("TRA")){
        cells_use = cells[which(cells %in% contig1)]
        clone_use = clone[cells_use]
        m_VDJ_TCR[match(cells_use, contig1),"clone1"] = clone_use
      }
      if(types[c] %in% c("TRB")){
        cells_use = cells[which(cells %in% contig2)]
        clone_use = clone[cells_use]
        m_VDJ_TCR[match(cells_use, contig2),"clone2"] = clone_use
      }
    }
  }
  print(concat(c("Number of matched TRA clones: ",length(which(m_VDJ_TCR[,"clone1"]!='-')))))
}

contig1 = apply(cbind(m_VDJ_BCR[, "contig1"], sample.vector), 1, paste, collapse="||")
contig2 = apply(cbind(m_VDJ_BCR[, "contig2"], sample.vector), 1, paste, collapse="||")

for(o in c(1:length(Overall_sample_groups))){
  cluster_file3 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_IGH.txt" )
  cluster_file4 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_IGL.txt" )
  cluster_files = c( cluster_file3, cluster_file4)
  types = c("IGH","IGL")
  info = file.info(cluster_files)
  w = which(info$size != 0)
  if(length(w)>0){
    for(c in c(1:length(cluster_files))){
      p1 <- as.matrix(read.csv(cluster_files[c], head=F, sep="\t"))
      p1=p1[-1,]
      cells = p1[,2]
      clone = p1[,1]
      clone = paste (as.numeric(clone), concat(c("||",Overall_sample_groups[o])))
      clone = gsub(" ", "", clone)
      cells = strsplit(cells,"||", fixed = T)
      cell=NULL
      sample_source = NULL
      for(i in c(1:length(cells))){
        cell = c(cell, cells[[i]][1])
        sample_source = c(sample_source, cells[[i]][2])
      }
      m = sample_id [match(sample_source, sample_id)]
      cells = apply(cbind(cell, m), 1, paste, collapse="||")
      names(clone) = cells
      if(types[c] %in% c("IGH")){
        cells_use = cells[which(cells %in% contig1)]
        clone_use = clone[cells_use]
        m_VDJ_BCR[match(cells_use, contig1),"clone1"] = clone_use
      }
      if(types[c] %in% c("IGL")){
        cells_use = cells[which(cells %in% contig2)]
        clone_use = clone[cells_use]
        m_VDJ_BCR[match(cells_use, contig2),"clone2"] = clone_use
      }
    }
  }
  print(concat(c("Number of matched IGH clones: ",length(which(m_VDJ_BCR[,"clone1"]!='-')))))
}


m_VDJ_BCR[which(m_VDJ_BCR[,"constant_region1"]=="None"),"constant_region1"] = '-'
m_VDJ_BCR[which(m_VDJ_BCR[,"constant_region2"]=="None"),"constant_region2"] = '-'
m_VDJ_TCR[which(m_VDJ_TCR[,"constant_region1"]=="None"),"constant_region1"] = '-'
m_VDJ_TCR[which(m_VDJ_TCR[,"constant_region2"]=="None"),"constant_region2"] = '-'

VDJ_object = c(list(m_VDJ_BCR),list(m_VDJ_TCR))
names(VDJ_object) = c("BCR","TCR")

saveRDS(file=concat(c(out_dir_raw,"/VDJ_information_", batch,".VDJ")), VDJ_object)

## output as matrix file
out_file_table = concat(c(out_dir_raw,"/VDJ_information_BCR_", batch,".txt"))
write.table(m_VDJ_BCR, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
out_file_table = concat(c(out_dir_raw,"/VDJ_information_TCR_", batch,".txt"))
write.table(m_VDJ_TCR, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
print(concat(c("scp -p mfj169@cluster2.bmrc.ox.ac.uk:", out_file_table," ./ " )))


##### check
w1 = which(m_VDJ_TCR[,"J_gene_10X1"]!="-")
w2 = which(m_VDJ_TCR[,"clone1"]!="-")
length(w1)
length(w2)
length(intersect(w1,w2))
head(m_VDJ_TCR[setdiff(w1,w2),])
head(m_VDJ_TCR[w2,])

table(sample.vector[setdiff(w1,w2)])
############## print our cell numbers per sample
VDJ_object = readRDS(file=concat(c(out_dir_raw,"/VDJ_information_", batch,".VDJ")))

TCR = VDJ_object$TCR
BCR = VDJ_object$BCR
sample = strsplit(rownames(BCR),"||", fixed = T)
for(i in c(1:length(sample))){sample[i] = sample[[i]][2]}
sample = unlist(sample)
names(sample) = rownames(BCR)

dIGH = which(BCR[,"chain1"]!='-')
dIGKL = which(BCR[,"chain2"]!='-')
dTRA = which(TCR[,"chain1"]!='-')
dTRB = which(TCR[,"chain2"]!='-')
dIGHKL = intersect(dIGH, dIGKL)
dTRAB = intersect(dTRA, dTRB)

n_cells = cbind(table(sample[dIGH]),table(sample[dIGKL]),table(sample[dIGHKL]),table(sample[dTRA]),table(sample[dTRB]),table(sample[dTRAB]))
colnames(n_cells) = c("IGH+", "IGK/L+", "IGH+IGK/L+", "TRA+","TRB+","TRA+TRB")
out_file_table = concat(c(out_dir_raw,"/VDJ_information_n_droplets_", batch,".txt"))
write.table(n_cells, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
print(concat(c("scp -p mfj169@cluster2.bmrc.ox.ac.uk:", out_file_table," ./ " )))


