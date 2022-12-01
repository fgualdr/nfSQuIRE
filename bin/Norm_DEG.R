#!/usr/bin/env Rscript 
# this scripts use: singularity shell -B /hpcnfs docker://fgualdr/envrnorm
# testing on: /hpcnfs/data/GN2/fgualdrini/Data_Projects/SMARCA_JULIA/RNAseq_SQuIPRE/count squire_normdeg
# opt= list()
# opt$file = "/hpcnfs/data/GN2/fgualdrini/Data_Projects/SMARCA_JULIA/RNAseq_SQuIPRE/count/"
# opt$out = "squire_normdeg"
# opt$processes = 4
# opt$sample_rm = "/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/SMARCA_JULIA/RNAseq_SQuIPRE/RM_SAMPLES"
# opt$deg_design = "/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/SMARCA_JULIA/RNAseq_SQuIPRE/DEG_DESIGN"
# opt$revel_conditions = "WT,UT"

library(optparse)
library(dplyr)
library(GeneralNormalizer)
# library(xlsx)

'%!in%' <- function(x,y)!('%in%'(x,y))

rows_in_lim <- function(m,ld,lu){
        m[m<ld] = 0
        m[m>lu] = 0
        m[m!=0] = 1
        return(rownames(m)[which(rowSums(m) == ncol(m))])
}
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])

option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, help="path to all count", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="squire_normdeg",help="directory name to be create", metavar="character"),
    make_option(c("-p", "--processes"), type="integer", default=NULL,help="number of processes", metavar="integer"),
    make_option(c("-r", "--sample_rm"), type="character", default=NULL,help="sampèles to be removed expect a text files with Sample_ID in first column", metavar="integer"),
    make_option(c("-d", "--deg_design"), type="character", default=NULL,help="path to design sample configuration must have: \n
                                                                                Sample_ID,Sample_Condition,Sample_Replicate\n
                                                                                Sample_Condition must be levelled conditions separated by an undescore :_", metavar="integer"),
    make_option(c("-c", "--revel_conditions"), type="character", default=NULL,help="comma separated list of Conditions to be set in Deseq2 design revels", metavar="integer")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}else{
  ## Within the "count" folder we expect 5 files per sample:
  # ".out_abund.txt"
  # ".out.gtf"
  # ".out_refGenecounts.txt"
  # ".out_subFcounts.txt"
  # ".out_TEcounts.txt"
  abund = list.files(opt$file,recursive=TRUE,full.names=TRUE,pattern="abund.txt")
  refGenecounts = list.files(opt$file,recursive=TRUE,full.names=TRUE,pattern="refGenecounts.txt")
  subFcounts = list.files(opt$file,recursive=TRUE,full.names=TRUE,pattern="subFcounts.txt")
  TEcounts = list.files(opt$file,recursive=TRUE,full.names=TRUE,pattern="TEcounts.txt")
}
if( length(abund) == 0 |  length(refGenecounts) == 0 |  length(subFcounts) == 0 |  length(TEcounts) == 0 ){
  stop("Check count folder as some files are missing", call.=FALSE)
}

cpus=opt$processes
if(!is.null(cpus)){
  param <- BiocParallel::SnowParam(workers = cpus,tasks=0,stop.on.error=TRUE,progressbar=TRUE, type = "SOCK")
}else{
  param <- NULL
}

save_folder = paste0(getwd(),"/",opt$out,"/")
dir.create(save_folder)

# 1) 
# Gene ID	Gene Name	Reference	Strand	Start	End	Coverage	FPKM	TPM

abund_tpm = lapply(1:length(abund),function(x){
  y = read.delim(abund[x],sep="\t",row.names=NULL,stringsAsFactors=FALSE)
  rownames(y) = paste0(y$Gene.ID,"_",y$Reference,"_",y$Start,"_",y$End,"_",y$Strand)
  z = cbind(rownames(y),y$TPM)
  z = as.data.frame(z,stringsAsFactors=FALSE)
  colnames(z) = c("id",gsub("\\..*","",basename(abund[x])))
  return(z)
})
abund_tpm_m = abund_tpm %>% Reduce( function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="id"), .)
rownames(abund_tpm_m) = abund_tpm_m$id
abund_tpm_m = abund_tpm_m[,!grepl("id",colnames(abund_tpm_m))]
abund_tpm_m[is.na(abund_tpm_m)] = 0

# 1) 
# StringTie parameters for gene expression counting (-M .95) for guided assembly

header = c( "chr", #chromosome of transcription
            "tx_start" ,# coordinate for start of transcription
            "tx_stop" ,# coordinate for end of transcription
            "Gene_ID" ,# gene name
            "fpkm" , # fragments per kilobase per million reads
            "strand" ,# + or - for stranded data
            "count" ,# read count, computed by transcript length * coverage / readlength # TPM!
            "transcript_ID" ) #transcript specifier

refGenecounts_count = lapply(1:length(refGenecounts),function(x){
  y = read.delim(refGenecounts[x],sep="\t",row.names=NULL,stringsAsFactors=FALSE,header=FALSE)
  colnames(y) = header
  rownames(y) = paste0(y$Gene_ID,"_",y$chr,"_",y$tx_start,"_",y$tx_stop,"_",y$strand)
  z = cbind(rownames(y),y$count)
  z = as.data.frame(z,stringsAsFactors=FALSE)
  colnames(z) = c("id",gsub("\\..*","",basename(refGenecounts[x])))
  return(z)
})
refGenecounts_count = refGenecounts_count %>% Reduce( function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="id"), .)
rownames(refGenecounts_count) = refGenecounts_count$id
refGenecounts_count = refGenecounts_count[,!grepl("id",colnames(refGenecounts_count))]
refGenecounts_count[is.na(refGenecounts_count)] = 0

# 2)
# Sample	aligned_libsize	Subfamily:Family:Class	copies	fpkm	uniq_counts	tot_counts	tot_reads	score

subFcounts_tot_reads = lapply(1:length(subFcounts),function(x){
  y = read.delim(subFcounts[x],sep="\t",row.names=NULL,stringsAsFactors=FALSE)
  rownames(y) = y$Subfamily.Family.Class
  z = cbind(rownames(y),y$tot_reads)
  z = as.data.frame(z,stringsAsFactors=FALSE)
  colnames(z) = c("id",gsub("\\..*","",basename(subFcounts[x])))
  return(z)
})
subFcounts_tot_reads = subFcounts_tot_reads %>% Reduce( function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="id"), .)
rownames(subFcounts_tot_reads) = subFcounts_tot_reads$id
subFcounts_tot_reads = subFcounts_tot_reads[,!grepl("id",colnames(subFcounts_tot_reads))]
subFcounts_tot_reads[is.na(subFcounts_tot_reads)] = 0

subFcounts_tot_counts = lapply(1:length(subFcounts),function(x){
  y = read.delim(subFcounts[x],sep="\t",row.names=NULL,stringsAsFactors=FALSE)
  rownames(y) = y$Subfamily.Family.Class
  z = cbind(rownames(y),y$tot_counts)
  z = as.data.frame(z,stringsAsFactors=FALSE)
  colnames(z) = c("id",gsub("\\..*","",basename(subFcounts[x])))
  return(z)
})

subFcounts_tot_counts = subFcounts_tot_counts %>% Reduce( function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="id"), .)
rownames(subFcounts_tot_counts) = subFcounts_tot_counts$id
subFcounts_tot_counts = subFcounts_tot_counts[,!grepl("id",colnames(subFcounts_tot_counts))]
subFcounts_tot_counts[is.na(subFcounts_tot_counts)] = 0

# 3) 
# Headers are:
# tx_chr	tx_start	tx_stop	TE_ID	fpkm	tx_strand	Sample	alignedsize	TE_chr	TE_start	TE_stop	TE_name	milliDiv	TE_strand	uniq_counts	tot_counts	tot_reads	score
# uniq_counts	: ..
# tot_counts : fractional counts assigned based on the Count algorithm 
# tot_reads : sum of unique and multi-mapping reads that align to that TE

TEcounts_tot_reads = lapply(1:length(TEcounts),function(x){
  y = read.delim(TEcounts[x],sep="\t",row.names=NULL,stringsAsFactors=FALSE)
  rownames(y) = y$TE_ID
  z = cbind(rownames(y),y$tot_reads)
  z = as.data.frame(z,stringsAsFactors=FALSE)
  colnames(z) = c("id",gsub("\\..*","",basename(subFcounts[x])))
  return(z)
})
TEcounts_tot_reads = TEcounts_tot_reads %>% Reduce( function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="id"), .)
rownames(TEcounts_tot_reads) = TEcounts_tot_reads$id
TEcounts_tot_reads = TEcounts_tot_reads[,!grepl("id",colnames(TEcounts_tot_reads))]
TEcounts_tot_reads[is.na(TEcounts_tot_reads)] = 0
TEcounts_tot_reads <- dplyr::mutate_all(TEcounts_tot_reads, function(x) as.numeric(as.character(x)))
TEcounts_tot_reads = TEcounts_tot_reads[rowSums(TEcounts_tot_reads)!=0,]

TEcounts_tot_counts = lapply(1:length(TEcounts),function(x){
  y = read.delim(TEcounts[x],sep="\t",row.names=NULL,stringsAsFactors=FALSE)
  rownames(y) = y$TE_ID
  z = cbind(rownames(y),y$tot_counts)
  z = as.data.frame(z,stringsAsFactors=FALSE)
  colnames(z) = c("id",gsub("\\..*","",basename(subFcounts[x])))
  return(z)
})
TEcounts_tot_counts = TEcounts_tot_counts %>% Reduce( function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="id"), .)
rownames(TEcounts_tot_counts) = TEcounts_tot_counts$id
TEcounts_tot_counts = TEcounts_tot_counts[,!grepl("id",colnames(TEcounts_tot_counts))]
TEcounts_tot_counts[is.na(TEcounts_tot_counts)] = 0
TEcounts_tot_counts <- dplyr::mutate_all(TEcounts_tot_counts, function(x) as.numeric(as.character(x)))
TEcounts_tot_counts = TEcounts_tot_counts[rowSums(TEcounts_tot_counts)!=0,]

### APPLY Skewed normalisation -
### Compute scaling factors:
## Re-construct the design based on Sample Names which needs to be separated by "_" for each conditional level
## Replicates will begin with "^R" followed by a number [0-9]

deg_design = read.delim(opt$deg_design,sep="\t",header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
rownames(deg_design) = deg_design$Sample_ID

sample_rm=opt$sample_rm
if(!is.null(sample_rm)){
  rmS = read.delim(sample_rm,sep="\t",header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
  rownames(rmS) = rmS$Sample_ID
  deg_design = deg_design[deg_design$Sample_ID %!in% rownames(rmS),]
}

dat_ll =list( "abund_tpm_m" = abund_tpm_m[,deg_design$Sample_ID],
              "refGenecounts_count" = refGenecounts_count[,deg_design$Sample_ID],
              "subFcounts_tot_reads" = subFcounts_tot_reads[,deg_design$Sample_ID],
              "subFcounts_tot_counts" = subFcounts_tot_counts[,deg_design$Sample_ID],
              "TEcounts_tot_reads" = TEcounts_tot_reads[,deg_design$Sample_ID],
              "TEcounts_tot_counts" = TEcounts_tot_counts[,deg_design$Sample_ID]
)

reVel = unlist(strsplit(opt$revel_conditions,","))
r = lapply(reVel,function(x){
    g = unique(deg_design$Sample_Condition)[grep(x,unique(deg_design$Sample_Condition))]
    return(as.character(g))
  })
r = Reduce(intersect,r)
if(length(r)==0){stop("!! elements in revel are not contained within levels of Sample_Condition")}

### Execute Normalisation on each:
for( j in 1:length(dat_ll)){
  cat(names(dat_ll)[j],"\n")
  x = dat_ll[[j]]
  
  save_folder_sub = paste0(save_folder,"/",names(dat_ll)[j],"_InternalNorm/")
  dir.create(save_folder_sub)

  Result = RunNorm(x,deg_design,fix_reference="random",row_name_index=1,saving_path=save_folder_sub,n_pop=1,n_pop_reference=1,BiocParam=param)
  QC_plot(Result,deg_design,saving_path=save_folder_sub)

  # Perform DEGs:
  # This will design on the assembled design table:

  colData = Result$scaling_factors
  mat = Result$norm_mat
  
  AV = as.data.frame(matrix(NA,ncol=length(unique(colData$Sample_Condition)),nrow= nrow(mat) ),stringsAsFactors=FALSE)
  colnames(AV) = unique(colData$Sample_Condition)
  rownames(AV) = rownames(mat)
  for(cc in unique(colData$Sample_Condition)){
      sel = rownames(colData)[colData$Sample_Condition %in% cc]
      w1 = which(colnames(mat) %in% sel)
      AV[,cc] = rowMeans(mat[,w1])
  }
  AV = round(AV)

  mat = round(Result$norm_mat,0)
  rn = lapply(split_tibble(colData,"Sample_Condition"),function(x){
      return( rows_in_lim( mat[,rownames(x)],5,max(AV)) )
  })
  rn =unique(unlist(rn))
  mat = round(mat[rn,rownames(colData)],0)

  cat("Run DESEQ2:\n")
  formula = ~  Sample_Replicate + Sample_Condition   
  dds_ex <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = colData, design = formula)
  dds_ex <- DESeq2::estimateSizeFactors(dds_ex)
  SummarizedExperiment::colData(dds_ex)$sizeFactor <- 1
  dds_ex$Sample_Condition = relevel(   dds_ex$Sample_Condition, r)
  dds_ex <- DESeq2::estimateDispersions(dds_ex,fitType="parametric",maxit=100000) # parametric local mean
  dds_ex <- DESeq2::nbinomWaldTest(dds_ex,maxit = 100000)       
  
  pdf(paste0(save_folder_sub,"/MAplot.pdf"))
        DESeq2::plotDispEsts(dds_ex)
  dev.off()

  contrasts = as.data.frame(t(combn(unique(colData$Sample_Condition),2)))

  res_l = list( Norm_average_expression=AV )

  save_folder_sub_stat = paste0(save_folder_sub,"/STAT/")
  dir.create(save_folder_sub_stat)

  cat("By contrasts:\n")
  
  for(rn in 1:nrow(contrasts)){

      res = DESeq2::results(dds_ex, contrast=c("Sample_Condition",contrasts[rn,1],contrasts[rn,2]))
      nm <- paste0(save_folder_sub_stat,"/",contrasts[rn,1],"_vs_",contrasts[rn,2],"_RESULTS.txt")
      pdf(paste0(save_folder_sub_stat,"/",contrasts[rn,1],"_vs_",contrasts[rn,2],"_RESULTS.pdf"))
              DESeq2::plotMA(res,alpha = 0.01, main = "", xlab = "mean of normalized counts", MLE = FALSE)
      dev.off()
      
      res = as.data.frame(res,stringsAsFactors=FALSE)

      cn_samp = colData[grep(paste0(contrasts[rn,2],"|",contrasts[rn,1]),colData$Sample_Condition), ]
      AV_sel = AV[,as.character(unique(cn_samp$Sample_Condition))]
      
      res = cbind(AV_sel[rownames(AV_sel),],res[rownames(AV_sel),])

      res_l[[ paste0(contrasts[rn,1],"_vs_",contrasts[rn,2] )  ]] = res

      write.table(res,paste0(save_folder_sub_stat,"/",contrasts[rn,1],"_vs_",contrasts[rn,2],"_DESEQ2result.txt"),sep="\t",col.names=NA)

  }
  
  # cat("Save results to xlsx")
  # output_file <- paste0(save_folder_sub_stat,"/MAT_Results.xlsx")
  # write.xlsx( res_l, output_file ,row.names = TRUE )

}

 