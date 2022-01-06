# Troubleshooting -- Error; return code from  pthread_create() is 22 when using affy::justRMA()
# https://support.bioconductor.org/p/122925/#124701 suggests that this issue comes from an affy dependency called preprocessCore library.
# The following removes affy and preprocessCore libraries, reinstalls preprocessCore with threading disabled.
remove.packages("preprocessCore")
BiocManager::install("preprocessCore", configure.args="--disable-threading")

library(stringr)
library(Biobase)
library(GEOquery)
library(affy)
library(oligo)
library(affycoretools)
library(pd.hta.2.0)
library(RUVnormalize)
library(RUVnormalizeData)
library(AnnotationDbi)
install.packages(paste0(input_dir, "hta20hsensgcdf_24.0.0.tar.gz"), sep="", repos = NULL, type = "source")

input_dir <- "/pfs/downloadGBMData/"
cell_dir <- "/pfs/getGBMCellData/"
out_dir <- "/pfs/out/" 
functions <- "/pfs/getGBMCellData/functions.R"

# input_dir <- "~/Documents/pfs/downloadGBMData/"
# cell_dir <- "~/Documents/pfs/getGBMCellData/"
# out_dir <- "~/Documents/pfs/getGBMGeneExpression/" 
# functions <- "./functions.R"

source(functions)
load(paste0(input_dir, "Ensembl.v99.annotation.RData"))
cell <- readRDS(paste0(cell_dir, "cell.rds"))

# ======================== GSM_maping ========================
# Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160
# which includes GSM ids of the samples mapped to their corresponding cell ids

GSM_html<-readLines("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160")

GSM_map<- data.frame( accession_id= GSM_html[grep(">GSM", GSM_html)] , cellid= GSM_html[grep(">GSM", GSM_html) + 1])
GSM_map $accession_id <- str_sub( GSM_map $accession_id, -19,-10)
GSM_map $cellid <- str_sub(sub('<td valign=\"top\">' ,'', GSM_map $cellid) , end = -6)

# Creating replicate column including replicate markers (e.g. "A") 
GSM_map$Replicate <- gsub(".*cells", "", GSM_map$cellid)
GSM_map$Replicate [GSM_map$Replicate =="" | GSM_map$Replicate =="human astrocytes"]<-NA
GSM_map$Replicate <- gsub(" ", "", GSM_map$Replicate)


# Creating Patient-id column consistant with cell names in the main paper
GSM_map$Patient_id <-gsub("cells.*", "", GSM_map$cellid)
GSM_map$Patient_id <- gsub(" ", "", GSM_map$Patient_id)
GSM_map$Patient_id [GSM_map$Patient_id =="humanastrocytes"]<-"human_astrocytes"

# Creating unique cell-ids based on patient ids and replicates
GSM_map$cellid<-sub("MG cells", "", GSM_map$cellid)
GSM_map$cellid<-sub(" ", "_", GSM_map$cellid)

# GSM ids are only used in the expression data
GSM_map<-GSM_map[, c("Patient_id","Replicate","cellid", "accession_id")]

phen_exp<-merge(GSM_map , cell, by="Patient_id" , all.x=TRUE)
phen_exp<-phen_exp[, c(1:4 , 11:15)]
phen_exp$batchid <- NA
rownames(phen_exp)<-phen_exp$accession_id

# ============= expression data =============
print("expression data")
#Creating eset from raw expression data 
# GSE152160_RAW.tar file to be downloaded from here "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152160&format=file" 
untar(paste(input_dir, "GSE152160_RAW.tar", sep=""), exdir="GSE152160_RAW")#Unpack the CEL files
cels<-list.files("GSE152160_RAW/", pattern = "CEL")
for(cel in cels){
  if(!file.exists(paste0("GSE152160_RAW/", gsub("[.]gz$", "", cel)))){
    GEOquery::gunzip(filename=paste0("GSE152160_RAW/", cel), overwrite = TRUE, remove = FALSE)
  }
}
# sapply(paste("GSE152160_RAW", cels, sep="/"), GEOquery::gunzip)

####Gene-level expression
cdf <- "hta20hsensgcdf"
cels <- list.celfiles("GSE152160_RAW/", full.names = TRUE)#Raw_expression Folder includes 145 CEL files
expr_cel <- affy::justRMA(filenames = cels, verbose = TRUE, cdfname = cdf)

# Assay data expression 
assay_exp<-as.data.frame(exprs(expr_cel))
rownames(assay_exp)<-sub("_at","",rownames(assay_exp))
colnames(assay_exp)<-sub("_PA.*","",colnames(assay_exp))

# Removing control genes (which start with "AFFX")
assay_exp<-assay_exp[-c(grep("AFFX", rownames(assay_exp))),] 


# Feature data expression
feat_exp<-fdata_builder(annotation=features_gene, assay=assay_exp,ID_column="gene_id")#features_gene from Ensembel.v99.annotation.RData
feat_exp$BEST<-NA

# Protocol data expression normalized based on RMA only
protocol_exp<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp)),row.names=colnames(assay_exp))
colnames(protocol_exp)<-"Array"
protocol_exp$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp$Folder_raw_data<-"GSE152160_RAW.tar" #Folder including the raw data in the above link

#Creating ExpressionSet normalized based on RMA only
assay_exp<-assay_exp[,rownames(phen_exp)] #rearranging the colnames so it is similar to pheno data
protocol_exp<-protocol_exp[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet<- ExpressionSet(assayData = as.matrix(assay_exp), phenoData = AnnotatedDataFrame(phen_exp), 
                         featureData = AnnotatedDataFrame(feat_exp),
                         protocolData=AnnotatedDataFrame(protocol_exp)) 
expression_SE<- eSetToSE(exp_eSet,annot_name="rna")
print("expression data - normalized based on RMA only: done")

# ============= Normalizing based on negative control genes =============
ctrl<-read.delim(paste(input_dir, "HK_genes.txt", sep=""),stringsAsFactors = FALSE, header = FALSE) # Negative control genes 
ctrl$gene_name<-gsub(" ","",ctrl$V1) #Removing spaces from the gene names
ctrl<-merge(ctrl, feat_exp[,c("gene_name","gene_id")], by="gene_name")
ctrl_ind<-which(rownames(assay_exp) %in% ctrl$gene_id) 

## Prepare control samples
# A table that has as many columns as the largest set of replicates for one sample. Each
# row corresponds to a set of replicates of the same sample and gives the row indices of the
# replicates in the gene expression matrix, padded with -1 entries.
# See https://www.bioconductor.org/packages/release/bioc/vignettes/RUVnormalize/inst/doc/RUVnormalize.pdf

YY <- t(assay_exp) 
rep_cIdx <- matrix(-1,145,2) # 145 is number of cell lines and 2 is max number of replicates a cell has
added_pat<-c()
count <-1#row number of rep_cIdx

for (i in 1:nrow(rep_cIdx)){
  pat = phen_exp$Patient_id[i]
  if (pat %in% added_pat)
    next
  GSM=phen_exp$accession_id[phen_exp$Patient_id==pat]
  if(length(GSM)==2){
    rep_cIdx[count,1]=which(rownames(YY)==GSM[1])
    rep_cIdx[count,2]=which(rownames(YY)==GSM[2])
  }
  else(rep_cIdx[count,1]=which(rownames(YY)==GSM))
  
  added_pat=c(added_pat,pat)  
  count=count+1
}

## Replicate-based
# Set k to the number of samples / 4 or to the number of replicates, if the latter is smaller than the former. 
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896498/

Res <- naiveReplicateRUV(YY, ctrl_ind, rep_cIdx, k=17)

# Assay expression RUV
assay_exp_ruv<-data.frame(t(Res$cY))

# Protocol data expression RUV
protocol_exp_ruv<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp_ruv)),row.names=colnames(assay_exp_ruv))
colnames(protocol_exp_ruv)<-"Array"
protocol_exp_ruv$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp_ruv$Folder_raw_data<-"GSE152160_RAW.tar" #Folder including the raw data in the above link
protocol_exp_ruv$Negative_control_genes<-"https://www.tau.ac.il/~elieis/HKG/"

# Creating ExpressionSet RUV
assay_exp_ruv<-assay_exp_ruv[,rownames(phen_exp)] #rearranging the colnames so it is similar to pheno data
protocol_exp_ruv<-protocol_exp_ruv[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet_ruv<- ExpressionSet(assayData = as.matrix(assay_exp_ruv), phenoData = AnnotatedDataFrame(phen_exp), 
                             featureData = AnnotatedDataFrame(feat_exp),
                             protocolData=AnnotatedDataFrame(protocol_exp_ruv)) 
expression_ruv_SE<- eSetToSE(exp_eSet_ruv,annot_name="rna_ruv")
print("expression data - creating ExpressionSet RUV: done")

# ============= Probe-level expression =============
# memory.limit(size=76000)
raw_data<- oligo::read.celfiles(cels) # If required run this: memory.limit(size=76000) , memory.limit()
norm_data<-rma(raw_data,target="core") # Perform RMA normalization
norm_main<- getMainProbes(norm_data, level = "core")#Remove the control transcripts (Controls do not match to any genes)
norm_main_annot<- annotateEset(norm_main, pd.hta.2.0)#Annotating the eset on transcript level : "https://support.bioconductor.org/p/89308/"

# Assay data expression-probe level
assay_exp_probe<-as.data.frame(exprs(norm_main_annot))
colnames(assay_exp_probe)<-sub("_PA.*","",colnames(assay_exp_probe))

# Feature data expression-probe level
feat_exp_probe<-fData(norm_main_annot)
feat_exp_probe$BEST<-NA

# Pheno data expression-probe level is same as Pheno data expression

# Protocol data expression-probe level
protocol_exp_probe<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp_probe)),row.names=colnames(assay_exp_probe))
colnames(protocol_exp_probe)<-"Array"
protocol_exp_probe$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp_probe$annotation<-"pd.hta.2.0"

#Creating ExpressionSet
assay_exp_probe<-assay_exp_probe[,rownames(phen_exp)]#rearranging the colnames so it is similar to pheno data
protocol_exp_probe<-protocol_exp_probe[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet_probe<- ExpressionSet(assayData = as.matrix(assay_exp_probe), phenoData = AnnotatedDataFrame(phen_exp),
                               featureData = AnnotatedDataFrame(feat_exp_probe),
                               protocolData=AnnotatedDataFrame(protocol_exp_probe))
expression_probe_SE<- eSetToSE(exp_eSet_probe,annot_name="rna_probe")
print("expression data - probe level expression: done")

saveRDS(expression_SE, paste0(out_dir, "expression_SE.rds"))
saveRDS(expression_ruv_SE, paste0(out_dir, "expression_ruv_SE.rds"))
saveRDS(expression_probe_SE, paste0(out_dir, "expression_probe_SE.rds"))
saveRDS(phen_exp, paste0(out_dir, "phen_exp.rds"))