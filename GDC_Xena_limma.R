rm(list=ls())
options(stringsAsFactors = F)
dir_Rdata <- 'c:/Users/jufei/Documents/TCGA/Rdata/'
dir_Results <- 'c:/Users/jufei/Documents/TCGA/Results/'

#提取GDC的临床数据
if(F){
  library("XML")
  library("methods")
  dir_XML='c:/Users/jufei/Documents/TCGA/LIHC/XML/'
  
  all_fiels=list.files(path = dir_XML ,pattern='*.xml$',recursive=T)
  all_fiels
  cl = lapply(all_fiels
              , function(x){
                #x=all_fiels[1]
                result <- xmlParse(file = file.path(dir,x)) 
                rootnode <- xmlRoot(result)  
                xmldataframe <- xmlToDataFrame( rootnode[2] ) 
                return(t(xmldataframe))
              })
  
  cl_df <- t(do.call(cbind,cl))
  save(cl_df,file = 'GDC_TCGA_LUAD_clinical_df.Rdata')
}
  
#提取GDC的TCGA数据并整理
if(F){
  dir_GDC='c:/Users/jufei/Documents/TCGA/LIHC/LIHC_TCGA'
  mi = lapply(list.files(path = dir_GDC ,pattern='*.htseq.counts.gz$',recursive=T)
              , function(x){
                result <- read.table(file = file.path(dir,x),sep = '\t',header = F)[,1:2]
                return( result )
              })
  
  mi_df <- t(do.call(cbind,mi))
  dim(mi_df)
  mi_df[1:4,1:4]
  colnames(mi_df)=mi_df[1,]
  mi_df=mi_df[seq(2,nrow(mi_df),by=2),]
  mi_df[1:4,1:4]
  
  library(rjson)
  json <- fromJSON(file = 'c:/Users/jufei/Documents/TCGA/LIHC/LIHC_files.2020-02-29.json')
  caseid <- as.data.frame(toupper(unlist(lapply(json,function(x){x$cases[[1]]$case_id}))))
  filename <- as.data.frame(unlist(lapply(json,function(x){x["file_name"]})))
  case_id <- cbind(caseid,filename)
  colnames(case_id) <- c('bcr_patient_uuid','filenames')
  
  filenames <- unlist(lapply(list.files(path = dir ,pattern='*.htseq.counts.gz$',recursive=T), function(x){return(strsplit(x, split = "/")[[1]][2])}))
  mi_df <- cbind(filenames,mi_df)
  mi_df[1:4,1:4]
  patient <- cl_df[,4:5]
  patient[,2] <- toupper(patient[,2])
  df_case_id <- merge(case_id,patient,all.x=T)
  expr <- t(merge(df_case_id,mi_df))
  colnames(expr) <- expr["bcr_patient_barcode",]
  expr[1:4,1:4]
  expr <- expr[-1:-3,]
  dim(expr)
  
  save(expr,cl_df,file = "LIHC_GDC_TCGA_rowdata.Rdata")
}

#提取xena数据并去除非肿瘤组织样本
if(F){
  mi_df <- read.table(file = 'c:/Users/jufei/Documents/TCGA/LIHC/TCGA-LIHC.htseq_counts.tsv.gz',sep = '\t',header = T,row.names = 1)
  mi_df <- (2^mi_df)-1
  names(mi_df) <- gsub('[.]','_',colnames(mi_df))
  save(mi_df,file = "LIHC_xena_TCGA_count_rowdata.Rdata")
  expr <- mi_df[,ifelse(as.numeric(substr(colnames(mi_df),14,15)) < 10,T,F)]
}

#提取xena生存期数据
if(F){
  survival_data <- read.table(file = "c:/Users/jufei/Documents/TCGA/LIHC/TCGA-LIHC.survival.tsv.gz",sep = '\t',header = T)
  names(survival_data) <- gsub('[.]','_',colnames(survival_data))
  survival_data$sample <- gsub('-','_',survival_data$sample)
}
  
load(file= 'LIHC_GDC_TCGA_rowdata.Rdata')

dim(expr)

#从GTEx数据中提取需要分析的正常样本
if(F){
  options(stringsAsFactors = F)
  GTEx=read.table('~/TCGA/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
                  ,header = T,sep = '\t',skip = 2)
  # SMTS Tissue Type, area from which the tissue sample was taken. 
  # SMTSD Tissue Type, more specific detail of tissue type
  b=read.table('~/TCGA/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
               header = T,sep = '\t',quote = '')
  table(b$SMTS) 
  breat_gtex=GTEx[,gsub('[.]','-',colnames(GTEx)) %in% b[b$SMTS=='Liver',1]]
  rownames(breat_gtex)=GTEx[,1]
  dim(breat_gtex)
  
  save(breat_gtex,file = 'Liver_gtex_geneID_counts.Rdata')
}

#把正常样品与肿瘤样品数据合并
if(F){
  load(file= 'c:/Users/jufei/Documents/TCGA/GTEx/Liver_gtex_geneID_counts.Rdata')
  load(file= 'c:/Users/jufei/Documents/TCGA/LIHC/Rdata/LIHC_GDC_TCGA_rowdata.Rdata')
  normal_data <- breat_gtex
  normal_data <- normal_data[-(grep("PAR_Y",rownames(normal_data))),]
  rownames(normal_data) <- do.call(rbind,lapply(strsplit(rownames(normal_data),split = '[.]'),function(x)return(x[[1]][1])))
  normal_data <- cbind(rownames(normal_data),normal_data)
  names(normal_data)[1] <- "Gene_ID"
  normal_data[1:4,1:4]
  tumor_data <- expr
  rownames(tumor_data) <- do.call(rbind,lapply(strsplit(rownames(tumor_data),split = '[.]'),function(x)return(x[[1]][1])))
  tumor_data <- as.data.frame(cbind(rownames(tumor_data),tumor_data))
  names(tumor_data)[1] <- "Gene_ID"
  tumor_data[1:4,1:4]
  data_for_limma <- merge(tumor_data,normal_data,by = "Gene_ID")
  
  rownames(data_for_limma) <- data_for_limma$Gene_ID
  data_for_limma <- as.matrix(data_for_limma[,-1])
  data_for_limma[1:4,1:4]
  
  #分组
  cnames <- substr(colnames(data_for_limma),1,4)
  group_list <- ifelse(apply(as.data.frame(cnames), 1, function(x){ifelse(x=='GTEX',1,0)}),'normal','tumor')
  table(group_list)
  exprSet=na.omit(data_for_limma)
  class(exprSet) <- "numeric"
  exprSet[1:4,1:4]
  save(exprSet,file = "LIHC_TCGA_vs_Liver_GTEx_rowdata.Rdata")
  
  load(file = "./LIHC/LIHC_TCGA_vs_Brain_GTEx_rowdata.Rdata")
  
  rm(breat_gtex,expr,cnames)
  source('./functions.R')
}

#差异分析
if(T){
  suppressMessages(library(limma))
  library(edgeR)
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('tumor-normal'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='tumor-normal', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom)
  nrDEG=DEG_limma_voom[,c(1,4)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
}

#作图
if(F){
  #function.R火山图标注基因名
  nrDEG$sign <- ifelse(nrDEG$pvalue < 0.0000001 & abs(nrDEG$log2FoldChange) > 5,rownames(nrDEG),NA)
  draw_h_v(exprSet,nrDEG,'limma',group_list,1)
  write.csv(nrDEG[,-3],'LIHC_TCGA_vs_GTEx_Total.csv')
  nrDEG_cut <- nrDEG[abs(nrDEG["log2FoldChange"])>=0.5 & nrDEG["pvalue"]<=0.001,-3]
  write.csv(nrDEG_cut,'LIHC_TCGA_vs_GTEx_Cut.csv')
}

#基因注释
if(F){
  ##Gene description
  library(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  Description <- getBM(attributes =c("ensembl_gene_id", "external_gene_name","description"),
                       filters = "ensembl_gene_id", values = rownames(nrDEG_cut),
                       mart = mart)
  FC <- nrDEG_cut
  FC <- cbind(FC,rownames(FC))
  names(FC)[3] <- "ensembl_gene_id"
  Descriptioned <- merge(FC,Description,by="ensembl_gene_id")
  write.csv(Descriptioned,"LIHC_Descriptioned_TCGA_vs_GTEx_FC1.5_P0.0001.csv",row.names = F)
}

#计算log2(TPM)
if(F){
  #获取每个基因去除冗余序列后的外显子长度
  suppressMessages(library(GenomicFeatures))
  library(parallel)#并行计算
  txdb <- makeTxDbFromGFF("gencode.v33.annotation.gtf.gz",format="gtf")
  exons_gene <- exonsBy(txdb, by = "gene")
  gc()
  cl <- makeCluster(getOption("cl.cores",6))
  exons_gene_lens <- parLapply(cl,exons_gene,function(x){sum(width(reduce(x)))})
  stopCluster(cl)
  exons_gene_lens <- do.call(rbind,exons_gene_lens)
  rownames(exons_gene_lens) <- do.call(rbind,strsplit(rownames(exons_gene_lens),split = '[.]'))[,1]
  names(exons_gene_lens)[1] <- "length"
  write.csv(Human_exons_gene_lens,"Human_exons_gene_lens.csv")
  save(exons_gene_lens,file = "Human_exons_gene_lens.Rdata")
  
  load(file = "Human_exons_gene_lens.Rdata")
  
  #TPM计算
  kb <- exons_gene_lens[intersect(rownames(exons_gene_lens),rownames(exprSet)),]/1000
  kb[1:4]
  countdata <- exprSet[intersect(rownames(exons_gene_lens),rownames(exprSet)),]
  rpk <- countdata / kb
  rpk[1:4,1:4]
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  save(tpm,file = "LIHC_TCGA_VS_GTEx_TPM.Rdata")
  
  load(file = paste0(dir_Rdata,'LIHC_TCGA_VS_GTEx_TPM.Rdata'))
  
  log2tpm <- log2(tpm+1)
  head(tpm)
}

#添加病人OS信息与注释
if(F){
  annotation_gene <- read.table(file = paste0(dir_Rdata,'gencode.v22.annotation.gene.probeMap'),sep = '\t',header = T)
  annotation_gene$id <- do.call(rbind,strsplit(as.character(annotation_gene$id),split = '[.]'))[,1]
  rownames(annotation_gene) <- annotation_gene$id
  temp <- as.data.frame(log2tpm)
  temp <- t(cbind(annotation_gene[rownames(temp),2],temp))
  rownames(temp)[1] <- 'Gene_name'
  temp <- as.data.frame(cbind(rownames(temp),temp))
  names(temp)[1] <- 'sample'
  temp <- merge(survival_data,temp,by='sample',all.y = T)
  temp[1:4,1:4]
  write.csv(temp,file = paste0(dir_Results,'LIHC_Log2TPM.csv'),row.names = F)
}
