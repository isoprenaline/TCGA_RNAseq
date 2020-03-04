### ---------------
rm(list=ls())
options(stringsAsFactors = F)
# 注意，并不是说使用 RTCGA.miRNASeq包的数据是最佳选择，只是因为这个演示起来最方便。
# 因为GDC官网下载数据具有一定门槛，也不是每个人都必须学会的。
getwd()
Rdata_dir='./Rdata/'
Figure_dir='./figures/'

# 如果开启下面代码，就会从RTCGA.rnaSeq包里面提取mRNA表达矩阵和对应的样本临床信息。 
if(T){
  library(RTCGA.rnaseq)
  s=GBM.rnaseq[,1]
  expr <- expressionsTCGA(GBM.rnaseq)
  dim(expr)
  expr[1:40,1:4]
  mi=colnames(expr)
  expr=apply(expr,1,as.numeric) 
  colnames(expr)=s
  rownames(expr)=mi
  expr[1:4,1:4]
  expr=na.omit(expr)
  expr=expr[apply(expr, 1,function(x){sum(x>1)>10}),]
  rownames(expr) <- gsub('\\|','_',rownames(expr))
  
  library(RTCGA.clinical) 
  meta <- GBM.clinical
  tmp=as.data.frame(colnames(meta))
  meta[(grepl('patient.bcr_patient_barcode',colnames(meta)))]
  meta[(grepl('patient.days_to_last_followup',colnames(meta)))]
  meta[(grepl('patient.days_to_death',colnames(meta)))]
  meta[(grepl('patient.vital_status',colnames(meta)))]
  ## patient.race  # patient.age_at_initial_pathologic_diagnosis # patient.gender 
  # patient.stage_event.clinical_stage
  meta=as.data.frame(meta[c('patient.bcr_patient_barcode','patient.vital_status',
                            'patient.days_to_death','patient.days_to_last_followup',
                            'patient.race',
                            'patient.age_at_initial_pathologic_diagnosis',
                            'patient.gender' ,
                            'patient.stage_event.pathologic_stage')])
  #meta[(grepl('patient.stage_event.pathologic_stage',colnames(meta)))]
  ## 每次运行代码，就会重新生成文件。
  save(expr,meta,file = './Rdata/TCGA-GBM-rnaseq-example.Rdata')
}

## 我们已经运行了上面被关闭的代码，而且保存了miRNA表达矩阵和对应的样本临床信息
# 现在直接加载即可。
load(file= './Rdata/TCGA-GBM-rnaseq-example.Rdata')

dim(expr)
dim(meta)
# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

# 这里需要解析TCGA数据库的ID规律，来判断样本归类问题。
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
table(group_list)
exprSet=na.omit(expr)

# 取出Normal样品，只留下Tumor样品
colData <- data.frame(row.names=colnames(exprSet), 
                      group_list=group_list)
expr_tumor = expr[,as.numeric(substr(colnames(expr),14,15)) < 10]
expr_tumor <- expr_tumor[-1:-21,]
expr_tumor[1:4,1:4]
rownames(expr_tumor)[which(rownames(expr_tumor)=="SLC35E2_9906")] <- "SLC35E2A_9906"
rownames(expr_tumor)[which(rownames(expr_tumor)=="SLC35E2_728661")] <- "SLC35E2B_728661"
expr_tumor["SLC35E2A_9906",]
expr_tumor["SLC35E2B_728661",]

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
  breat_gtex=GTEx[,gsub('[.]','-',colnames(GTEx)) %in% b[b$SMTS=='Pancreas',1]]
  rownames(breat_gtex)=GTEx[,1]
  dat=breat_gtex
  
  #把基因ID转换为symbol
  dat=breat_gtex
  ids=GTEx[,1:2]
  head(ids)
  colnames(ids)=c('probe_id','symbol')
  dat=dat[ids$probe_id,]
  dat[1:4,1:4] 
  ids$median=apply(dat,1,median)
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  ids=ids[!duplicated(ids$symbol),]
  dat=dat[ids$probe_id,]
  rownames(dat)=ids$symbol
  dat[1:4,1:4] 
  breat_gtex=dat
  save(breat_gtex,file = './GTEx/Pancreas_gtex_counts.Rdata')
}

#把正常样品与肿瘤样品数据合并
if(F){
  load(file= './GTEx/Pancreas_gtex_counts.Rdata')
  normal_data <- breat_gtex[unlist(lapply(strsplit(rownames(expr_tumor),split = '_'), function(x){which(rownames(breat_gtex)==x[[1]][1])})),]
  normal_data <- cbind(rownames(normal_data),normal_data)
  names(normal_data)[1] <- "symbol"
  tumor_data <- as.data.frame(cbind(unlist(lapply(strsplit(rownames(expr_tumor),split = '_'), function(x){x[[1]][1]})),expr_tumor))
  names(tumor_data)[1] <- "symbol"
  rownames(tumor_data) <- tumor_data$symbol
  data_for_limma <- merge(normal_data,tumor_data,by.x = "symbol")
  rownames(data_for_limma) <- data_for_limma[,'symbol']
  data_for_limma <- as.matrix(data_for_limma[,-1])
  
  #分组
  cnames <- substr(colnames(data_for_limma),1,4)
  group_list <- ifelse(apply(as.data.frame(cnames), 1, function(x){ifelse(x=='GTEX',1,0)}),'normal','tumor')
  table(group_list)
  exprSet=na.omit(data_for_limma)
  class(exprSet) <- "numeric"
  exprSet[1:4,1:4]
  source('./functions.R')
}

#筛选出三阴性乳腺癌
if(F){
  patient_brca <- read.table('c:/Users/jufei/Documents/TCGA/BRCA/nationwidechildrens.org_clinical_patient_brca.txt',sep = '\t',fill = T,header = T)[-c(1,2),]
  patient_brca <- patient_brca[,c('bcr_patient_uuid','bcr_patient_barcode','er_status_by_ihc','pr_status_by_ihc','her2_status_by_ihc')]
  TNBC_patient_brca <- patient_brca[apply(patient_brca, 1, function(x){x['er_status_by_ihc']=='Negative'&x['pr_status_by_ihc']=='Negative'&x['her2_status_by_ihc']=='Negative'}),]
  expr_tumor = expr_tumor[,unlist(lapply(TNBC_patient_brca$bcr_patient_barcode, function(x){grep(x,colnames(expr_tumor),fixed = T)}))]
  write.csv(TNBC_patient_brca,'TNBC_patient_brca.csv')
  write.csv(expr_tumor,'TNBC_row_data.csv')}

# 判断是否是altered样品
if(F){
  cnames <- substr(colnames(expr_tumor),1,15)
  altered <- read.table(file = 'c:/Users/jufei/Documents/TCGA/GBM/GBM_altered_samples.txt',sep = ':')[,2]
  group_list <- ifelse(apply(as.data.frame(cnames), 1, function(x){ifelse(sum(grepl(x,altered)),1,0)}),'altered','unaltered')

  table(group_list)
  exprSet=na.omit(expr_tumor)
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
  cont.matrix=makeContrasts(contrasts=c('altered-unaltered'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='altered-unaltered', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom)
  nrDEG=DEG_limma_voom[,c(1,4)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
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
  exons_gene_lens_1 <- do.call(rbind,exons_gene_lens)
  exons_gene_lens_1 <- cbind(do.call(rbind,lapply(strsplit(rownames(exons_gene_lens_1),split = '[.]'), function(x){x[[1]][1]})),exons_gene_lens_1)
  colnames(exons_gene_lens_1) <- c("ensembl_gene_id","length")
  rownames(exons_gene_lens_1) <- exons_gene_lens_1[,1]
  Human_exons_gene <- getBM(attributes =c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id", 
                       values = do.call(rbind,lapply(strsplit(rownames(exons_gene_lens_1),split = '[.]'), function(x){x[[1]][1]})),
                       mart = mart_h)
  Human_exons_gene_lens <- merge(exons_gene_lens_1,Human_exons_gene,by="ensembl_gene_id")
  write.csv(Human_exons_gene_lens,"Human_exons_gene_lens.csv")
  #TPM计算
  kb <- Human_exons_gene_lens$length
  class(kb) <- "numeric"
  kb <- kb/1000
  kb
  countdata <- exprSet[,1:9]
  rownames <- as.data.frame(Human_exons_gene_lens[unlist(lapply(rownames(exprSet), function(x){which(Human_exons_gene_lens$external_gene_name == x)})),1])
  rpk <- countdata / kb
  rownames <- lapply(rownames(exprSet), function(x){which(Human_exons_gene_lens$external_gene_name == x)})
  rpk
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  head(tpm)
  write.table(tpm,file="2020武汉加油_tpm.xls",sep="\t",quote=F)
}

#作图
if(F){
  #function.R火山图标注基因名
  nrDEG$sign <- ifelse(nrDEG$pvalue < 0.01 & abs(nrDEG$log2FoldChange) > 1,rownames(nrDEG),NA)
  draw_h_v(exprSet,nrDEG,'limma',group_list,1)
  write.csv(nrDEG[,-3],'GBM_altered_vs_unaltered.csv')
  nrDEG_cut <- nrDEG[abs(nrDEG["log2FoldChange"])>=0.5 & nrDEG["pvalue"]<=0.001,-3]
  write.csv(nrDEG_cut,'GBM_altered_vs_unaltered_cut.csv')
  
  #计算相关性矩阵
  library(Hmisc)
  set <- rownames(nrDEG_cut)
  data <- t(expr[apply(as.data.frame(set), 1, function(x){grep(x,rownames(expr),fixed = T)}),])
  
  #把PDCD1，CTLA4，LAG3和TIGIT的表达数据添加进去
  data <- cbind(data, t(expr[c("PDCD1_5133", "CTLA4_1493","TIGIT_201633", "LAG3_3902"),]))
  
  res2 <- rcorr(as.matrix(data),type = 'spearman')
  res2r <- res2$r
  res2r4i <- res2r[,c("PDCD1_5133", "CTLA4_1493", "TIGIT_201633", "LAG3_3902")]
  res2r4i <- res2r4i[!duplicated(res2r4i),]
  res2r4i_cut4 <- res2r4i[apply(res2r4i,1,function(x) {ifelse(sum(x>=0.4)==4,T,F)}),]
  write.csv(res2r4i,'GBM_res2r4i.csv')
  write.csv(res2r4i_cut4,'GBM_res2r4i_cut4.csv')
  
  #绘制相关性散点图
  library(ggpubr)
  data <- as.data.frame(data)
  data_r <- data[,unlist(apply(as.data.frame(rownames(res2r4i_cut4)), 1, function(x){grep(x,colnames(data),fixed = T)}))]
  data_r <- cbind(data_r,data[,c("CTLA4_1493","TIGIT_201633", "LAG3_3902")])
  
  #向量方式输出y轴变量
  paste0('c(',paste(paste0('\'',colnames(data_r)[-match(c("PDCD1_5133", "CTLA4_1493","TIGIT_201633", "LAG3_3902"),colnames(data_r))],'\''),collapse = ','),')')
  
  ggscatter(data_r, x = "LAG3_3902", 
            y = c('SIGLEC1_6614','PTPRCAP_5790','MX2_4600','ADAM6_8755','IL2RG_3561','CD27_939','CD3E_916','SAMD3_154075','CD8A_925','IGJ_3512','SLA2_84174','ANO9_338440','MGC29506_51237','PDCD1_5133.1','HSH2D_84941','CD247_919','SIRPG_55423','AIM2_9447','ZBP1_81030','SLAMF7_57823','ODF3B_440836','KLHDC7B_113730'), 
            conf.int = TRUE, combine = TRUE,
            color = "red", shape = 21, size = 1, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
            cor.coef = TRUE, cor.method = "spearman",cor.coef.size =2,
            xlab = "LAG3_3902", ylab = "Genes", 
            xscale = "log10", yscale="log10")
  
  ggscatter(data_r, x = "CD300A_11314", y = c("PDCD1_5133", "CTLA4_1493", "TIGIT_201633", "LAG3_3902"), 
            conf.int = TRUE, combine = TRUE,
            color = "red", shape = 21, size = 1, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
            cor.coef = TRUE, cor.method = "spearman",cor.coef.size =2,
            xlab = "CD300A", ylab = "Genes", 
            xscale = "log10", yscale="log10")
  
}

#基因注释
if(F){
  ##Gene description
  library(biomaRt)
  mart_h <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  Description <- getBM(attributes =c("ensembl_gene_id", "external_gene_name","description"),
          filters = "external_gene_name", values = rownames(exprSet),
          mart = mart_h)
  FC <- nrDEG
  FC <- cbind(FC,rownames(FC))
  names(FC)[3] <- "external_gene_name"
  Descriptioned <- merge(FC,Description,by="external_gene_name")
  write.csv(Descriptioned,"GBM_Descriptioned_limma_DEG_Tumor_vs_Normal_cut.csv")
}


