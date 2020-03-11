rm(list=ls())
options(stringsAsFactors = F)
source('./functions.R')

library(openxlsx)
wb <- createWorkbook()
my_style <- createStyle(fontName = '等线')
addWorksheet(wb, 'new sheet')

dir_Rdata <- 'c:/Users/jufei/Documents/TCGA/Rdata/'
dir_Results <- 'c:/Users/jufei/Documents/TCGA/Results/Altered_Correlation_Analysis/'

#提取xena数据并去除非肿瘤组织样本
if(F){
  mi_df <- read.table(file = 'c:/Users/jufei/Documents/TCGA/LIHC/TCGA-LIHC.htseq_counts.tsv.gz',sep = '\t',header = T,row.names = 1)
  mi_df <- mi_df[-(grep('_',rownames(mi_df))),]
  mi_df <- (2^mi_df)-1
  names(mi_df) <- gsub('[.]','_',colnames(mi_df))
  
  save(mi_df,file = paste0(dir_Rdata,"LIHC_xena_TCGA_count_rowdata.Rdata"))
  
  expr <- mi_df[,ifelse(as.numeric(substr(colnames(mi_df),14,15)) < 10,T,F)] #只保留肿瘤组织数据
  
  #注释
  annotation_gene <- read.table(file = 'c:/Users/jufei/Documents/TCGA/Rdata/gencode.v22.annotation.gene.probeMap',
                                sep = '\t', header = T)
  rownames(annotation_gene) <- annotation_gene$id
  rownames(expr) <- paste(annotation_gene[rownames(expr),2],rownames(expr),sep = '_')
}

# 判断是否是altered样品
if(F){
  cnames <- substr(colnames(expr),1,15)
  alter_data <- read.table(file = 'c:/Users/jufei/Documents/TCGA/LIHC/LIHC_sample_matrix_FC2.txt',sep = '\t',header = T)[c(1,3)]
  alter_data$studyID.sampleId <- gsub('-','_',alter_data$studyID.sampleId)
  colnames(alter_data)
  
  altered <- alter_data[which(alter_data[,2]==1),]

  group_list <- ifelse(lapply(cnames, function(x){ifelse(sum(grepl(x,altered[,1])),1,0)}),'altered','unaltered')
  
  table(group_list)
  exprSet=na.omit(expr)
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
  rm(mi_df,design,dge,logCPM,v,fit,cont.matrix,fit2,tempOutput,DEG_limma_voom)
}

#作图
if(F){
  #function.R火山图标注基因名
  nrDEG$sign <- ifelse(nrDEG$pvalue < 1e-10 & abs(nrDEG$log2FoldChange) > 0.5,rownames(nrDEG),NA)
  draw_h_v(exprSet,nrDEG,dir_Results,'LIHC_CD45_FC2_altered_vs_unaltered_limma',group_list,1)
  nrDEG_cut <- nrDEG[abs(nrDEG["log2FoldChange"])>=0.5 & nrDEG["pvalue"]<=0.001,-3]
  
  writeData(wb, sheet = "new sheet", nrDEG[,-3], rowNames = T)
  addStyle(wb, 'new sheet', my_style, rows = 1:nrow(nrDEG[,-3]), cols = 1:(ncol(nrDEG[,-3])+1),gridExpand = T)
  saveWorkbook(wb, file = paste0(dir_Results,'LIHC_CD45_FC2_altered_vs_unaltered.xlsx'))
  
  writeData(wb, sheet = "new sheet", nrDEG_cut, rowNames = T)
  addStyle(wb, 'new sheet', my_style, rows = 1:nrow(nrDEG_cut), cols = 1:(ncol(nrDEG[,-3])+1),gridExpand = T)
  saveWorkbook(wb, file = paste0(dir_Results,'LIHC_CD45_FC2_altered_vs_unaltered_cut.xlsx'))
  
  #计算相关性矩阵
  library(Hmisc)
  set <- rownames(nrDEG_cut)
  data <- t(expr[apply(as.data.frame(set), 1, function(x){grep(x,rownames(expr),fixed = T)}),])
  
  #把PDCD1，CTLA4，LAG3和TIGIT的表达数据添加进去
  #data <- cbind(data, t(expr[c("PDCD1_ENSG00000188389.9", "CTLA4_ENSG00000163599.13","TIGIT_ENSG00000181847.10", "LAG3_ENSG00000089692.7"),]))

  #把CD3E，CD4,CD8B和CD8A的表达数据添加进去
  data <- cbind(data, t(expr[c("CD3E_ENSG00000198851.8", "CD4_ENSG00000010610.8","CD8A_ENSG00000153563.14", "CD8B_ENSG00000172116.20"),]))
  
  res2 <- rcorr(as.matrix(data),type = 'spearman')
  res2r <- res2$r
  res2r4i <- res2r[,c("CD3E_ENSG00000198851.8", "CD4_ENSG00000010610.8","CD8A_ENSG00000153563.14", "CD8B_ENSG00000172116.20")]
  res2r4i <- res2r4i[!duplicated(res2r4i),]
  res2r4i_cut4 <- res2r4i[apply(res2r4i,1,function(x) {ifelse(sum(x>=0.4)==4,T,F)}),]
  
  writeData(wb, sheet = "new sheet", res2r4i, rowNames = T)
  addStyle(wb, 'new sheet', my_style, rows = 1:nrow(res2r4i), cols = 1:(ncol(res2r4i)+1),gridExpand = T)
  saveWorkbook(wb, file = paste0(dir_Results,'LIHC_CD45_FC2_Correlation_Analysis_Total.xlsx'))
  
  writeData(wb, sheet = "new sheet", res2r4i_cut4, rowNames = T)
  addStyle(wb, 'new sheet', my_style, rows = 1:nrow(res2r4i_cut4), cols = 1:(ncol(res2r4i_cut4)+1),gridExpand = T)
  saveWorkbook(wb, file = paste0(dir_Results,'LIHC_CD45_FC2_Correlation_Analysis_cut.xlsx'))
  
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
}

#基因注释
if(F){
  ##Gene description
  library(biomaRt)
  FC <- res2r4i_cut4
  FC <- as.data.frame(cbind(do.call(rbind,strsplit(rownames(FC),split = '_')),FC))
  names(FC)[1:2] <- c("external_gene_name","ensembl_gene_id")
  FC$ensembl_gene_id <- do.call(rbind,strsplit(FC$ensembl_gene_id,split = '[.]'))[,1]
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  Description <- getBM(attributes =c("ensembl_gene_id","description"),
                       filters = "ensembl_gene_id", values = FC$ensembl_gene_id,
                       mart = mart)
  Descriptioned <- merge(FC,Description,by="ensembl_gene_id")
  
  writeData(wb, sheet = "new sheet", Descriptioned)
  addStyle(wb, 'new sheet', my_style, rows = 1:nrow(Descriptioned), cols = 1:(ncol(Descriptioned)+1),gridExpand = T)
  saveWorkbook(wb, file = paste0(dir_Results,'LIHC_CD45_FC2_Correlation_Analysis_cut_Descripted.xlsx'))
}