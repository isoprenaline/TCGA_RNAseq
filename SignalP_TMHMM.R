#信号肽序列与跨膜区段分析

#基因ID转氨基酸序列
if(F){
  library(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  EnsemblGeneID <- read.csv('C:/Users/jufei/Documents/EnsemblGeneID.csv', header=F, stringsAsFactors = FALSE)
  EnsemblGeneID <- gene_out
  
  GeneID2PeptideID <- getBM(attributes = "ensembl_peptide_id",
                            filters = "ensembl_gene_id", values = EnsemblGeneID,  #"external_gene_name","ensembl_gene_id"
                            mart = mart)
  protein <- getSequence(id= GeneID2PeptideID,
                         type="ensembl_peptide_id",
                         seqType="peptide", 
                         mart= mart)
  protein_out <- protein
  
  save(protein,file = "Ensembl_Peptides_sequence_data.Rdata")
  exportFASTA(protein, file = "C:/Users/jufei/Documents/Datas/AA_sequecnes_for_SignalP_out.fa")
  exportFASTA(protein, file = "C:/Users/jufei/Documents/Datas/AA_sequecnes_for_TMHMM.fa")
}

#文件分割
if(F){
  
  data_split <- as.data.frame(protein)
  nprotein <- 5000
  file_name <- "C:/Users/jufei/Documents/Datas/AA_sequecnes_for_SP_"
  
  
  file_spilt <- function(x){
    a <- (x-1)*nprotein+1
    b <- x*nprotein
    c <- data_split[a:b,]
    d <- paste0(file_name,b,".fa")
    exportFASTA(c, file = d)
  }

  file_spilt_end <- function(x){
    a <- (x-1)*nprotein+1
    b <- nrow(data_split)
    c <- data_split[a:b,]
    d <- paste0(file_name,b,".fa")
    exportFASTA(c, file = d)
  }
  
  for (i in 1:ceiling(nrow(data_split)/nprotein)) {
    if (i==ceiling(nrow(data_split)/nprotein)) {
      file_spilt_end(i)
    }
    else{
      file_spilt(i)
    }
    print(i)
  }
}

#信号肽数据整理
if(F){
  dir_sp <- "c:/Users/jufei/Documents"
  SP_results = lapply(list.files(path = dir_sp ,pattern='Results_of_SP*')
                      , function(x){
                        result <- read.table(file = file.path(dir_sp,x),sep = '\t',header = F)
                        return( result )
                      })
  SP_results <- do.call(rbind,SP_results)
  names(SP_results) <- c("ensembl_peptide_id", "Prediction","SP(Sec/SPI)",	"OTHER",	"CS_Position")
  save(SP_results,file = "SP_results.Rdata")
  
  load(file = "SP_results.Rdata")
  
  geneID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","ensembl_peptide_id", "description"),
                  filters = "ensembl_peptide_id", values = SP_results$ensembl_peptide_id,
                  mart = mart)
  
  co_sp <- merge(SP_results, geneID, by = "ensembl_peptide_id")
}


#去除信号肽序列进行TMHMM分析
if(F){
  protein_sp <- co_sp[grep("SP(",co_sp$Prediction,fixed = T),c(1,5)]
  cs_position <- lapply(as.character(protein_sp[,2]), function(x){return(strsplit(x,split = '-',fixed = T)[[1]][2])})
  cs_position <- do.call(rbind,lapply(cs_position, function(x){return(strsplit(x,split = '.',fixed = T)[[1]][1])}))
  protein_sp[,2] <- cs_position
  
  protein_cut_sp <- lapply(protein_sp$ensembl_peptide_id, function(x){substr(protein[which(x==protein$ensembl_peptide_id),1],protein_sp[which(x==protein_sp$ensembl_peptide_id),2],nchar(protein[which(x==protein$ensembl_peptide_id),1]))})
  protein_cut_sp_TM <- as.data.frame(cbind(do.call(rbind,protein_cut_sp),as.character(protein_sp$ensembl_peptide_id)))
  names(protein_cut_sp_TM) <- c('peptide','ensembl_peptide_id')
  protein_no_sp_TM <- protein[-(do.call(rbind,lapply(protein_cut_sp_TM[,2],function(x){return(which(x==protein$ensembl_peptide_id))}))),]

  exportFASTA(protein_cut_sp_TM, file = "C:/Users/jufei/Documents/Datas/AA_sequecnes_for_TM_cutsp.fa")
  exportFASTA(protein_no_sp_TM, file = "C:/Users/jufei/Documents/Datas/AA_sequecnes_for_TM_nosp.fa")
  
  #读取TMHMM分析数据
  dir_TM <- "c:/Users/jufei/Documents"
  TM_results = lapply(list.files(path = dir_TM ,pattern='Result_of_TM_*')
                      , function(x){
                        result <- read.table(file = file.path(dir_TM,x),sep = '\t',header = F,fill = T)
                        return( result )
                      })
  TM_results <- do.call(rbind,TM_results)
  names(TM_results) <- c("ensembl_peptide_id","Len","ExpAA","First60","PredHel","Topology")
  save(TM_results,file = "TM_results.Rdata")
  
  load(file = "TM_results.Rdata")
  
  Full_results_TM <- merge(Full_results,TM_results[,-c(2,3,4)],by="ensembl_peptide_id",all = T)
  
  write.csv(Full_results_TM,"PAAD_SP_TMHMM_limma_DEG_Tumor_TCGA_vs_Normal_GTEx_Cut.csv")
  
}

co_sp_tm <- merge(co_sp, TM_results, by = "ensembl_peptide_id",all=T)
save(co_sp_tm,file = "CO_SignalIP_TMHMM.Rdata")
load(file = "CO_SignalIP_TMHMM.Rdata")





