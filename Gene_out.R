a <- read.csv(file = "c:/Users/jufei/Documents/TCGA/GBM/GBM_altered_vs_unaltered_FC1.5_P0.01_Total.csv",header = T)
names(a)[1] <- "external_gene_name"
cnames <- do.call(rbind,strsplit(as.character(a$external_gene_name),split = "_"))
a$external_gene_name <- cnames[,1]

load(file = "CO_SignalIP_TMHMM.Rdata")

descritption_res <- merge(a,co_sp_tm,by="external_gene_name",all.x = T)
write.csv(descritption_res,"GBM_altered_vs_unaltered_SP_TM.csv",row.names = F)

gene_out <- descritption_res[is.na(descritption_res$external_gene_name),1]
View(gene_out)

a[which(a$external_gene_name=="KIAA0748"),1] <- "TESPA1"

gene_out <- FC[unlist(lapply(FC$ensembl_gene_id, function(x){return(ifelse(sum(grepl(x,co_sp_tm$ensembl_gene_id))==0,T,F))})),3]
co_sp_tm_out <- merge(SP_results, TM_results, by = "ensembl_peptide_id",all=T)
Description <- getBM(attributes =c("ensembl_gene_id", "external_gene_name","ensembl_peptide_id","description"),
                     filters = "ensembl_peptide_id", values = co_sp_tm_out$ensembl_peptide_id,
                     mart = mart)
co_sp_tm_out <- merge(co_sp_tm_out, Description, by = "ensembl_peptide_id",all=T)

load(file = "CO_SignalIP_TMHMM.Rdata")
co_sp_tm <- rbind(co_sp_tm,co_sp_tm_out)
save(co_sp_tm,file = "CO_SignalIP_TMHMM.Rdata")
Total_result <- merge(Descriptioned,co_sp_tm[,-7:-8],by="ensembl_gene_id",all.x=T)
write.csv(Total_result,file = "LIHC_TCGA_VS_GTEX_FC0.5_P0.01_Total.csv",row.names = F)
