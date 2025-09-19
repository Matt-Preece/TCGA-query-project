#set directory you will be working from
##NOTE: when writing directory paths in R use "/" not "\"
setwd("G:/Desktop/Corr. Download")

#state which genes are to be investigated under goi
##NOTE: Only use NCBI gene names
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

#this code is written specifically for breast cancer data do not change
cancer <- "BRCA"

libraries <- c("TCGAbiolinks","SummarizedExperiment","dplyr","ggplot2","ggpubr","cowplot","rstatix","DESeq2","BiocParallel")
for (lib in libraries) {
  if (require(package = lib, character.only = TRUE)) {
    successful <- "Successful"
  } else {
    installing <- "Installing"
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkgs = lib, suppressUpdates = T)
    library(lib, character.only = TRUE )
  }
}

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
						  
if (any(list.files() %in% paste(cancer,"_subtype.rds", sep = ""))) {

print(paste(cancer,"SUBTYPE DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_path_stage.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[!is.na(clinical_data$paper_BRCA_Subtype_PAM50),]
clinical_data$paper_BRCA_Subtype_PAM50 <- factor(data$paper_BRCA_Subtype_PAM50, levels = c("Normal","Basal","Her2","LumA","LumB"))
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	prep <- prep[,!is.na(prep$paper_BRCA_Subtype_PAM50)]
	
	gc()

	data <- prep[,!is.na(prep$paper_BRCA_Subtype_PAM50)]
	data$paper_BRCA_Subtype_PAM50 <- factor(data$paper_BRCA_Subtype_PAM50, levels = c("Normal","Basal","Her2","LumA","LumB"))
	data <- DESeqDataSet(data, design = ~ paper_BRCA_Subtype_PAM50)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	data <- DESeq(data)
	gc()
	saveRDS(data,paste(cancer, "_subtype.rds", sep = ""))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	rm(prep)
	gc()

	} else {

	print(paste("STARTING DOWNLOAD",cancer, sep = " - "))

	proj <- paste("TCGA-",cancer,sep = "")
	data <- GDCquery(project = proj,
                      data.category = "Transcriptome Profiling",
                      experimental.strategy = "RNA-Seq",
                      data.type = "Gene Expression Quantification",
                      sample.type = c("Solid Tissue Normal","Primary Tumor"))

	GDCdownload(data, files.per.chunk = 75)
	prep <- GDCprepare(data)
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	saveRDS(paste(cancer,".rds", sep = ""))
	gc()
	
	print("RUNNING DESEQ")
	data <- prep[,!is.na(prep$paper_BRCA_Subtype_PAM50)]
	data$paper_BRCA_Subtype_PAM50 <- factor(data$paper_BRCA_Subtype_PAM50, levels = c("Normal","Basal","Her2","LumA","LumB"))
	data <- DESeqDataSet(data, design = ~ paper_BRCA_Subtype_PAM50)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	gc()
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_subtype.rds", sep = ""))
	gc()

}

}

keep <- which(rowRanges(data)$gene_name %in% goi)
gene_id <- rowRanges(data)$gene_id[keep]

res_01 <- results(data, contrast = c("paper_BRCA_Subtype_PAM50","Normal","Basal"))
res_02 <- results(data, contrast = c("paper_BRCA_Subtype_PAM50","Normal","Her2"))
res_03 <- results(data, contrast = c("paper_BRCA_Subtype_PAM50","Normal","LumA"))
res_04 <- results(data, contrast = c("paper_BRCA_Subtype_PAM50","Normal","LumB"))


#individual tables are combined into a list and then removed to preserve memory
res <- list(res_01,res_02,res_03,res_04)
names <- c("Norm_v_Bas","Norm_v_Her2","Norm_v_LumA","Norm_v_LumB")
rm(res_01,res_02,res_03,res_04)
gc()

#narrow down results to only those from genes of interest
#results for each individual gene are combined into a single table of all the pairwise comparisons and p.signif symbols are added according to padj
#anaylsis done on full transcriptome to maintian statistical power and reliability of results
for(gene in goi) {
  keep <- which(rowRanges(data)$gene_name %in% gene)
  tmp <- sapply(res,"[",keep,,simplify = F)
  tmp <- do.call(rbind,tmp)
  rownames(tmp) <- names
  tmp$gene_name <- rowRanges(data)[keep,]$gene_name
  tmp$p.signif <- ifelse(tmp$padj > 0.05, "ns",
                         ifelse(tmp$padj > 0.01, "*",
                                ifelse(tmp$padj > 0.001, "**",
                                       ifelse(tmp$padj > 0.0001, "***","****"
                                       ))))
  assign(paste("subtype",gene,sep = "_"), value = tmp)
  
}

#extract normalised gene counts from data set for plotting
vsd <- vst(data, blind = F)
gene_count <- assay(vsd)
keep <- which(rowRanges(data)$gene_name %in% goi)
gene_id <- rowRanges(data)$gene_id[keep]
keep <- which(rownames(gene_count) %in% gene_id)
gene_count <- gene_count[keep,,drop = T]
gene_count <- as.data.frame(t(gene_count))

#replace column names from ENSG code to NCBI name
for(gene in goi) {
  keep <- which(rowRanges(data)$gene_name %in% gene)
  tmp <- which(colnames(gene_count) %in% rowRanges(data)$gene_id[keep])
  colnames(gene_count)[tmp] <- gene
}

#no log2 required as the counts have already been normalised
#combine BRCA.subj and gene_count to create data frame for plotting
gene_count$barcode <- rownames(gene_count)
gene_data <- left_join(clinical_data, gene_count, by = "barcode")
gene_data <- gene_data[!is.na(gene_data$paper_BRCA_Subtype_PAM50),]
gene_data$paper_BRCA_Subtype_PAM50 <- factor(gene_data$paper_BRCA_Subtype_PAM50, levels = c("Normal","Basal","Her2","LumA","LumB"))
                                             
for(gene in goi) {
  tmp <- gene_data[c("paper_BRCA_Subtype_PAM50",gene)]
  tmp$gene <- gene
  colnames(tmp)[2] <- "counts"
  
  keep <- !is.na(tmp$paper_BRCA_Subtype_PAM50)
  subtype_vs_counts <- tmp[keep,]
  subtype_vs_counts$paper_BRCA_Subtype_PAM50 <- factor(subtype_vs_counts$paper_BRCA_Subtype_PAM50)
  
  #use compare_means to produce data in such a way for stat_pvalue_manual to plot p.signif
  stat_test <- compare_means(counts ~ paper_BRCA_Subtype_PAM50, data = subtype_vs_counts, group.by = "gene", method = "wilcox")
  stat_test <- as.data.frame(stat_test[1:4,])
  res <- get(paste("subtype",gene,sep = "_"))
  stat_test$p.signif <- res$p.signif
  non.sig <- which(stat_test$p.signif %in% "ns")
  
  #plot chart, problems arise with hide.ns = T and all ns so using if{} else{} as work around
  if(length(non.sig) == 4) {
    p1 <- ggplot(subtype_vs_counts) +
      geom_boxplot(aes(x = paper_BRCA_Subtype_PAM50, y = counts, fill = paper_BRCA_Subtype_PAM50))
  } else {
    sig <- which(stat_test$p.signif != "ns")
    stat_test[sig,5:7] <- 0.0001
      p1 <- ggplot(subtype_vs_counts) +
      geom_boxplot(aes(x = paper_BRCA_Subtype_PAM50, y = counts, fill = paper_BRCA_Subtype_PAM50)) +
      stat_pvalue_manual(stat_test, label = "p.signif", hide.ns = T,
                         y.position = max(subtype_vs_counts$counts) + 0.5, step.increase = 0.15) 
  }
  ggsave(file = paste(cancer,"_subtype",gene,".png", sep = ""), p1)
  subtype_df <- data.frame()
  keep <- which(subtype_vs_count$paper_BRCA_Subtype_PAM50 %in% "Normal")
  subtype_df <- cbind.fill(subtype_df,subtype_vs_count[keep,"counts"])
  keep <- which(subtype_vs_count$paper_BRCA_Subtype_PAM50 %in% "Basal")
  subtype_df <- cbind.fill(subtype_df,subtype_vs_count[keep,"counts"])
  keep <- which(subtype_vs_count$paper_BRCA_Subtype_PAM50 %in% "Her2")
  subtype_df <- cbind.fill(subtype_df,subtype_vs_count[keep,"counts"])
  keep <- which(subtype_vs_count$paper_BRCA_Subtype_PAM50 %in% "LumA")
  subtype_df <- cbind.fill(subtype_df,subtype_vs_count[keep,"counts"])
  keep <- which(subtype_vs_count$paper_BRCA_Subtype_PAM50 %in% "LumB")
  subtype_df <- cbind.fill(subtype_df,subtype_vs_count[keep,"counts"])
  colnames(subtype_df) <- c("Normal","Basal","Her2","LumA","LumB")
  subtype_df[is.na(subtype_df)] <- ""
  write.csv(subtype_df, file = paste0("../results/",cancer,"_path_stage_",gene,".csv"))
	
}



