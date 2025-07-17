setwd("G:/Desktop/Corr. Download")
cancer <- "BRCA"
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)


libraries <- c("TCGAbiolinks","SummarizedExperiment","dplyr","ggplot2","ggpubr","rstatix","DESeq2")
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

if (any(list.files() %in% paste(cancer,"_norm_vs_tum.rds", sep = ""))) {

print(paste(cancer,"NORMAL VS TUMOUR DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_norm_vs_tum.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[!is.na(clinical_data$shortLetterCode),]
clinical_data$shortLetterCode <- factor(clinical_data$shortLetterCode, levels = c("NT","TP"))
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "TCGA DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	prep <- prep[,!is.na(prep$shortLetterCode)]
	prep$shortLetterCode <- factor(prep$shortLetterCode, levels = c("NT","TP"))
	
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)

	data <- DESeqDataSet(prep, design = ~ shortLetterCode)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	gc()
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_norm_vs_tum.rds", sep = ""))
	gc()

	} else {

	print(paste("STARTING DOWNLOAD",cancer, sep = " - "))

	proj <- paste("TCGA-",cancer,sep = "")
	query <- GDCquery(project = proj,
                      data.category = "Transcriptome Profiling",
                      experimental.strategy = "RNA-Seq",
                      data.type = "Gene Expression Quantification",
                      sample.type = c("Solid Tissue Normal","Primary Tumor"))

	GDCdownload(query, files.per.chunk = 75)
	prep <- GDCprepare(query)
	saveRDS(prep, paste(cancer,".rds", sep = ""))
	prep <- prep[,!is.na(prep$shortLetterCode)]
	prep$shortLetterCode <- factor(prep$shortLetterCode, levels = c("NT","TP"))
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	
	
	print("RUNNING DESEQ")
	data <- DESeqDataSet(prep, design = ~ shortLetterCode)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	gc()
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_norm_vs_tum.rds", sep = ""))
	gc()

}
}


#filter to extract results for goi then create data.frame with results
keep <- which(rowRanges(data)$gene_name %in% goi)
res <- results(data[keep,])
res <- as.data.frame(res)
res$gene_name <- rowRanges(data)[keep,]$gene_name
res <- res[order(res$gene_name),]
res$p.signif <- ifelse(res$padj > 0.05, "ns",
			ifelse(res$padj > 0.01, "*",
			ifelse(res$padj > 0.001, "**",
			ifelse(res$padj > 0.0001, "***","****"
			))))

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

gene_count$barcode <- rownames(gene_count)
gene_data <- left_join(gene_count, clinical_data, by = "barcode" )

#create data.frame suitable for plotting
for(gene in goi) {
	tmp <- gene_data[c("shortLetterCode",gene)]
	tmp$gene <- gene
	colnames(tmp)[2] <- "counts"
	assign(gene, value = tmp)
	}

#join data frames for all goi so they can be plotted together
tmp <- mget(goi[1:length(goi)])
NT_TP_vs_counts <- do.call(rbind, tmp)
keep <- !is.na(NT_TP_vs_counts$shortLetterCode) 
NT_TP_vs_counts <- NT_TP_vs_counts[keep,]

#create stat_test data.frame to be used to add p values to final plot, add p values calculated from DESeq
stat_test <- compare_means(counts ~ shortLetterCode, data = NT_TP_vs_counts, group.by = "gene", method = "wilcox")
stat_test <- as.data.frame(stat_test)
stat_test$p.signif <- res$p.signif

#plot graph and save
p1 <- ggplot(NT_TP_vs_counts) +
  geom_boxplot(aes(x = gene, y = counts, fill = shortLetterCode)) +
  stat_pvalue_manual(stat_test, x = "gene", label = "p.signif", y.position = max(NT_TP_vs_counts$counts) + 0.5)

name <- as.character()
for(gene in goi) {
  name <- paste(name, gene, sep = "_")
}

ggsave(file = paste(cancer,"_N_v_T_",name,".png", sep = ""), p1)

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

count_csv <- data.frame()
for(gene in goi) {
  
  tmp <- get(gene)
  keep <- which(tmp$vital_status %in% "NT")
  tmp_nt <- tmp[keep,]
  keep <- which(tmp$vital_status %in% "TP")
  tmp_tp <- tmp[keep,]
  count_csv <- cbind.fill(count_csv,tmp_nt[,2],tmp_tp[,2])
  count_csv <- as.data.frame(count_csv)
  colnames(count_csv)[c(i*2-1,i*2)] <- c(paste(gene,"vst_counts","NT",sep = "_"),paste(gene,"vst_counts","TP",sep = "_"))
  
}

count_csv[is.na(count_csv)] <- ""
write.csv(count_csv, file = paste(cancer,"_N_v_T_",name,"_vst_counts.csv", sep = ""), row.names = F)
