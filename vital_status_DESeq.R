setwd("G:/Desktop/Corr. Download")
cancer <- "BRCA"
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

libraries <- c("TCGAbiolinks","SummarizedExperiment","dplyr","ggplot2","ggpubr","cowplot","rstatix","DESeq2","BiocParallel","parallel")
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




if (any(list.files() %in% paste(cancer,"_vital_stat.rds", sep = ""))) {

print(paste(cancer,"VITAL STATUS DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_vital_stat.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[!is.na(clinical_data$vital_status),]
clinical_data$vital_status <- factor(clinical_data$vita;_status, levels = c("Alive","Dead"))
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	prep <- prep[,!is.na(prep$vital_status)]
	prep$vital_status <- factor(prep$vital_status, levels = c("Alive","Dead"))

	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)

	data <- DESeqDataSet(prep, design = ~ vital_status)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	gc()
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_vital_stat.rds", sep = ""))
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
	prep <- prep[,!is.na(prep$vital_status)]
	prep$vital_status <- factor(prep$vital_status, levels = c("Alive","Dead"))
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	
	print("RUNNING DESEQ")
	data <- DESeqDataSet(prep, design = ~ vital_status)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	gc()
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_vital_stat.rds", sep = ""))
	gc()

}

}

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
gene_data <- left_join(clinical_data, gene_count, by = "barcode")

for(gene in goi) {
	tmp <- gene_data[c("vital_status",gene)]
	tmp$gene <- gene
	colnames(tmp)[2] <- "counts"
	assign(gene, value = tmp)
	}

tmp <- mget(goi[1:length(goi)])
vital_vs_counts <- do.call(rbind, tmp)
keep <- !is.na(vital_vs_counts$vital_status) 
vital_vs_counts <- vital_vs_counts[keep,]

stat_test <- compare_means(counts ~ vital_status, data = vital_vs_counts, group.by = "gene", method = "wilcox")
stat_test <- as.data.frame(stat_test)
stat_test$p.signif <- res$p.signif

p1 <- ggplot(vital_vs_counts) +
  geom_boxplot(aes(x = gene, y = counts, fill = vital_status)) +
  stat_pvalue_manual(stat_test, x = "gene", label = "p.signif", y.position = max(vital_vs_counts$counts) + 0.5)

name <- as.character()
for(gene in goi) {
  name <- paste(name, gene, sep = "_")
}

ggsave(file = paste(cancer,"_vital_stat_",name,".png", sep = ""), p1)

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
  keep <- which(tmp$vital_status %in% "Alive")
  tmp_a <- tmp[keep,]
  keep <- which(tmp$vital_status %in% "Dead")
  tmp_d <- tmp[keep,]
  
  count_csv <- cbind.fill(count_csv,tmp_a[,2],tmp_d[,2])
  count_csv <- as.data.frame(count_csv)
  colnames(count_csv)[c(i*2-1,i*2)] <- c(paste(gene,"vst_counts","Alive",sep = "_"),paste(gene,"vst_counts","Dead",sep = "_"))
  
}

count_csv[is.na(count_csv)] <- ""
write.csv(count_csv, file = paste(cancer,"_vital_stat_",name,"_vst_counts.csv", sep = ""), row.names = F)
