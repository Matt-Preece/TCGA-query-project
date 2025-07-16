#setwd("/data/home/qp241207/TCGA_query/raw_data/")
setwd("G:/Desktop/Corr. Download")

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


cancer <- "BRCA"
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

if (any(list.files() %in% paste(cancer,"_vital_stat.rds", sep = ""))) {

print(paste(cancer,"VITAL STATUS DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_vital_stat.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	gc()

	data <- prep[,!is.na(prep$vital_status)]
	data$vital_status <- factor(data$vital_status, levels = c("Dead","Alive"))
	data <- DESeqDataSet(data, design = ~ vital_status)
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
	data <- prep[,!is.na(prep$vital_status)]
	data$vital_status <- factor(data$vital_status, levels = c("Dead","Alive"))
	data <- DESeqDataSet(data, design = ~ vital_status)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	gc()
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "#DATA TYPE#.rds", sep = ""))
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
gene.count <- assay(vsd)
keep <- which(rowRanges(data)$gene_name %in% goi)
gene_id <- rowRanges(data)$gene_id[keep]
keep <- which(rownames(gene.count) %in% gene_id)
gene.count <- gene.count[keep,,drop = T]
gene.count <- as.data.frame(t(gene.count))

#replace column names from ENSG code to NCBI name
for(i in 1:length(goi)) {
  keep <- which(rowRanges(data)$gene_name %in% goi[i])
  tmp <- which(colnames(gene.count) %in% rowRanges(data)$gene_id[keep])
  colnames(gene.count)[tmp] <- goi[i]
}

gene.count$barcode <- rownames(gene.count)
gene.data <- left_join(clinical_data, gene.count, by = "barcode")

for(i in 1:length(goi)) {
	tmp <- gene.data[c("vital_status",goi[i])]
	tmp$gene <- goi[i]
	colnames(tmp)[2] <- "counts"
	assign(goi[i], value = tmp)
	}

tmp <- mget(goi[1:length(goi)])
vital_vs_tpm <- do.call(rbind, tmp)
keep <- !is.na(vital_vs_tpm$vital_status) 
vital_vs_tpm <- vital_vs_tpm[keep,]

stat.test <- compare_means(counts ~ vital_status, data = vital_vs_tpm, group.by = "gene", method = "wilcox")
stat.test <- as.data.frame(stat.test)
stat.test$p.signif <- res$p.signif

p1 <- ggplot(vital_vs_tpm) +
  geom_boxplot(aes(x = gene, y = counts, fill = vital_status)) +
  stat_pvalue_manual(stat.test, x = "gene", label = "p.signif", y.position = max(vital_vs_tpm$counts) + 0.5)

name <- as.character()
for(i in 1:length(goi)) {
  name <- paste(name, goi[i], sep = "_")
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

for(i in 1:length(goi)) {
  
  tmp <- get(goi[i])
  keep <- which(tmp$vital_status %in% "Alive")
  tmp_a <- tmp[keep,]
  keep <- which(tmp$vital_status %in% "Dead")
  tmp_d <- tmp[keep,]
  
  count_csv <- cbind.fill(count_csv,tmp_a[,2],tmp_d[,2])
  count_csv <- as.data.frame(count_csv)
  colnames(count_csv)[c(i*2-1,i*2)] <- c(paste(goi[i],"vst_counts","Alive",sep = "_"),paste(goi[i],"vst_counts","Dead",sep = "_"))
  
}

count_csv[is.na(count_csv)] <- ""

write.csv(count_csv, file = paste(cancer,"_vital_stat_",name,"_vst_counts.csv", sep = ""), row.names = F)
