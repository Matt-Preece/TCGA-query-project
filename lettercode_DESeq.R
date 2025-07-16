setwd("G:/Desktop/Corr. Download")

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

cancer <- "BRCA"
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

if (any(list.files() %in% paste(cancer,"_norm_vs_tum.rds", sep = ""))) {

print(paste(cancer,"NORMAL VS TUMOUR DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_norm_vs_tum.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "TCGA DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	gc()

	data <- prep[,!is.na(prep$shortLetterCode)]
	data$shortLetterCode <- factor(data$shortLetterCode, levels = c("NT","TP"))
	data <- DESeqDataSet(data, design = ~ shortLetterCode)
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
	data <- prep[,!is.na(prep$shortLetterCode)]
	data$shortLetterCode <- factor(data$shortLetterCode, levels = c("NT","TP"))
	data <- DESeqDataSet(data, design = ~ shortLetterCode)
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
gene.data <- left_join(gene.count, clinical_data, by = "barcode" )

#create data.frame suitable for plotting
for(i in 1:length(goi)) {
	tmp <- gene.data[c("shortLetterCode",goi[i])]
	tmp$gene <- goi[i]
	colnames(tmp)[2] <- "counts"
	assign(goi[i], value = tmp)
	}

#join data frames for all goi so they can be plotted together
tmp <- mget(goi[1:length(goi)])
feature_vs_counts <- do.call(rbind, tmp)
keep <- !is.na(feature_vs_counts$shortLetterCode) 
feature_vs_counts <- feature_vs_counts[keep,]

#create stat.test data.frame to be used to add p values to final plot, add p values calculated from DESeq
stat.test <- compare_means(counts ~ shortLetterCode, data = feature_vs_counts, group.by = "gene", method = "wilcox")
stat.test <- as.data.frame(stat.test)
stat.test$p.signif <- res$p.signif

#plot graph and save
p1 <- ggplot(feature_vs_counts) +
  geom_boxplot(aes(x = gene, y = counts, fill = shortLetterCode)) +
  stat_pvalue_manual(stat.test, x = "gene", label = "p.signif", y.position = max(feature_vs_counts$counts) + 0.5)

name <- as.character()
for(i in 1:length(goi)) {
  name <- paste(name, goi[i], sep = "_")
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

for(i in 1:length(goi)) {
  
  tmp <- get(goi[i])
  keep <- which(tmp$vital_status %in% "NT")
  tmp_nt <- tmp[keep,]
  keep <- which(tmp$vital_status %in% "TP")
  tmp_tp <- tmp[keep,]
  
  count_csv <- cbind.fill(count_csv,tmp_nt[,2],tmp_tp[,2])
  count_csv <- as.data.frame(count_csv)
  colnames(count_csv)[c(i*2-1,i*2)] <- c(paste(goi[i],"vst_counts","NT",sep = "_"),paste(goi[i],"vst_counts","TP",sep = "_"))
  
}

count_csv[is.na(count_csv)] <- ""

write.csv(count_csv, file = paste(cancer,"_N_v_T_",name,"_vst_counts.csv", sep = ""), row.names = F)
