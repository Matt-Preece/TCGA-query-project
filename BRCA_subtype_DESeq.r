setwd("G:/Desktop/Corr. Download")

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

cancer <- "BRCA"
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

if (any(list.files() %in% paste(cancer,"_subtype.rds", sep = ""))) {

print(paste(cancer,"SUBTYPE DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_subtype.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
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
for(i in 1:length(goi)) {
  keep <- which(rowRanges(data)$gene_name %in% goi[i])
  tmp <- sapply(res,"[",keep,,simplify = F)
  tmp <- do.call(rbind,tmp)
  rownames(tmp) <- names
  tmp$gene_name <- rowRanges(data)[keep,]$gene_name
  tmp$p.signif <- ifelse(tmp$padj > 0.05, "ns",
                         ifelse(tmp$padj > 0.01, "*",
                                ifelse(tmp$padj > 0.001, "**",
                                       ifelse(tmp$padj > 0.0001, "***","****"
                                       ))))
  assign(paste("subtype",goi[i],sep = "_"), value = tmp)
  
}

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

#no log2 required as the counts have already been normalised
#combine BRCA.subj and gene.count to create data frame for plotting
gene.count$barcode <- rownames(gene.count)
gene.data <- left_join(clinical_data, gene.count, by = "barcode")
gene.data <- gene.data[!is.na(gene.data$paper_BRCA_Subtype_PAM50),]
gene.data$paper_BRCA_Subtype_PAM50 <- factor(gene.data$paper_BRCA_Subtype_PAM50, levels = c("Normal","Basal","Her2","LumA","LumB"))
                                             
for(i in 1:length(goi)) {
  tmp <- gene.data[c("paper_BRCA_Subtype_PAM50",goi[i])]
  tmp$gene <- goi[i]
  colnames(tmp)[2] <- "counts"
  
  keep <- !is.na(tmp$paper_BRCA_Subtype_PAM50)
  pathstage_vs_tpm <- tmp[keep,]
  pathstage_vs_tpm$paper_BRCA_Subtype_PAM50 <- factor(pathstage_vs_tpm$paper_BRCA_Subtype_PAM50)
  
  #use compare_means to produce data in such a way for stat_pvalue_manual to plot p.signif
  stat.test <- compare_means(counts ~ paper_BRCA_Subtype_PAM50, data = pathstage_vs_tpm, group.by = "gene", method = "wilcox")
  stat.test <- as.data.frame(stat.test[1:4,])
  res <- get(paste("subtype",goi[i],sep = "_"))
  stat.test$p.signif <- res$p.signif
  non.sig <- which(stat.test$p.signif %in% "ns")
  
  #plot chart, problems arise with hide.ns = T and all ns so using if{} else{} as work around
  if(length(non.sig) == 4) {
    p1 <- ggplot(pathstage_vs_tpm) +
      geom_boxplot(aes(x = paper_BRCA_Subtype_PAM50, y = counts, fill = paper_BRCA_Subtype_PAM50))
  } else {
    p1 <- ggplot(pathstage_vs_tpm) +
      geom_boxplot(aes(x = paper_BRCA_Subtype_PAM50, y = counts, fill = paper_BRCA_Subtype_PAM50)) +
      stat_pvalue_manual(stat.test, label = "p.signif", hide.ns = T,
                         y.position = max(pathstage_vs_tpm$counts) + 0.5, step.increase = 0.15) 
  }
  assign(goi[i], value = p1)
  
}

#plot graph and save
tmp <- mget(goi[1:length(goi)])
p2 <- plot_grid(plotlist = tmp, ncol = 1, nrow = length(goi))

name <- as.character()
for(i in 1:length(goi)) {
  name <- paste(name, goi[i], sep = "_")
}

ggsave(file = paste(cancer,"_subtype",name,".png", sep = ""), p2)

