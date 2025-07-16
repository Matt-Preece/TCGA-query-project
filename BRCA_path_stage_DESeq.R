setwd("G:/Desktop/Corr. Download")

libraries <- c("TCGAbiolinks","SummarizedExperiment","dplyr","ggplot2","ggpubr","rstatix","DESeq2","cowplot")

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

#state which genes are to be investigated under goi
##NOTE: Only use NCBI gene names

cancer <- "BRCA"
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

if (any(list.files() %in% paste(cancer,"_path_stage.rds", sep = ""))) {

print(paste(cancer,"PATHOLOGICAL STAGE DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_path_stage.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[,!is.na(clinical_data$ajcc_pathologic_stage)]
clinical_data$ajcc_pathologic_stage <- gsub("A|B|C","",clinical_data$ajcc_pathologic_stage)
clinical_data$ajcc_pathologic_stage <- gsub(" ","_",clinical_data$ajcc_pathologic_stage)
keep <- which(clinical_data$ajcc_pathologic_stage != "Stage_X")
clinical_data <- clinical_data[,keep]
clinical_data$ajcc_pathologic_stage <- factor(clinical_data$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	clinical_data <- clinical_data[,!is.na(clinical_data$ajcc_pathologic_stage)]
	clinical_data$ajcc_pathologic_stage <- gsub("A|B|C","",clinical_data$ajcc_pathologic_stage)
	clinical_data$ajcc_pathologic_stage <- gsub(" ","_",clinical_data$ajcc_pathologic_stage)
	clinical_data$ajcc_pathologic_stage <- factor(clinical_data$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
	keep <- which(clinical_data$ajcc_pathologic_stage != "Stage_X")
	clinical_data <- clinical_data[,keep]
	gc()

	data <- prep[,!is.na(prep$ajcc_pathologic_stage)]
	data$ajcc_pathologic_stage <- gsub("A|B|C","",data$ajcc_pathologic_stage)
	data$ajcc_pathologic_stage <- gsub(" ","_",data$ajcc_pathologic_stage)
	keep <- which(data$ajcc_pathologic_stage != "Stage_X")
	data <- data[,keep]
	data$ajcc_pathologic_stage <- factor(data$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
	data <- DESeqDataSet(data, design = ~ ajcc_pathologic_stage)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_path_stage.rds", sep = ""))
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
	saveRDS(prep,paste(cancer,".rds", sep = ""))
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	clinical_data <- clinical_data[,!is.na(clinical_data$ajcc_pathologic_stage)]
	clinical_data$ajcc_pathologic_stage <- gsub("A|B|C","",clinical_data$ajcc_pathologic_stage)
	clinical_data$ajcc_pathologic_stage <- gsub(" ","_",clinical_data$ajcc_pathologic_stage)
	keep <- which(clinical_data$ajcc_pathologic_stage != "Stage_X")
	clinical_data <- clinical_data[,keep]
	clinical_data$ajcc_pathologic_stage <- factor(clinical_data$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
	gc()
	
	print("RUNNING DESEQ")
	data <- prep[,!is.na(prep$ajcc_pathologic_stage)]
	data$ajcc_pathologic_stage <- gsub("A|B|C","",data$ajcc_pathologic_stage)
	data$ajcc_pathologic_stage <- gsub(" ","_",data$ajcc_pathologic_stage)
	keep <- which(data$ajcc_pathologic_stage != "Stage_X")
	data <- data[,keep]
	data$ajcc_pathologic_stage <- factor(data$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
	data <- DESeqDataSet(data, design = ~ ajcc_pathologic_stage)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_path_stage.rds", sep = ""))
	gc()

}

}

#extract ENSG code for each of the genes of interest
keep <- which(rowRanges(data)$gene_name %in% goi)
gene_id <- rowRanges(data)$gene_id[keep]

#divide results into individual tables based on pairwise comparisons
res_01 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_0","Stage_I"))
res_02 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_0","Stage_II"))
res_03 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_0","Stage_III"))
res_04 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_0","Stage_IV"))
res_12 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_I","Stage_II"))
res_13 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_I","Stage_III"))
res_14 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_I","Stage_IV"))
res_23 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_II","Stage_III"))
res_24 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_II","Stage_IV"))
res_34 <- results(data, contrast = c("ajcc_pathologic_stage","Stage_III","Stage_IV"))

#individual tables are combined into a list and then removed to preserve memory
res <- list(res_01,res_02,res_03,res_04,res_12,res_13,res_14,res_23,res_24,res_34)
names <- c("0_v_1","0_v_2","0_v_3","0_v_4",
           "1_v_2","1_v_3","1_v_4",
           "2_v_3","2_v_4","3_v_4")
rm(res_01,res_02,res_03,res_04,res_12,res_13,res_14,res_23,res_24,res_34)
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
	assign(paste("ajcc",goi[i],sep = "_"), value = tmp)
 
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




for(i in 1:length(goi)) {
  tmp <- gene.data[c("ajcc_pathologic_stage",goi[i])]
  tmp$gene <- goi[i]
  colnames(tmp)[2] <- "counts"
  
  keep <- !is.na(tmp$ajcc_pathologic_stage)
  pathstage_vs_tpm <- tmp[keep,]
  pathstage_vs_tpm$ajcc_pathologic_stage <- factor(pathstage_vs_tpm$ajcc_pathologic_stage)
  
  #use compare_means to produce data in such a way for stat_pvalue_manual to plot p.signif
  stat.test <- compare_means(counts ~ ajcc_pathologic_stage, data = pathstage_vs_tpm, group.by = "gene", method = "wilcox")
  stat.test <- as.data.frame(stat.test)
  res <- get(paste("ajcc",goi[i],sep = "_"))
  stat.test$p.signif <- res$p.signif
  non.sig <- which(stat.test$p.signif %in% "ns")
  
  #plot chart, problems arise with hide.ns = T and all ns so using if{} else{} as work around
if(length(non.sig) == 10) {
p1 <- ggplot(pathstage_vs_tpm) +
  geom_boxplot(aes(x = ajcc_pathologic_stage, y = counts, fill = ajcc_pathologic_stage))
} else {
  p1 <- ggplot(pathstage_vs_tpm) +
    geom_boxplot(aes(x = ajcc_pathologic_stage, y = counts, fill = ajcc_pathologic_stage)) +
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

ggsave(file = paste(cancer,"_path_stage",name,".png", sep = ""), p2)
