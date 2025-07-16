#set directory you will be working from
##NOTE: when writing directory paths in R use "/" not "\"
setwd("G:/Desktop/Corr. Download")

#state which genes are to be investigated under goi
##NOTE: Only use NCBI gene names
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

#this code is written specifically for breast cancer data do not change
cancer <- "BRCA"

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

if (any(list.files() %in% paste(cancer,"_path_stage.rds", sep = ""))) {

print(paste(cancer,"PATHOLOGICAL STAGE DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_path_stage.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[!is.na(clinical_data$ajcc_pathologic_stage),]
clinical_data$ajcc_pathologic_stage <- gsub("A|B|C","",clinical_data$ajcc_pathologic_stage)
clinical_data$ajcc_pathologic_stage <- gsub(" ","_",clinical_data$ajcc_pathologic_stage)
keep <- which(clinical_data$ajcc_pathologic_stage != "Stage_X")
clinical_data <- clinical_data[keep,]
clinical_data$ajcc_pathologic_stage <- factor(clinical_data$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
rm(prep)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "TCGA DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer, ".rds", sep = ""))
	prep <- prep[!is.na(prep$ajcc_pathologic_stage),]
	prep$ajcc_pathologic_stage <- gsub("A|B|C","",prep$ajcc_pathologic_stage)
	prep$ajcc_pathologic_stage <- gsub(" ","_",prep$ajcc_pathologic_stage)
	prep$ajcc_pathologic_stage <- factor(prep$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
	keep <- which(prep$ajcc_pathologic_stage != "Stage_X")
	prep <- prep[,keep]
		
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)

	data <- DESeqDataSet(prep, design = ~ ajcc_pathologic_stage)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_path_stage.rds", sep = ""))
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
	saveRDS(prep,paste(cancer,".rds", sep = ""))
	prep <- prep[!is.na(prep$ajcc_pathologic_stage),]
	prep$ajcc_pathologic_stage <- gsub("A|B|C","",prep$ajcc_pathologic_stage)
	prep$ajcc_pathologic_stage <- gsub(" ","_",prep$ajcc_pathologic_stage)
	prep$ajcc_pathologic_stage <- factor(prep$ajcc_pathologic_stage, levels = c("Stage_0","Stage_I","Stage_II","Stage_III","Stage_IV"))
	keep <- which(prep$ajcc_pathologic_stage != "Stage_X")
	prep <- prep[,keep]
		
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	gc()
	
	print("RUNNING DESEQ")
	data <- DESeqDataSet(prep, design = ~ ajcc_pathologic_stage)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	rm(prep)
	data <- DESeq(data)
	saveRDS(data,paste(cancer, "_path_stage.rds", sep = ""))
	gc()

}
}

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
	assign(paste("ajcc",gene,sep = "_"), value = tmp)
 
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

#combine clinical_data and gene_count to create data frame for plotting
gene_count$barcode <- rownames(gene_count)
gene_data <- left_join(clinical_data, gene_count, by = "barcode")

for(gene in goi) {
  pathstage_vs_count <- gene_data[c("ajcc_pathologic_stage",gene)]
  pathstage_vs_count$gene <- gene
  colnames(pathstage_vs_count)[2] <- "counts"
  
  #use compare_means to produce data in such a way for stat_pvalue_manual to plot p.signif
  stat_test <- compare_means(counts ~ ajcc_pathologic_stage, data = pathstage_vs_count, group.by = "gene", method = "wilcox")
  stat_test <- as.data.frame(stat_test)
  res <- get(paste("ajcc",gene,sep = "_"))
  stat_test$p.signif <- res$p.signif
  non_sig <- which(stat_test$p.signif %in% "ns")
  
#plot chart, problems arise with hide.ns = T and all ns so using if{} else{} as work around, also stat_pvalue_manual() uses values from all p, p.adj and p.format columns, as the p.signif is determined using
# DESeq() rather than compare_means() there can be discrepancies hence the short sig <- lines to ensure there is no conflict

	if(length(non_sig) == 10) {
	
	p1 <- ggplot(pathstage_vs_count) +
  	geom_boxplot(aes(x = ajcc_pathologic_stage, y = counts, fill = ajcc_pathologic_stage))
	
	} else {

	sig <- which(stat_test$p.signif != "ns")
	stat_test[sig,5:7] <- 0.00001
	p1 <- ggplot(pathstage_vs_count) +
    		geom_boxplot(aes(x = ajcc_pathologic_stage, y = counts, fill = ajcc_pathologic_stage)) +
    		stat_pvalue_manual(stat_test, label = "p.signif", hide.ns = T,
                       y.position = max(pathstage_vs_count$counts) + 0.5, step.increase = 0.15) 
	}
	
  assign(gene, value = p1)

}

#plot graph and save
tmp <- mget(goi[1:length(goi)])
p2 <- plot_grid(plotlist = tmp, ncol = 1, nrow = length(goi))

name <- as.character()
for(gene in goi) {
  name <- paste(name, gene, sep = "_")
}

ggsave(file = paste(cancer,"_path_stage",name,".png", sep = ""), p2)
