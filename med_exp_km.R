setwd("G:/Desktop/Corr. Download")

libraries <-c("survival","TCGAbiolinks","SummarizedExperiment","dplyr","ggplot2",
              "ggpubr","ggfortify","rstatix","DESeq2","survminer","cowplot")


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

#set genes to be investigated
cancer <- "BRCA"
goi <- c("ATAT1","HDAC6","SIRT2")
goi <- toupper(goi)

if (any(list.files() %in% paste(cancer,"_surv_anal.rds", sep = ""))) {

print(paste(cancer,"SURVIVAL ANALYSIS DESEQ COMPLETE", sep = " "))
data <- readRDS(paste(cancer, "_surv_anal.rds", sep = ""))
prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_data <- as.data.frame(clinical_data)
gc()

} else {

	if (any(list.files() %in% paste(cancer, ".rds", sep = ""))) {

	print(paste(cancer, "DATA AVAILABLE - RUNNING DESEQ", sep = " "))
	prep <- readRDS(paste(cancer,".rds", sep = ""))
	exp_data <- assay(prep, "unstranded")
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	gc()
	
#create DESeqDataSet and then normalise by variance stabilising transformation which gives the best form of the data for non-DE based analysis
	data <- DESeqDataSetFromMatrix(countData = exp_data,
                               colData = clinical_data,
                               design = ~1)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	vsd <- vst(data, blind = F)
	data <- assay(vsd)
	data <- as.data.frame(data)
	saveRDS(data, paste(cancer,"_surv_anal.rds", sep = ""))
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
	exp_data <- assay(prep, "unstranded")
	clinical_data <- colData(prep)
	clinical_data <- as.data.frame(clinical_data)
	saveRDS(prep, paste(cancer,".rds", sep = ""))
	gc()
	
	print("RUNNING DESEQ")
	data <- DESeqDataSetFromMatrix(countData = exp_data,
                               colData = clinical_data,
                               design = ~1)
	keep <- rowSums(counts(data)) >= 10
	data <- data[keep,]
	vsd <- vst(data, blind = F)
	data <- assay(vsd)
	data <- as.data.frame(data)
	saveRDS(data, paste(cancer,"_surv_anal.rds", sep = ""))
	gc()
}

}


#filter the dataset to only include goi
keep <- which(rowRanges(prep)$gene_name %in% goi)
gene_id <- rowRanges(prep)$gene_id[keep]
keep <- which(rownames(data) %in% gene_id)
data <- data[keep,]
data <- t(data)
data <- as.data.frame(data)
data$sample_id <- rownames(data)

#extract necessary data for kaplan-meier analysis, alive/censored and days to death/last follow up and attach expression data
km_plot <- clinical_data[,c("vital_status","days_to_death","paper_days_to_last_followup")]
km_plot <- as.data.frame(km_plot)
km_plot[,c(2:3)] <- sapply(km_plot[,c(2:3)], as.numeric)
km_plot$status <- ifelse(km_plot$vital_status == "Alive",1,2)
km_plot$time <- ifelse(km_plot$vital_status == "Alive",
                       km_plot$paper_days_to_last_followup,
                       km_plot$days_to_death)
keep <- !is.na(km_plot$time)
km_plot <- km_plot[keep,]
keep <- !is.na(km_plot$status)
km_plot <- km_plot[keep,]
km_plot$sample_id <- rownames(km_plot)
km_plot <- left_join(km_plot, data, by = "sample_id")

#rename expression columns from gene ids to gene symbols
for(i in 1:length(goi)) {
  keep <- which(rowRanges(prep)$gene_name %in% goi[i])
  tmp <- which(colnames(km_plot) %in% rowRanges(prep)$gene_id[keep])
  colnames(km_plot)[tmp] <- goi[i]
}

#create the survivor object which can then be used for furher analysis
surv.obj <- Surv(km_plot$time, km_plot$status)

#create two groups based on the expression data, in this case LOW =< median, HIGH > median other methods are possible
for(i in 1:length(goi)) {
  km_plot[,paste(goi[i],"strat",sep = "_")] <- ntile(km_plot[goi[i]],2)
  km_plot[,paste(goi[i],"strat",sep = "_")] <- ifelse(km_plot[,paste(goi[i],"strat",sep = "_")] == 1,"LOW","HIGH")
  
  km_plot$tmp <- km_plot[,paste(goi[i],"strat",sep = "_")]
  
  s1 <- survfit(surv.obj ~ tmp, data = km_plot)
  p1 <- ggsurvplot(s1,
                   data = km_plot,
                   pval = T,
                   risk.table = T,
                   conf.int = T)
  assign(goi[i], value = p1)
}


#stitch all plots together and save $path is necessary as ggsurvplots can't be saved. $plot extracts just the plot and not the additional data
tmp <- mget(goi)
for(i in 1:length(goi)){
  tmp[[i]] <- tmp[[i]]$plot  
}

p2 <- plot_grid(plotlist = tmp, ncol = 1, nrow = length(goi))

name <- as.character()
for(i in 1:length(goi)) {
  name <- paste(name, goi[i], sep = "_")
}

ggsave(file = paste(cancer,"_kaplan_meier",name,".png", sep = ""), p2)
