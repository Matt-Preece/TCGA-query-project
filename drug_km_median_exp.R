setwd("/data/home/qp241207/TCGA_query/raw_data/")

libraries <- c("survival","TCGAbiolinks","SummarizedExperiment","dplyr","ggplot2",
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

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

tmp <- read.delim("cancer.txt", header = F)
cancer <- tmp[,1]
cancer <- toupper(cancer)
tmp <- read.delim("genes.txt", header = F)
goi <-  tmp[,1]
goi <- toupper(goi)
tmp <- read.delim("therapeutics.txt", header = F)
therapeutics <- read.delim("therapeutics.txt", header = T)


prep <- readRDS(paste(cancer, ".rds", sep = ""))
clinical_data <- colData(prep)
clinical_treatment <- clinical_data$treatments
n = length(clinical_treatment)

for(i in 1:dim(therapeutics)[2]) {
tmp <- therapeutics[,i]
keep <- which(tmp != "")
treat <- tmp[keep]
part_tmp <- c()
treat_title <- colnames(therapeutics)[i]

	for(j in 1:n) {
	treat_tmp <- clinical_treatment[[j]]
	participant <- treat_tmp$submitter_id[1]
	participant <- substr(participant,1,12)
		
		if(any(treat_tmp$treatment_type %in% treat | treat_tmp$therapeutic_agents %in% treat)) {
		part_tmp <- c(part_tmp,participant)
		}
	}
	
keep <- which(prep$submitter_id %in% part_tmp)
subset_prep <- prep[,keep]
exp_data <- assay(subset_prep, "unstranded")
subset_clin_data <- colData(subset_prep)
subset_clin_data <- as.data.frame(subset_clin_data)

drug_data <- DESeqDataSetFromMatrix(countData = exp_data,
									colData = subset_clin_data,
									design = ~1)
keep <- rowSums(counts(drug_data)) >= 10
drug_data <- drug_data[keep,]
vsd <- vst(drug_data, blind = F)
drug_data <- assay(vsd)
drug_data <- as.data.frame(drug_data)
gc()


keep <- which(rowRanges(subset_prep)$gene_name %in% goi)
gene_id <- rowRanges(subset_prep)$gene_id[keep]
keep <- which(rownames(drug_data) %in% gene_id)
drug_data <- drug_data[keep,]
drug_data <- t(drug_data)
drug_data <- as.data.frame(drug_data)
drug_data$sample_id <- rownames(drug_data)

#extract necessary data for kaplan-meier analysis, alive/censored and days to death/last follow up and attach expression data
km_plot <- subset_clin_data[,c("vital_status","days_to_death","paper_days_to_last_followup")]
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
km_plot <- left_join(km_plot, drug_data, by = "sample_id")

#rename expression columns from gene ids to gene symbols
for(i in 1:length(goi)) {
  keep <- which(rowRanges(subset_prep)$gene_name %in% goi[i])
  tmp <- which(colnames(km_plot) %in% rowRanges(subset_prep)$gene_id[keep])
  colnames(km_plot)[tmp] <- goi[i]
}

#create the survivor object which can then be used for furher analysis
surv.obj <- Surv(km_plot$time, km_plot$status)

#create two groups based on the expression data, in this case LOW =< median, HIGH > median other methods are possible
for(gene in goi) {
  km_plot[,paste(gene,"strat",sep = "_")] <- ntile(km_plot[gene],2)
  km_plot[,paste(gene,"strat",sep = "_")] <- ifelse(km_plot[,paste(gene,"strat",sep = "_")] == 1,"LOW","HIGH")
  
  km_plot$tmp <- km_plot[,paste(gene,"strat",sep = "_")]
  
  s1 <- survfit(surv.obj ~ tmp, data = km_plot)

suppressWarnings(
  p1 <- ggsurvplot(s1,
                   data = km_plot,
                   pval = T,
                   conf.int = T,
				   palette = c("red","blue"),
				   surv.median.line = "hv",
				   axes.offset = F,
				   legend.labs = c(paste(gene,"High"), paste(gene, "Low")),
				   xlim = c(0,max(km_plot$time)*1.05),
				   title = paste(cancer,gene,treat_title,"Survival Plot")
	)
	)
ggsave(file = paste0("../results/",cancer,"_km_",treat_title,"_",gene,".png"), p1, height = 6, width = 12)
print(paste(cancer,gene,treat_title,"Survival Plot Complete"))
}
}
