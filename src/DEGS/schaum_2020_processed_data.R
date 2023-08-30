if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Load the data
setwd("C:/Users/Francisca/Desktop/TeseDeMestrado")
load("MACA_Bulk_deseq2_results.bin")

# View all the tissues
View(names(results_list_MACA_bulk))

# View all comparisons for the Liver tissue
head(results_list_MACA_bulk$Liver)
head(results_list_MACA_bulk$Liver$`15_vs_6`)


# get table of deferentially expressed genes(ressing) and corresponding information
liver_sig_24vs27 <- as.data.frame(results_list_MACA_bulk$Liver$'24_vs_27'$ressig)
brain_sig_24vs27 <- as.data.frame(results_list_MACA_bulk$Brain$'24_vs_27'$ressig)


# Extract tables of all genes and corresponding information
#liver_all_24vs27 <- as.data.frame(results_list_MACA_bulk$Liver$`3_vs_6`$resall)
#brain_all_24vs27 <- as.data.frame(results_list_MACA_bulk$Brain$`3_vs_6`$resall)

# get table of differentially expressed genes(ressing) and corresponding information
#get table with all genes (resall) and corresponding information 
genes_liver_sig_24vs27 <- liver_sig_24vs27$gene_symbol[liver_sig_24vs27$padj<0.05] 
genes_brain_sig_24vs27<- brain_sig_24vs27$gene_symbol[brain_sig_24vs27$padj<0.05] 
#genes_brain_all_24vs27 <- brain_all_28vs3$gene_symbol[liver_sig_28vs3$padj<0.05] 
#genes_liver_all_28vs3 <- liver_all_28vs3$gene_symbol[liver_sig_28vs3$padj<0.05]

genes_liver_sig_24vs27_df <- as.data.frame(genes_liver_sig_24vs27) 
genes_brain_sig_24vs27_df <- as.data.frame(genes_brain_sig_24vs27) 
#genes_brain_all_28vs3_df <- as.data.frame(genes_brain_all_28vs3)
#genes_liver_all_28vs3_df <- as.data.frame(genes_liver_sig_28vs3)

# Save genes name 
write.csv(genes_liver_sig_24vs27_df, "genes_liver_sig_24vs27.csv", row.names = FALSE)
write.csv(genes_brain_sig_24vs27_df, "genes_brain_sig_24vs27.csv", row.names = FALSE)
#write.csv(genes_brain_all_27vs3_df, "genes_brain_all_27vs3.csv", row.names = FALSE)
#write.csv(genes_liver_all_27vs3_df, "genes_liver_all_27vs3.csv", row.names = FALSE)

