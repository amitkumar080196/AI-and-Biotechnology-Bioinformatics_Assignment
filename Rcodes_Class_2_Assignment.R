getwd()
dir.create("Raw_Data")
dir.create("Script")
dir.create("Results")

input_dir <- "Raw_Data"
output_dir <- "Results"
if (!dir.exists(output_dir)){dir.create(output_dir)}

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")
result_list <- list()



classify_gene <- function(logFC, padj){
  ifelse(padj < 0.05 & logFC < -1, "Dpwnregulation",
         ifelse(padj < 0.05 & logFC > 1, "Upregulated",
                "Not significant"))
}

for(files in files_to_process){
  cat("\nprocessing :", files, "\n")
  file_path <- file.path(input_dir, files)
  data <- read.csv(file_path, header = TRUE)
  cat("Files imported ! , cheking for missing values.....\n")
}

if("padj" %in% names(data)){
  missing_count <- sum(is.na(data$padj))
  cat("Missing values in padj column : " , missing_count , "\n")
  data$padj[is.na(data$padj)] <- mean(data$padj, na.rm = TRUE)
}

if("logFC" %in% names(data)){
  missing_values <- sum(is.na(data$logFC))
  cat("Missing values in logFC : " , missing_count , "\n")
  data$logFC[is.na(data$logFC)] <- mean(data$logFC, na.rm = TRUE)
  
}

data$gene_class <- classify_gene(data$logFC , data$padj)
cat("Gene has been classified succesfully !\n")
result_list[[files]] <- data

ouput_file_path <- file.path(output_dir, paste0("classification"))
write.csv(data, output_file_path, row.names = FALSE)
cat("Result saved to:", ouput_file_path, "\n")

gene_counts <- table(data$gene_class)
cat("Summary counts for", files, ":\n")
print(gene_counts)

result_1 = result_list[1]
result_2 = result_list[2]

save.image(file = "AmitKumar_2_Assignment.RData")

View(result_1)
View(result_2)






