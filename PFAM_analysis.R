data <- read.csv("E:/R work/PFAM data/At_stressand backmain_domain.csv", header = TRUE)
# Display the first few rows of the dataset
head(data)
# Count unique gene IDs
unique_gene_count <- length(unique(data$gene_id))
# Count unique PFAM IDs for stress and background
unique_pfam_stress_count <- length(unique(data$pfam_id_stress))
unique_pfam_back_count <- length(unique(data$pfam_id_back))

# Print summary information
cat("Total Unique Gene IDs:", unique_gene_count, "\n")
cat("Unique PFAM IDs for Stress:", unique_pfam_stress_count, "\n")
cat("Unique PFAM IDs for Background:", unique_pfam_back_count, "\n")
# Calculate the total number of unique PFAM IDs in both conditions
total_unique_pfam_ids <- unique(c(data$pfam_id_stress, data$pfam_id_back))

# Calculate the number of unique PFAM IDs in the "Stress" condition
unique_pfam_ids_stress <- unique(data$pfam_id_stress)

# Create a function to perform enrichment analysis
enrichment_analysis <- function(pfam_id) {
  # Calculate the number of genes with the PFAM ID in the "Stress" condition
  genes_with_pfam_stress <- sum(data$pfam_id_stress == pfam_id)
  
  # Calculate the number of genes without the PFAM ID in the "Stress" condition
  genes_without_pfam_stress <- sum(data$pfam_id_stress != pfam_id)
  
  # Calculate the number of genes with the PFAM ID in the "Background" condition
  genes_with_pfam_back <- sum(data$pfam_id_back == pfam_id)
  
  # Calculate the number of genes without the PFAM ID in the "Background" condition
  genes_without_pfam_back <- sum(data$pfam_id_back != pfam_id)
  
  # Perform Fisher's exact test for enrichment analysis
  enrichment_result <- fisher.test(matrix(c(genes_with_pfam_stress, genes_without_pfam_stress,
                                            genes_with_pfam_back, genes_without_pfam_back),
                                          ncol = 2))
  
  # Return the PFAM ID and the p-value
  return(c(pfam_id, enrichment_result$p.value))
}

# Perform enrichment analysis for each unique PFAM ID
enrichment_results <- sapply(unique_pfam_ids_stress, enrichment_analysis)

# Create a data frame to store the results
enrichment_df <- data.frame(PFAM_ID = enrichment_results[1, ],
                            P_Value = enrichment_results[2, ])

# Filter significant results (e.g., p-value < 0.05)
significant_enrichment <- enrichment_df[enrichment_df$P_Value < 0.05, ]

# Sort the results by p-value
significant_enrichment <- significant_enrichment[order(significant_enrichment$P_Value), ]

# Print the significant enrichment results
print(significant_enrichment)
# Assuming you have your data in a data frame named 'data'

# Perform association analysis for the "Stress" condition
association_stress <- table(data$gene_id, data$pfam_id_stress)

# Perform association analysis for the "Background" condition
association_background <- table(data$gene_id, data$pfam_id_back)

# View the association analysis results for the "Stress" condition
View(association_stress)

# View the association analysis results for the "Background" condition
View(association_background)
# Assuming you have your data in a data frame named 'data'

# Calculate the frequency of each PFAM ID in the "Stress" condition for each gene
pfam_freq_stress <- table(data$gene_id, data$pfam_id_stress)

# Calculate the frequency of each PFAM ID in the "Background" condition for each gene
pfam_freq_background <- table(data$gene_id, data$pfam_id_back)

# Calculate the correlation between gene frequencies in the "Stress" and "Background" conditions
correlation <- cor(pfam_freq_stress, pfam_freq_background)

# View the correlation matrix
View(correlation)
# Create a heatmap for the correlation matrix
heatmap(correlation,
        col = colorRampPalette(c("blue", "white", "red"))(100),  # Define color palette
        main = "Correlation Heatmap between Stress and Background Conditions",
        xlab = "PFAM IDs (Stress)",
        ylab = "PFAM IDs (Background)")