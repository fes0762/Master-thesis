# Read common_genes.csv file
data <- read.csv("common_genes.csv", header = TRUE)

# Select the relevant columns for correlation calculation
time_points <- c("Induced_5hrs", "Induced_9hrs", "Induced_12hrs")
locations <- c("PS_cEpiblast", "NP_Epiblast", "aAO")

# Extract common genes
common_genes <- data[, 3]

# Create an empty correlation matrix
cor_matrix <- matrix(0, nrow = length(common_genes), ncol = length(common_genes))
rownames(cor_matrix) <- colnames(cor_matrix) <- common_genes

# Calculate the pairwise correlations
for (i in 1:length(common_genes)) {
  gene1 <- as.numeric(data[data[, 3] %in% common_genes[i], c(time_points, locations)])
  
  for (j in 1:(i+1)) {
    gene2 <- as.numeric(data[data[, 3] %in% common_genes[j], c(time_points, locations)])
    
    # Check if both vectors are not empty and have non-zero standard deviation
    if (all(!is.na(gene1) & !is.na(gene2)) && sd(gene1) > 0 && sd(gene2) > 0) {
      cor_value <- cor(gene1, gene2, method = "pearson", use = "pairwise.complete.obs")
      cor_matrix[i, j] <- cor_value
      cor_matrix[j, i] <- cor_value
    }
  }
  
  # Assign the diagonal element (the correlation of the gene with itself)
  cor_matrix[i, i] <- 1
}

#Fixing the vector error
# Calculate the pairwise correlations
for (i in 1:length(common_genes)) {
  gene1 <- unlist(data[data[, 3] %in% common_genes[i], c(time_points, locations)])
  gene1 <- as.numeric(gene1)
  
  for (j in 1:length(common_genes)) {
    gene2 <- unlist(data[data[, 3] %in% common_genes[j], c(time_points, locations)])
    gene2 <- as.numeric(gene2)
    
    # Remove any rows with NA values from both vectors
    gene1_complete <- gene1[!is.na(gene1) & !is.na(gene2)]
    gene2_complete <- gene2[!is.na(gene1) & !is.na(gene2)]
    
    # Check if both vectors are not empty, have non-zero standard deviation, and do not contain NA values
    if (length(gene1_complete) > 0 && length(gene2_complete) > 0 && !any(is.na(gene1_complete)) && !any(is.na(gene2_complete)) && sd(gene1_complete) > 0 && sd(gene2_complete) > 0) {
      cor_value <- cor(gene1_complete, gene2_complete, method = "pearson", use = "pairwise.complete.obs")
      cor_matrix[i, j] <- cor_value
      cor_matrix[j, i] <- cor_value
    }
  }
  
  # Assign the diagonal element (the correlation of the gene with itself)
  cor_matrix[i, i] <- 1
}


# Print the correlation matrix
print(cor_matrix)


#Create csv file
write.csv(cor_matrix, file = "correlation_matrix.csv", row.names = TRUE)
    
# Create a copy of the correlation matrix
thresholded_matrix <- cor_matrix

# Set values below the threshold (0.7) to 0
thresholded_matrix[thresholded_matrix < 0.7] <- 0

# Set diagonal elements to 0 (correlation of a variable with itself)
diag(thresholded_matrix) <- 0

# Create a CSV file from the binary matrix
write.csv(thresholded_matrix, file = "thresholded_matrix.csv", row.names = TRUE)

