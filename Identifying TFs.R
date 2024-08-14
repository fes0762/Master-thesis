# Read the CSV files into RStudio
  protein_df <- read.csv("protein_class_transcription.csv", header = TRUE)
fpkm_df <- read.csv("FPKM_dataset.csv", header = TRUE)

# Extract the relevant columns
fpkm_genes <- fpkm_df[, 3]
protein_genes_col1 <- unique(protein_df[, 1])
protein_genes_col2 <- unlist(strsplit(as.character(protein_df[, 2]), ","))

# Remove leading/trailing spaces from gene synonyms
protein_genes_col2 <- trimws(protein_genes_col2)

# Combine the gene names from both columns
protein_genes <- unique(c(protein_genes_col1, protein_genes_col2))

# Find the intersection of genes
common_genes <- intersect(fpkm_genes, protein_genes)

# Print the common genes
print(common_genes)


