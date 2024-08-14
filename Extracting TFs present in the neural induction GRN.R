#Using the Trevers et al 2023 neural induction GRN as a seed list
#Write common interactions to a file
write.csv(common_interactions, "common_interactions.csv", row.names = FALSE)

#Split the interactions back into separate columns:
common_interactions_split <- common_interactions %>%
  separate(Interaction, into = c("Gene1", "Gene2"), sep = "-")

#Write the split common interactions to a file
write.csv(common_interactions_split, "common_interactions_split.csv", row.names = FALSE)

#Attempt to read GRN file as a text file
file_contents <- readLines("Suppl_File_S1_NI_GRN.btp")

#Print the first few lines to see the structure
print(head(file_contents))

#Print the number of lines in the file
print(paste("Number of lines:", length(file_contents)))

library(xml2)
library(dplyr)
biotapestry_data <- read_xml("Suppl_File_S1_NI_GRN.btp")

#Extract all gene nodes
genes <- xml_find_all(biotapestry_data, "//gene")

#Extract the names of the genes
gene_names <- xml_attr(genes, "name")
print(head(gene_names))

#Check if Gene1 and Gene2 are in the gene_names list
common_interactions_split <- common_interactions_split %>%
  mutate(Gene1_Highlight = ifelse(Gene1 %in% gene_names, "Yes", "No"),
         Gene2_Highlight = ifelse(Gene2 %in% gene_names, "Yes", "No"))
print(head(common_interactions_split))
write.csv(common_interactions_split, "common_interactions_highlighted.csv", row.names = FALSE)