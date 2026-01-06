# ------------------------------------------------------------
# Preparation and setup
# ------------------------------------------------------------

# Set the working directory (optional)
# setwd("path/to/R_for_GSN-SPCA")

# Load required packages
# install.packages("reticulate")

# This script invokes Python.
# It is recommended to use the same Python environment as for GSN-SPCA.
# The environment configuration can be modified in fun_GSN-SPCA.R.

# ------------------------------------------------------------
# Remove the result_max_values.txt file if it exists
# ------------------------------------------------------------
file_path <- "result_max_values.txt"
if (file.exists(file_path)) {
  file.remove(file_path)
  print(paste("The file has been deleted: ", file_path))
} else {
  print(paste("The file does not exist, no deletion required: ", file_path))
}

# ------------------------------------------------------------
# Load input data
# ------------------------------------------------------------
gene_new_o <- read.delim("gene_new_o.txt", header = FALSE) 
result.1_p2 <- read.table("result-1_p2.txt", sep = "\t", quote = "", comment.char = "")

# ------------------------------------------------------------
# Load required functions
# ------------------------------------------------------------
source('fun_GSN-SPCA.R', local = TRUE)

# ------------------------------------------------------------
# Run ESPCA analysis
# ------------------------------------------------------------
n = 500
myString <- "start!"
print(myString)

edges = list()

gene1 = data.matrix(gene_new_o)
gene = t(gene1)   # transpose the gene data matrix

# Convert result.1_p2 into edge list
for (i in 1:nrow(result.1_p2)) {
  edges[[i]] = result.1_p2[i, ]
}

# Run ESPCA
out3 = ESPCA(gene, k = 4, edges, k.group = 80, niter = 20, we = 0, t = 0.1, err = 0.005)

myString <- "end!"
print(myString)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

# Save U matrix (components for samples)
write.table(out3$U, file = "U_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# Alternatively, to save as CSV:
# write.csv(out3$U, file = "U_matrix.csv", row.names = FALSE)

# Save V matrix (components for features)
write.table(out3$V, file = "V_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# Alternatively, to save as CSV:
# write.csv(out3$V, file = "V_matrix.csv", row.names = FALSE)

# Save D matrix (singular values / eigenvalues)
write.table(out3$D, file = "D_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
