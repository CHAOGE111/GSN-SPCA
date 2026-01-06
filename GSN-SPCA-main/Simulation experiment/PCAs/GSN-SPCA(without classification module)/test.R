# Set the Working Directory
# setwd("(your workng Directory ))")
#install.packages("reticulate")

file_path <- "result_max_values.txt"
if (file.exists(file_path)) {
  file.remove(file_path)
  print(paste("The file has been deleted: ", file_path))
} else {
  print(paste("The file does not exist. There is no need to delete it.: ", file_path))
}

gene_new_o <- read.delim("gene_new_o.txt", header = FALSE) 

result.1_p2 <- read.table("result-1_p2.txt", sep = "\t", quote = "", comment.char = "")

source('fun_GSN-SPCA.R', local = TRUE)

n = 500
myString <- "start!"
print(myString)

edges = list()

gene1 = data.matrix(gene_new_o)
gene = t(gene1)

for (i in 1:nrow(result.1_p2)) {
  edges[[i]] = result.1_p2[i, ]
  print(edges[[i]])
}

out3 =  ESPCA(gene, k = 2, edges, k.group=3, niter=30, we=0, t = 0.1)
myString <- "end!"
print(myString)

file_path <- "result_max_values.txt"

# Save U matrix to txt or csv file
write.table(out3$U, file = "U_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# If you want to save U matrix to csv file, use this line instead
# write.csv(out3$U, file = "U_matrix.csv", row.names = FALSE)

# Save V matrix to txt or csv file
write.table(out3$V, file = "V_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# If you want to save V matrix to csv file, use this line instead
# write.csv(out3$V, file = "V_matrix.csv", row.names = FALSE)

write.table(out3$D, file = "D_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)

