setwd("?")
# Set the Working Directory
gene_new_o <- read.delim("gene_new_o.txt", header = FALSE)
result.1_p2 <- read.table("result-1_p2.txt", quote = "\"", comment.char = "")
source('fun_SPCA.R')
source('fun_ESPCA.R')
n = 500
myString <- "start!"
print(myString)
edges = list()
gene1 = data.matrix(gene_new_o)
gene = t(gene1)
#print(gene)
for (i in 1:1043168){
  edges[[i]] = c((result.1_p2[i,1]),(result.1_p2[i,2]))
  #print(edges[[i]])
}
out3 =  ESPCA(gene, k = 2, edges, k.group=25,we=0, t = 0.1,niter=10)
myString <- "end!"
#print(myString)

# Save U matrix to txt or csv file
write.table(out3$U, file = "U_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# If you want to save U matrix to csv file, use this line instead
# write.csv(out3$U, file = "U_matrix.csv", row.names = FALSE)

# Save V matrix to txt or csv file
write.table(out3$V, file = "V_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# If you want to save V matrix to csv file, use this line instead
# write.csv(out3$V, file = "V_matrix.csv", row.names = FALSE)

write.table(out3$D, file = "D_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)

