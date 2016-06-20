setwd("/home/arpit.mishra/Network_analysis/alloted_proj/GSE9006_RAW/unp_files")
name = read.csv("sample.csv")
gene_data = read.delim("filtered_genes_t_test.txt")

names_b <- name[grep("GSM25",name$Accession),]

group_one <- names_b[grep("1MO",names_b$Title),]
group_two <- names_b[grep("4MO",names_b$Title),]

ids <- names_b[which(substr(names_b$Title,23,28) %in% substr(group_two$Title,23,28)),] 
group_zero <- as.vector(ids[grep("New",ids$Title),][,1])
group_one <- as.vector(group_one[,1])
group_two <- as.vector(group_two[,1])


group_zero_data <- gene_data[, which(substr(names(gene_data),1,9) %in% group_zero)]
group_one_data <- gene_data[, which(substr(names(gene_data),1,9) %in% group_one)]
group_two_data <- gene_data[, which(substr(names(gene_data),1,9) %in% group_two)]

diff_zero_one <- group_zero_data-group_one_data
diff_zero_two <- group_zero_data-group_two_data
diff_one_two <- group_one_data-group_two_data

gene_data_fold <- cbind(group_zero_data,group_one_data,group_two_data)

n <- nrow(diff_one_two)
print(n)

for (i in 1:n){
  if(max(abs(diff_one_two[i,]))<1.5&&max(abs(diff_zero_one[i,]))<1.5&&max(abs(diff_zero_two[i,]))<1.5){
    gene_data_fold <- gene_data_fold[-i,]
  }
    
}

write.table(gene_data_fold, file="Month_fold_change.txt", quote=F, sep="\t")

