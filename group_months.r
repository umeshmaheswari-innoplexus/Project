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

cl_value <- rep(0,ncol(group_zero_data))
group_zero_data <- rbind(group_zero_data,cl_value)

cl_value <- rep(1,ncol(group_one_data))
group_one_data <- rbind(group_one_data,cl_value)

cl_value <- rep(2,ncol(group_two_data))
group_two_data <- rbind(group_two_data,cl_value)

gene_data_grouped <- cbind(group_zero_data,group_one_data,group_two_data)


cl <- as.numeric(gene_data_grouped[16134,])

gene_data_sam <- gene_data_grouped[1:16133,]
sam.out = sam(gene_data_sam,cl,rand=123,gene.names = dimnames(gene_data_sam)[[1]])
plot(sam.out)
identify(sam.out)
