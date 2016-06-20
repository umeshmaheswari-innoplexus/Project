setwd("/home/arpit.mishra/Network_analysis/alloted_proj/GSE9006_RAW/unp_files")
name = read.csv("sample.csv")
gene_data = read.delim("filtered_genes_t_test.txt")

names_b <- name[grep("GSM25",name$Accession),]


healthy_names <- names_b[grep("Healthy",names_b$Title),][,1]
rem_names <- names_b[!grepl("Healthy",names_b$Title),][,1]

# Initially it is frame cannot be used with which
healthy_names <- as.vector(healthy_names)
rem_names <- as.vector(rem_names)

#Taking specific columns only
#substr used as names is gene_id.CEl
health_gene_data <- gene_data[, which(substr(names(gene_data),1,9) %in% healthy_names)]
rem_gene_data <- gene_data[,which(substr(names(gene_data),1,9) %in% rem_names                            )]
dim(health_gene_data)
dim(rem_gene_data)

cl_value <- rep(0,ncol(health_gene_data))
health_gene_data <- rbind(health_gene_data,cl_value)

cl_value <- rep(1,ncol(rem_gene_data))
rem_gene_data <- rbind(rem_gene_data,cl_value)

gene_data_grouped <- cbind(health_gene_data,rem_gene_data)

group <- as.numeric(gene_data_grouped[16134,])

gene_data_fold <- gene_data_grouped[1:16133,]

library(edgeR)

cds <- DGEList( gene_data_fold , group = group )

cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
dim( cds )

cds <- calcNormFactors( cds )

plotMDS( cds , main = "MDS Plot for Gene Data", labels = colnames( cds$counts ) )

cds <- estimateCommonDisp( cds )
