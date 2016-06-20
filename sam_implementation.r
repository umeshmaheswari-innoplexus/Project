setwd("/home/arpit.mishra/Network_analysis/alloted_proj/GSE9006_RAW/unp_files")
gene_data = read.delim("filtered_genes_t_test.txt")
library(siggenes)
library(multtest)
library(splines)

cl = rep(1,ncol(gene_data))
View(cl)

gene.names = dimnames(gene_data)[[1]]

sam.out = sam(gene_data,cl,rand=123,gene.names = dimnames(gene_data)[[1]])
plot(sam.out)
0
sam.out = sam(gene_data,cl,method = d.stat,control = samControl(),gene.names = dimnames(gene_data)[[1]])
plot(sam.out,0.296)
identify(sam.out)
