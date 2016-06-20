setwd("/home/arpit.mishra/Network_analysis/alloted_proj/GSE9006_RAW/unp_files")
gene_data = read.delim("Healthy_fold_change.txt")
x = gene_data[,1:24]
y = gene_data[,25:117]

n = nrow(x)

p.value_x <- rep(0,n)
for (i in 1:n){
  p.value_x[i] <- t.test(x[i,]) $ p.value
  
}


p.value_y <- rep(0,n)
for (i in 1:n){
  p.value_y[i] <- t.test(y[i,]) $ p.value
  
}



p.value_x_y <- rep(0,n)
for (i in 1:n){
  p.value_x_y[i] <- t.test(x[i,],y[i,]) $ p.value
  
}





df<- cbind(p.value_x,p.value_y,p.value_x_y)
colnames(df) <- c('x','y','x_y')
rownames(df) <- rownames(gene_data)

write.table(df, "p_value_healthy.txt", sep="\t", quote=F)
write.csv(df,"p_value_healthy.csv")

gene_id_yes <- NULL
gene_id_no <- NULL
for (i in 1:nrow(df)){
  if(min(df[i,3])>0.05){
    gene_id_no <- rbind(gene_id_no,i)
  }
  else{
    gene_id_yes <- rbind(gene_id_yes,i)
  }
}

df <- as.data.frame(df)

final_genes_expression <- NULL
final_genes_pvalue <- NULL
for (i in 1:nrow(gene_id_yes)){
  final_genes_pvalue <- rbind(final_genes_pvalue,df[gene_id_yes[i],])
  final_genes_expression <- rbind(final_genes_expression,gene_data[gene_id_yes[i],])
}

write.table(final_genes_pvalue, "final_genes_p_value_healthy.txt", sep="\t", quote=F)
write.csv(final_genes_pvalue,"final_genes_p_value_healthy.csv")

write.table(final_genes_expression, "final_genes_expr_value_healthy.txt", sep="\t", quote=F)

nba_heatmap <- heatmap(as.matrix(final_genes_pvalue)[1:50,], Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
