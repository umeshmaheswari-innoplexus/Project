setwd("/home/arpit.mishra/Network_analysis/alloted_proj/GSE9006_RAW/unp_files")
library(affy)
affy.data = ReadAffy()
eset.mas5 = mas5(affy.data)
exprSet.nologs = exprs(eset.mas5)

dim(exprSet.nologs)

exprSet = log(exprSet.nologs,2)

write.table(exprSet, file="Matrix_after_log.txt", quote=F, sep="\t")

data.mas5calls = mas5calls(affy.data)

data.mas5calls.calls = exprs(data.mas5calls)
write.table(data.mas5calls.calls, file="AP_calls_Matrix.txt", quote=F, sep="\t")

AP = apply(data.mas5calls.calls, 1, paste, collapse="")

genes.present = names(AP[AP != "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"])

length(genes.present)

exprSet.present = exprSet[genes.present,]
dim(exprSet.present)
p.value.all.genes = apply(exprSet, 1, function(x) { t.test(x[1:117]) $p.value } )
p.value.all.genes[1:5]

raw.pvals.present = p.value.all.genes[genes.present]

fdr.pvals.present = p.adjust(raw.pvals.present, method="fdr")

fdr.pvals.present.sorted = fdr.pvals.present[order(fdr.pvals.present)]

View(raw.pvals.present)
View(fdr.pvals.present)
expression.plus.pvals = cbind(exprSet.present,raw.pvals.present,fdr.pvals.present)

dim(expression.plus.pvals)


write.table(expression.plus.pvals, "Su_mas5_DE_analysis.txt", sep="\t", quote=F)

DE.probesets = names(raw.pvals.present[raw.pvals.present < 0.05])


all.data = read.delim("Su_mas5_DE_analysis.txt")
dim(all.data)
all.data = all.data[,1:117]
dim(all.data)

filtered_genes = all.data[DE.probesets,]
dim(filtered_genes)

write.table(filtered_genes, "filtered_genes_t_test.txt", sep="\t", quote=F)







