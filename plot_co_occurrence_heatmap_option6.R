
#library(e1071)
#use hamming.distance() 

library(gplots)

#chrom = "chr10"
chrom = "chr19"

stage = "CM"

#motif_count = read.table(paste("table_s4_union_enhancer_ESC_", chrom, "_co_occurrence_counts.txt", sep=""), stringsAsFactors=F)
#tf = read.table(paste("table_s4_union_enhancer_ESC_", chrom, "_co_occurrence_counts_tfname_header.txt", sep=""), stringsAsFactors=F, fill=T)

motif_count = read.table(paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence.txt", sep=""), stringsAsFactors=F)
tf = read.table(paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence_tfname_header.txt", sep=""), stringsAsFactors=F, fill=T)
colnames(motif_count) = tf[2,]
rownames(motif_count) = tf[1,!is.na(tf[1,])]
motif_count_rmCol = motif_count[,which(colSums(motif_count) != 0)]
motif_count2 = cor(t(as.matrix(motif_count_rmCol/sum(motif_count_rmCol))), method=c("spearman"))
#there is no difference hamming distance calcualted on binary and numerical data 
#motif_count2 = hamming.distance(as.matrix(motif_count))
pdf(paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence.pdf", sep=""))
heatmap.2(motif_count2,scale="n",trace="n",col=topo.colors(100), cexRow =0.2, cexCol =0.2)
dev.off()

#permute for pvalue
N=1000
motif_perm = list()
motif_perm2 = list()
for (x in seq(1, N)) {
	motif_perm[[x]] = t(apply(motif_count_rmCol, 1, sample))
	motif_perm2[[x]] = cor(t(as.matrix( motif_perm[[x]]/sum(motif_perm[[x]]) )), method=c("spearman"))
}

row_names = rownames(motif_count_rmCol)
#initialize the list of vectors within a list 
for (x in seq(1,N) ) {
	values = list()
	for (index1 in 1:length(row_names)) {
		values[[index1]] = list()	
		for (index2 in 1:length(row_names)) {
			values[[ index1 ]] [[ index2 ]]  = vector()
		}
	}
}

#set value of each vector within each list
for (x in seq(1,N) ) {
	for (index1 in 1:length(row_names)) {
		for (index2 in 1:length(row_names)) {			
			values[[ index1 ]][[ index2 ]] = c( values[[ index1 ]] [[ index2 ]], motif_perm2[[x]][ row_names[index1], row_names[index2] ] )
		}
	}
}

store_rows = data.frame()
pdf(paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence_permutations.pdf", sep=""))
par(mfrow=c(2,4))
for (index1 in 1:length(row_names)) {
	for (index2 in 1:length(row_names)) {
		x = motif_count2[ row_names[index1], row_names[index2] ]
		set = values[[ index1 ]][[ index2 ]]
		plot(set, ylim=c(-1,1), main=paste(row_names[index1], "\n", row_names[index2], sep=""))
		points(x, col="red")
		count = length(which(abs(set) >= abs(x)))
		keep_row = cbind(row_names[index1], row_names[index2], count, N, count/N)
		store_rows = rbind(store_rows, keep_row)

}} 
dev.off()

write.table(store_rows, paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence_permute_option6.txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)



library(igraph)	
t = which(motif_count2 > 0.70 & lower.tri(motif_count2), arr.ind=TRUE)	
t.graph	= graph.data.frame(t,directed=F)	
t.names	<-	colnames(motif_count2)[as.numeric(V(t.graph)$name)]	
pdf(paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence_network.pdf", sep="") )
plot(t.graph, vertex.size=5, vertex.shape="circle", vertex.label.color="red",	vertex.label=t.names, vertex.label.cex=0.2,	edge.width=0.5 ,	 layout=layout.fruchterman.reingold)	
dev.off()

library(reshape)
write.table(melt(motif_count2), paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence_network_count.txt", sep="") , row.names =F, col.names=F, quote=F)

write.table(motif_count2, paste("table_s4_union_enhancer_", stage, "_", chrom, "_numerical_cooccurrence_matrix.xls", sep=""), sep="\t", row.names=T, col.names=T, quote=F)
