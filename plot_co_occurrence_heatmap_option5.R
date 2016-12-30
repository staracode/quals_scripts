
library(e1071)
#use hamming.distance() 

library(gplots)

#chrom = "chr10"
chrom = "chr19"

stage = "CM"

#motif_present_enhancer = read.table(paste("table_s4_union_enhancer_ESC_", chrom, "_motif_present_by_enhancer.txt", sep=""), stringsAsFactors=F)
#tf = read.table(paste("table_s4_union_enhancer_ESC_", chrom, "_motif_present_by_enhancer_tfname_header.txt", sep=""), stringsAsFactors=F, fill=T)
motif_present_enhancer = read.table(paste("table_s4_union_enhancer_", stage, "_", chrom, "_binary_cooccurrence.txt", sep=""), stringsAsFactors=F)
tf = read.table(paste("table_s4_union_enhancer_", stage, "_", chrom, "_binary_cooccurrence_tfname_header.txt", sep=""), stringsAsFactors=F, fill=T)
colnames(motif_present_enhancer) = tf[2,]
rownames(motif_present_enhancer) = tf[1,!is.na(tf[1,])]
motif_present_enhancer2 = hamming.distance(as.matrix(motif_present_enhancer))
pdf(paste("table_s4_union_enhancer_", stage, "_", chrom, "_binary_cooccurrence.pdf", sep=""))
heatmap.2(motif_present_enhancer2,scale="n",trace="n",col=topo.colors(1000), cexRow =0.2, cexCol =0.2)
dev.off()

#permute for pvalue
N=1000
motif_perm = list()
motif_perm2 = list()
for (x in seq(1, N)) {
	#motif_perm[[x]]=motif_present_enhancer[sample(nrow(motif_present_enhancer)),]
	#motif_perm[[x]]=motif_present_enhancer
	#rownames(motif_perm[[x]]) = sample (row_names)
	motif_perm[[x]] = t(apply(motif_present_enhancer, 1, sample))
	motif_perm2[[x]] = hamming.distance(as.matrix(motif_perm[[x]]))
}


row_names = rownames(motif_present_enhancer)
#initialize the list of vectors within a list 
for (x in seq(1,N) ) {
	values = list()
	for (index1 in 1:length(row_names)) {
		values[[index1]] = list()	
		for (index2 in 1:length(row_names)) {
			values[[ index1 ]] [[ index2 ]]  = vector()
			#values[[ row_names[index1] ]][[ row_names[index2] ]] = c( values[[ row_names[index1] ]] [[ row_names[index2] ]], motif_perm2[[x]][ row_names[index1], row_names[index2] ] )
			#print (motif_perm2[[x]][ row_names[index1], row_names[index2] ])
		}
	}
}

for (x in seq(1,N) ) {
	for (index1 in 1:length(row_names)) {
		for (index2 in 1:length(row_names)) {			
			values[[ index1 ]][[ index2 ]] = c( values[[ index1 ]] [[ index2 ]], motif_perm2[[x]][ row_names[index1], row_names[index2] ] )
		}
	}
}

store_rows = data.frame()
for (index1 in 1:length(row_names)) {
	for (index2 in 1:length(row_names)) {
		x = motif_present_enhancer2[ row_names[index1], row_names[index2] ]
		set = values[[ index1 ]][[ index2 ]]
		count = sum(set>=x)
		#print (set)
		#print (x)
		#print (count)
		#print (paste(row_names[index1], row_names[index2], count/N, sep=" "))
		keep_row = cbind(row_names[index1], row_names[index2], count, N, count/N)
		store_rows = rbind(store_rows, keep_row)
	}
}

write.table(store_rows, paste("table_s4_union_enhancer_", stage, "_", chrom, "_binary_cooccurrence_permute_option5.txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
write.table(motif_present_enhancer2, paste("table_s4_union_enhancer_", stage, "_", chrom, "_binary_cooccurrence_matrix.xls", sep=""), sep="\t", quote=F, col.names=T, row.names=T)
