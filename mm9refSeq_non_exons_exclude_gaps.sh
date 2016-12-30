for i in $(seq 1 10)  
	do
	lineEnd=$((255*$i))
	lineStart=$(($lineEnd-254))
	echo $lineStart, $lineEnd
	#sed -n "$lineStart , $lineEnd p" co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_chr19.bed >  co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.bed

	#shuffle regions
	shuffleBed -excl co/cardiac/version5/exclude_regions.bed  -i co/cardiac/table_s4_union_enhancer_CM_chr19.bed  -g genomes/mm9_chr19.sizes | awk 'BEGIN{count=0} {count=count+1; print $1"\t"$2"\t"$3"\t"count  }' > co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.bed 

	#fasta
	fastaFromBed  -fi  genomes/mm9/chr19.fa -bed co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.bed -fo co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.fa -name

	#search for expressed tf motifs
	/data/transfac/1_transfac/transfac_2012.4/match/bin/match /data/transfac/1_transfac/transfac_2012.4/match/data/matrix.dat  co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.fa co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.match co/cardiac/transfac_CM_TF.prf

	#occurrences in null set -binary
	python  co/scripts/runTFModels.py  co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.bed co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.match  5 co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}_binary_cooccurrence.txt

	#numerical occurrence
	python  co/scripts/runTFModels.py  co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.bed co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}.match 6 co/cardiac/version5/mm9refSeq_non_exons_exclude_gaps_CM_chr19_null_${i}_numerical_cooccurrence.txt


done
