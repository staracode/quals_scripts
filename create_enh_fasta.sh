for stage in 'ESC' 'MES' 'CP' 'CM'
do 
	echo $stage
	for x in $(seq 1 19) 
		do 
		chrom='chr'${x}
		echo $chrom
		awk -v chr=$chrom '{ if ( $1 == chr ) { print  $1 "\t" $2 "\t" $3 "\t" chr "_" $4 } }' co/cardiac/k27ac_regions/k27ac_${stage}_regions.bed > co/cardiac/k27ac_regions/k27ac_${stage}_regions_${chrom}.bed
		fastaFromBed  -fi  genomes/mm9/${chrom}.fa -bed co/cardiac/k27ac_regions/k27ac_${stage}_regions_${chrom}.bed -fo co/cardiac/k27ac_regions/k27ac_${stage}_regions_${chrom}.fa -name
		done
	cat co/cardiac/k27ac_regions/k27ac_${stage}_regions_chr*.fa >  co/cardiac/k27ac_regions/k27ac_${stage}_regions.fa
	cat co/cardiac/k27ac_regions/k27ac_${stage}_regions_chr*.bed >  co/cardiac/k27ac_regions/k27ac_${stage}_regions.bed
done
