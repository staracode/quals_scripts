import sys

combo_tf = sys.argv[1]
combinatorial = sys.argv[2]

link_tf_to_name = "/data/transfac/1_transfac/transfac_2012.4/match/data/mappingsTransfac.txt"
profile = "/data/transfac/1_transfac/transfac_2012.4/match/data/prfs/vertebrate_non_redundant_minSUM.prf"
#profile = "/data/transfac/1_transfac/transfac_2012.4/match/data/minSUM_good.prf"
#profile = "/data/transfac/1_transfac/transfac_2012.4/match/data/prfs/vertebrates.prf"

#map tf motif to gene name
gene_to_tf = {}
link_tf_to_name_handle = open(link_tf_to_name, "r")
for line in link_tf_to_name_handle: 
	fields = line.rstrip("\n").split()
	tf = fields[0].upper()
	gene = fields[2].upper()
	if gene_to_tf.has_key(gene) :
		gene_to_tf[gene].append(tf)
	else:
		gene_to_tf[gene] = [tf]
	
#convert expressed genes to TF to non-redundant TF name
redundant_to_non_redundant_dict = {}
#redundant_to_non_redundant = "redundant_to_non_redundant_key_threshold.txt"
redundant_to_non_redundant = "/home/tara/co/cardiac/link_vertebrate_redundant_to_non_redundant_motifs/redundant_to_non_redundant_key_threshold.txt"
redundant_to_non_redundant_handle = open(redundant_to_non_redundant, "r")
for line in redundant_to_non_redundant_handle:
	if line.startswith("#"): continue
	fields = line.rstrip(" \n" ).split ()
	redundant = fields[0].replace("_", "$", 1)
	non_redundant = fields[1].replace("_", "$", 1)
	redundant_to_non_redundant_dict[redundant] = non_redundant 
	#print redundant, non_redundant

#list of expressed genes 
combo_tf_handle = open (combo_tf, "r")
combo_tfs_motifs = []
combo_tfs_gene_names = []
for line in combo_tf_handle: 
	if  line.startswith("Gene") : continue
	fields = line.rstrip("\n").split("\t")
	gene = fields[0].upper()
	if gene_to_tf.has_key(gene):
		tf_list = gene_to_tf[gene]
		#print "found", gene
		#convert tf to non-redundant tf name
		for tf in tf_list:
			#print gene, tf
			if tf in redundant_to_non_redundant_dict.keys():
				non_redundant_tf = redundant_to_non_redundant_dict[tf]
				combo_tfs_motifs.append(non_redundant_tf)
				#print gene, tf, non_redundant_tf
			else:
				combo_tfs_motifs.append(tf)
	#else: 
		#print "could not match motif with gene name (hopefully it is not a tf)", gene
	combo_tfs_gene_names.append(gene)
combo_tf_handle.close()

#write prf output to :
combinatorial_handle = open ( combinatorial, "w")
#filter the prf file for expressed genes using keys  :  
profile_handle = open ( profile, "r")
for line in profile_handle: 
	fields = line.rstrip(" \n" ).split ()  
	if len(fields) != 5: 
		combinatorial_handle.write( line )
	 	continue
	tf_motif = fields[4]
	tf_gene_name = fields[4].replace("V$", "").split("_")[0]
	if tf_motif in combo_tfs_motifs: 
		combinatorial_handle.write( line )
	#if exact motif name not found, filter using gene name
	elif tf_gene_name in combo_tfs_gene_names:
		combinatorial_handle.write( line )


profile_handle.close()
combinatorial_handle.close()
