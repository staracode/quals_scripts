import sys
combo_tf = sys.argv[1]
combinatorial = sys.argv[2]
profile = "/data/transfac/1_transfac/transfac_2012.4/match/data/minSUM_good.prf"

profile_handle = open ( profile, "r")
#combo_tf = "dev_gene_lists/ravasi_10_cell_tableS1.txt" 

combo_tf_handle  =open (combo_tf, "r")

combo_tfs = []
for line in combo_tf_handle: 
	fields = line.rstrip("\n").split("\t")
	tf = fields[0].upper()
	combo_tfs.append(tf)

combo_tf_handle.close()

#combinatorial = "minSUM_good_ravasi_10_cell_tableS1.prf"
combinatorial_handle = open ( combinatorial, "w")

for line in profile_handle: 
	fields = line.rstrip("\n").split()
	if len(fields) != 5: 
		combinatorial_handle.write( line)
		continue

	tf = fields[4].replace("V$", "").split("_")[0]
	if tf in combo_tfs: 
		combinatorial_handle.write( line)
	
combinatorial_handle.close()
profile_handle.close()
