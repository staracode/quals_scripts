import sys


f = open("final.xls", "r")
line1 = f.readline() 
tf_motifs = line1.rstrip("\n").split("\t")

tf_tf_pairs = {}
for tf1 in tf_motifs:
	tf_tf_pairs[tf1] = {}
	for tf2 in tf_motifs:
		tf_tf_pairs[tf1][tf2] = 0

tf_tf_pairs2 = {}
for tf1 in tf_motifs:
	tf_tf_pairs2[tf1] = {}
	for tf2 in tf_motifs:
		tf_tf_pairs2[tf1][tf2] = {}

motifs_by_enh = {}

fileName = "output1.txt"
fileHandle = open(fileName, "r")

for line in fileHandle:
	fields= line.rstrip("\n").split("\t")
	enhId = fields[1]
	length = fields[5]
	length_iter = fields[6]
	motif = fields[7]
	key = enhId + "|" + length + "|" + length_iter 
	"""
	if motifs_by_enh.has_key(enhId): 
		motifs_by_enh[enhId].append(motif)
	else: 
		motifs_by_enh[enhId] = [motif]
	"""
	if motifs_by_enh.has_key(key): 
		motifs_by_enh[key].append(motif)
	else: 
		motifs_by_enh[key] = [motif]

#tf_tf_pairs is symmetric
for key in motifs_by_enh.keys():
	motifs = motifs_by_enh[key]
	key_fields = key.split("|")
	regId = key_fields[1] + "-" + key_fields[2]

	for motif1 in motifs: 
			for motif2 in motifs: 
				if tf_tf_pairs.has_key(motif1): 
					if tf_tf_pairs[motif1].has_key(motif2): 
						count = tf_tf_pairs[motif1][motif2]
						tf_tf_pairs[motif1][motif2] = count + 1
						if tf_tf_pairs2[motif1][motif2].has_key(regId):
							count2 = tf_tf_pairs2[motif1][motif2][regId]
							tf_tf_pairs2[motif1][motif2][regId] = count2 + 1
						else: 
							tf_tf_pairs2[motif1][motif2][regId] = 1
					'''
					else: 
						tf_tf_pairs[motif1][motif2] = 1
						tf_tf_pairs2[motif1][motif2] = {}
				else: 
					tf_tf_pairs[motif1] = {}
					tf_tf_pairs2[motif1] = {}
				'''

				#reverse, motif2 before motif1
				if tf_tf_pairs.has_key(motif2):
					if tf_tf_pairs[motif2].has_key(motif1):
						count = tf_tf_pairs[motif2][motif1]
						tf_tf_pairs[motif2][motif1] = count + 1
						if tf_tf_pairs2[motif2][motif1].has_key(regId):
							count2 = tf_tf_pairs2[motif2][motif1][regId]
							tf_tf_pairs2[motif2][motif1][regId] = count2 + 1
						else: 
							tf_tf_pairs2[motif2][motif1][regId] = 1
					'''
					else: 
						tf_tf_pairs[motif2][motif1] = 1
						tf_tf_pairs2[motif2][motif1] = {}
				else: 
					tf_tf_pairs[motif2] = {}
					tf_tf_pairs2[motif2] = {}
				'''

output1 = open("random_set_counts1.xls", "w")
output2 = open("random_set_counts2.xls", "w")
output1.write("%s\t%s\t%s\n" %("tf1", "tf2", "occurrences"))
output2.write("%s\t%s\t%s\t%s\n" %("tf1", "tf2", "matched_length_id", "occurrences"))
for motif1 in tf_tf_pairs.keys():
	for motif2 in tf_tf_pairs[motif1].keys():
		output1.write("%s\t%s\t%s\n" % (motif1, motif2, tf_tf_pairs[motif1][motif2]/2 ))
		for regId in tf_tf_pairs2[motif1][motif2].keys():
			output2.write( "%s\t%s\t%s\t%s\n" % (motif1, motif2, regId, tf_tf_pairs2[motif1][motif2][regId]/2))

output1.close()
output2.close()
