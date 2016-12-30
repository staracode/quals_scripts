import sys
sys.path.append('/home/tara/bin/mine')
from network_classes import * 
from maf_functions import *
import numpy as np
import os.path
import random

###########################################################################################################################################################
chrom = "chr19"
stage = "CM"
fasta = "genomes/mm9/" + chrom + ".fa"
#option to re-create file containing  motif mean occurrences and sd in enhancers
option = 1
#debug="False"
debug="True"
#total_loops=100
max_bg_match = 100
outputResultsFile = "output2.txt"
#outputResults2File = "output2.txt"
###########################################################################################################################################################

###########################################################################################################################################################
#file input
###########################################################################################################################################################

###########################################################################################################################################################
#functions
###########################################################################################################################################################
#remove duplicates from sorted list
def uniq(lst):
    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item

#sort list and  removed duplicates from list
def sort_and_deduplicate(l):
    return list(uniq(sorted(l, reverse=True)))

def pickRandomRegion(regionsOfInterest3, length): 
	randStart = int("".join(map(str, random.sample(regionsOfInterest3, 1))))
	randEnd = randStart + length
	if randEnd in regionsOfInterest3:
		return (randStart, randEnd)
	#print length, randStart, randEnd
	return pickRandomRegion(regionsOfInterest3, length)

###########################################################################################################################################################
#enhancer set
###########################################################################################################################################################

#bed file containing enhancers
#retrieve distributions of enhancer lengths
bed = "co/cardiac/table_s4_union_enhancer_CM_chr19.bed"
regionsOfInterest=[]
bedHandle = open(bed, "r")
bed_lengths = {} 
for line in bedHandle:
	fields = line.rstrip("\n").split("\t")
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	name = fields[3]
	regionsOfInterest.append(region(name, chrom, start, end, "+"))

	#length of enhancers
	length = int(fields[2]) - int(fields[1])	
	if bed_lengths.has_key(length):
		count = bed_lengths[length]
		bed_lengths[length] = count + 1
	else: 
		bed_lengths[length] = 1

#create file if it does not exist (option == 1)
#otherwise read a file containing mean motif occurrences from enhancer set
if option == 0: 
	#retrieve motifs found within enhancers
	#move to version6
	motifFile = "co/cardiac/version5/table_s4_union_enhancer_" + stage + "_" + chrom + ".match"
	#motif positions in regions of interest	
	motifHandle = open(motifFile, "r")
	for line in motifHandle:
		if line.startswith("Inspecting sequence ID"):
			fields = line.rstrip("\n").split()
			name = fields[3]
			#new region 
			for index, regionObject in enumerate(regionsOfInterest):
				if regionObject.name == name: 
					index_roi = index
		elif line.startswith(" V$"):
			fields = line.rstrip("\n").split("|")
			tfname = fields[0].split()[0]
			position = int(fields[1].split()[0])
			strand = fields[1].split()[1]
			motifSeq = fields[4].split()[0]
			motifObject = motif(tfname, position, strand, motifSeq)
			regionsOfInterest[index_roi].add_motif(motifObject)
	
	#mean and sd motif occurrence across all enhancers 
	#list of all tfnames present in all reigons of interest
	all_seen_tfnames = []
	for regionObject in regionsOfInterest:
		row_tfnames = [motifObject.tfname for motifObject in regionObject.motifs]
		all_seen_tfnames.extend(row_tfnames)

	unique_all_seen_tfnames = sort_and_deduplicate( all_seen_tfnames)
	
	#matrix of all unique tf seen
	matrix  = np.zeros(shape=(len(unique_all_seen_tfnames),len(regionsOfInterest)) , dtype=np.int)

	#now count number of times tf occurs within a region of interest
	for index1,regionObject in enumerate(regionsOfInterest): 
		row_tfnames = [str(motifObject.tfname) for motifObject in regionObject.motifs]
		for tfname1 in row_tfnames: 
			index2 = [index for index,name in enumerate(unique_all_seen_tfnames) if str(name) == str(tfname1)]
			if len(index2) > 1: 
				print "this tf name should be unique"
				sys.exit()
			elif len(index2) == 0:
				print "tf not found"
				sys.exit()
			#print index1, index2, len(unique_all_seen_tfnames)
			count = matrix[np.array( index2 ), np.array( index1 )]
			count = count + 1
			matrix[np.array( index2 ), np.array(index1)] = count
	
	mean_motif_occur = np.mean(matrix, axis=1)
	sd_motif_occur = np.std(matrix, axis=1)

	#write to file
	with open('final.xls', 'wb') as f:
		#np.savetxt(f, mean_motif_occur, delimiter="\t")
		f.write("\t".join(unique_all_seen_tfnames))
		f.write("\n")
		f.write("\t".join(map (str, mean_motif_occur.tolist() ) ))
		f.write("\n")
		f.write("\t".join(map (str, sd_motif_occur.tolist() ) ))
		f.write("\n")
	f.close()

elif option == 1:
	f = open("final.xls", "r")
	line1 = f.readline() 
	tf_motifs = line1.rstrip("\n").split("\t")
	line2 = f.readline() 
	tf_motif_means = line2.rstrip("\n").split("\t")
	line3 = f.readline() 
	tf_motif_sds = line3.rstrip("\n").split("\t")



#keep track of null sets for each tf-tf pair
tf_tf_pairs = {}
for tf1 in tf_motifs:
	tf_tf_pairs[tf1] = {}
	for tf2 in tf_motifs:
		tf_tf_pairs[tf1][tf2] = {}
		for length in bed_lengths: 
			tf_tf_pairs[tf1][tf2][length] = 0

#indices of regionsOfInterest2 that fall within regionsOfInterest3 for each motif
tf_motifs_within_roi = {}
for tf1 in tf_motifs:
	tf_motifs_within_roi[tf1] = []

###########################################################################################################################################################
#chr19
###########################################################################################################################################################
print "chr19"
chromSizeHandle = open("genomes/mm9_chr19.sizes", "r")
( chrom, end) = chromSizeHandle.readline().split("\t")
chromSizeHandle.close()
name = "chr19"
start = 0
regionsOfInterest2 = [region(name, chrom, start, int(end), "+")]

filterMatch=0
if filterMatch==1:
	#regions that exlcude coding and gaps
	regionsToKeep = open("regionsToKeep.bed", "r")
	regionsOfInterest3 = []
	for line in regionsToKeep: 
		fields = line.rstrip("\n").split("\t")
		chrom1 = fields[0]
		if chrom1 == chrom:
			start = int(fields[1])
			end = int(fields[2])
			regionsOfInterest3.extend(range(start, end))
	regionsOfInterest4 = set(regionsOfInterest3)

	motifFile2 = "co/cardiac/version6/motifs_" + stage + "_" + chrom + ".match"
	motifFile3 = "co/cardiac/version6/motifs_" + stage + "_" + chrom + "_filtered_regions.match"
	motifHandle2 = open(motifFile2, "r")
	motifHandle3 = open(motifFile3, "w")
	count=0
	for line in motifHandle2:
		if line.startswith("Inspecting sequence ID"):
			fields = line.rstrip("\n").split()
			name = fields[3]
			#new region 
			for index, regionObject in enumerate(regionsOfInterest2):
				if regionObject.name == name: 
					index_roi = index
			motifHandle3.write(line)
		elif line.startswith(" V$"):
			fields = line.rstrip("\n").split("|")
			tfname = fields[0].split()[0]
			position = int(fields[1].split()[0])
			strand = fields[1].split()[1]
			motifSeq = fields[4].split()[0]
			if position in regionsOfInterest4:
				#print tfname, position
				motifObject = motif(tfname, position, strand, motifSeq)
				regionsOfInterest2[index_roi].add_motif(motifObject)
				motifHandle3.write(line)
		else: 
			motifHandle3.write(line)
		if debug== "True":
			if count ==1000000:break
			else: count = count+1
	motifHandle2.close()
	motifHandle3.close()

#filtered match file already created
if filterMatch==0:
	motifFile3 = "co/cardiac/version6/motifs_" + stage + "_" + chrom + "_filtered_regions.match"
	motifHandle3 = open(motifFile3, "r")
	count=0
	for line in motifHandle3:
		if line.startswith("Inspecting sequence ID"):
			fields = line.rstrip("\n").split()
			name = fields[3]
			#new region 
			for index, regionObject in enumerate(regionsOfInterest2):
				if regionObject.name == name: 
					index_roi = index
			index2 = 0 #use this index to keep track of psotion in regionsOfInterest2
		elif line.startswith(" V$"):
			fields = line.rstrip("\n").split("|")
			tfname = fields[0].split()[0]
			position = int(fields[1].split()[0])
			strand = fields[1].split()[1]
			motifSeq = fields[4].split()[0]
			motifObject = motif(tfname, position, strand, motifSeq)
			regionsOfInterest2[index_roi].add_motif(motifObject)
			if tf_motifs_within_roi.has_key(tfname): 
				tf_motifs_within_roi[tfname].append(index2)
			index2 = index2 + 1 
			
		if debug== "True":
			if count ==1000000:break
			else: count = count+1
	motifHandle3.close()


def tf_tf_pairs_length_matched ( tf_tf_pairs, tf, length, max_bg_match ):
	done=0
	for tf2 in tf_tf_pairs[tf].keys()
		for length in tf_tf_pairs[tf][tf2].keys():
			count = tf_tf_pairs[tf][tf2][length]
			if count == max_bg_match:
				done = 1
				return done



#outputResults = open(outputResultsFile, "w")
print "Looking for TF-TF background regions"
#pick random motif
for tf in tf_tf_pairs.keys():
	if len(tf_motifs_within_roi[tf]) == 0: continue
	for length in bed_lengths: 
		#while tf_tf_pairs[tf][tf2][length] < max_bg_match
		done_with_length = 	tf_tf_pairs_length_matched ( tf_tf_pairs, tf, length, max_bg_match )
		if done_with_length == 1 : continue

		#keep looking for regions containing tf with this length
		#but give up after a certain amount of tries
		

		#
		#pick random location for that motif stored in
		regionObject2 = regionsOfInterest2[0]
		#randomMotifIndex = np.random.randint( range(0, len(tf_motifs_within_roi[tf] )))
		#motifObject = regionsOfInterest2.motifs[ tf_motifs_within_roi[tf][randomMotifIndex] ]
		randomMotifIndex = random.sample( tf_motifs_within_roi[tf], 1 )
		motifObject = regionObject2.motifs[randomMotifIndex[0]]
		motifStart = regionObject2.return_motif_start_position(motifObject)
		motifEnd = regionObject2.return_motif_end_position(motifObject)
		#randomPlacement of motif within length of region
		randStart = np.random.randint(0, length)
		#assign genomic coordinates to roi
		regStart = motifStart - randStart
		regEnd = regStart + length	

		#look for other tf found within that region
		#count occurrences of each tf within this region and store 
		count_occurrences = {}
		for tf2 in tf_tf_pairs.keys():
			#don't bother if we've reach the max number of matched tf-tf-length regions
			if tf_tf_pairs[tf][tf2][length] == max_bg_match: continue
			for index2 in tf_motifs_within_roi[tf2]:
				motifObject2 = regionObject2.motifs[ index2 ]
				motifStart2 = regionObject2.return_motif_start_position(motifObject2)
				motifEnd2 = regionObject2.return_motif_end_position(motifObject2)
				
				#does motif intersect roi
				overlap = intersectF( (chrom, regStart, regEnd, "+"), (chrom, motifStart2, motifEnd2, "+") ) 
				if overlap != 0: 
					#print tf, tf2, length, chrom, regStart, regEnd, chrom, motifStart, motifEnd, chrom, motifStart2, motifEnd2, overlap
					if count_occurrences.has_key(tf2):
						count_occurrences[tf2] += 1
					else: 
						count_occurrences[tf2] = 1

		#do the occurrences for each tf fall within two sd of the mean for each tf
		motif1_count = count_occurrences[tf]
		for tf2 in tf_tf_pairs.keys():
			if count_occurrences.has_key(tf2):
				motif2_count = count_occurrences[tf2]
			else: 
				#this tf was not seen in this roi
				continue
			#calculate number of occurrences for tf within roi
			index3 = tf_motifs.index(tf)
			tf1_mean = float(tf_motif_means[index3])
			tf1_sd = float(tf_motif_sds[index3])
			#calculate number of occurrences for other tf within roi
			index4 = tf_motifs.index(tf2)
			tf2_mean = float(tf_motif_means[index4])
			tf2_sd = float(tf_motif_sds[index4])
			#calculate zscore for each 
			zscore1 = (tf1_mean-motif1_count)/tf1_sd
			abs_zscore1 = abs(zscore1)
			zscore2 = (tf2_mean-motif2_count)/tf2_sd
			abs_zscore2 = abs(zscore2)
			#keep region if zscore less than 2
			if abs_zscore1 <= 2 and abs_zscore2 <=2:
				#increment  to keep track of this TF-TF-length matched region 
				tf_tf_pairs[tf][tf2][length] += 1
				print tf, tf2, length, chrom, regStart, regEnd, motif1_count, motif2_count, tf1_mean, tf2_mean, tf1_sd, tf2_sd

#outputResults.close()


'''
#create a dictionary for each tf listing out the indexes in regionsOfInterset2.motifs that contain motifs falling outside of excluded regions
regionObject2 = regionsOfInterest2[0]
for index, motifObject in enumerate(regionObject2.motifs): 
	motifStart = regionObject2.return_motif_start_position(motifObject)
	motifEnd = regionObject2.return_motif_end_position(motifObject)
	#check to see if motif falls within included roi (no gaps or coding regions)
	if motifStart in regionsOfInterest3 and motifEnd in regionsOfInterest3:
		key = motifObject.tfname
		#add to dictionary
		tf_motifs_within_roi[key].append(index)

store_tf_index_handle = open("final_store_tf_index.txt", "w")
for key in tf_motifs_within_roi.keys(): 
		store_tf_index_handle.write(key + "\t" + "\t".join(map(str,tf_motifs_within_roi[key])) + "\n")
store_tf_index_handle.close()
'''
'''
outputResults.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"  %("loop", "regId", "regChr", "regStart", "regEnd", "length", "length_iter", "motif", "motifChr", "motifStart", "motifEnd", "motif_occ", "enh_mean", "enh_sd", "zscore"))
enhancer_id = 1
#create intervals from bed lengths
for loop in range(1, total_loops) :
		for length in bed_lengths: 
			counts = bed_lengths[length]
			for count in range(1, counts): 
				#pick random length-matched intervalfrom starts (end = start + length)
				regionObject2 = regionsOfInterest2[0]
				#randStart = np.random.randint(regionObject2.start, regionObject2.end )
				(randStart, randEnd) = pickRandomRegion(regionsOfInterest3, length)
				#print "loop", length, randStart, randEnd

				#loop through motifs on chr and store motif occurrences that overlap with random region
				motif_occ = {}
				for motifObject in regionObject2.motifs: 
					motifStart = regionObject2.return_motif_start_position(motifObject)
					motifEnd = regionObject2.return_motif_end_position(motifObject)
					#count motif_occurrence for each random interval
					intersection = intersectF ( ("chr19", motifStart, motifEnd, "+") , ("chr19", randStart, randEnd, "+") )  
					if intersection != 0 : 
						#print randStart, randEnd, motifObject.tfname, motifStart, motifEnd
						#store motif occurrences
						if motif_occ.has_key(motifObject.tfname): 
							x = motif_occ[motifObject.tfname] 
							motif_occ[motifObject.tfname] = x + 1
						else:  
							motif_occ[motifObject.tfname] = 1
						#if occurrences falls within range (for all motifs), keep interval and motif occurrences
				
				for motif in motif_occ.keys():
					#print "%s\t%s\t%s" % (length, motif, motif_occ[motif] )
					motif_count = motif_occ[motif]
					if motif in tf_motifs: 
						index = tf_motifs.index(motiy6f)
						tf_mean = float(tf_motif_means[index])
						tf_sd = float(tf_motif_sds[index])
						zscore = (tf_mean-motif_count)/tf_sd
						abs_zscore = abs(zscore)
						#if motif_count falls within one std dev of the enhancer motif mean, keep
						if abs_zscore <= 2:
							outputResults.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %  (loop, enhancer_id, chrom, randStart, randEnd, length, count, motif, chrom, motifStart, motifEnd, motif_count, tf_mean, tf_sd, zscore))

				enhancer_id = enhancer_id + 1
'''

#stop when each tf-tf has 10 matched-length sets 
#for each tf-tf
#store coordinates and tf counts (minimum)

					
