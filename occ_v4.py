import sys
sys.path.append('/home/tara/bin/mine')
from network_classes import * 
from maf_functions import *
import numpy as np
import os.path
import random

###########################################################################################################################################################
#current tf to focus on
tf_name = sys.argv[1]
tf = "V$" + tf_name
#current chrom
chrom = "chr19"
#curr stage
stage = "CM"

folder1 = "co/cardiac/" 
folder2 = "co/cardiac/version5/"
folder3 = "co/cardiac/version5/match_bg_sets/"
folder4 = "co/cardiac/version6/"
folder5 = "genomes"

##############################################
#enhancers
##############################################
#enhancers
#bed file containing enhancers
bed = folder1 + "/table_s4_union_enhancer_" + stage + "_" + chrom + ".bed"
#motif matches
motifFile = folder2 + "/table_s4_union_enhancer_" + stage + "_" + chrom + ".match"
#motif matches in which non-coding and repeat regions have been filtered out
#enh_motif_mean_sd = "final.xls"
enh_motif_mean_sd = folder2 + "/table_s4_union_enhancer_" + stage + "_" + chrom + "_motif_mean_sd.txt" 
print bed, motifFile, enh_motif_mean_sd

##############################################
#chrom 
##############################################
#size of current chrom
genome_size_file=folder5 + "/mm9_" + chrom + ".sizes"
#fasta of current chrom
#fasta = "genomes/mm9/" + chrom + ".fa"
#non-coding, non-repetitive roi
keepRegions = folder1 + "regionsToKeep.bed"
motifFile2 = folder4 + "/motifs_" + stage + "_" + chrom + ".match"
motifFile3 = folder4 + "/motifs_" + stage + "_" + chrom + "_filtered_regions.match"
print genome_size_file, keepRegions, motifFile2, motifFile3

##############################################

##############################################
#output
##############################################
#tf_name = tf.replace("V$", "")
outputResultsFile = folder3 + "/table_s4_union_enhancer_" + stage + "_" + chrom + "_" +  tf_name + "_results.txt"
outputResultsFile2 = folder3 + "/table_s4_union_enhancer_" + stage + "_" + chrom + "_" + tf_name + "_length_pairs_not_found.txt"
print outputResultsFile, outputResultsFile2
##############################################

#if option set to 0:  re-create file containing  motif mean occurrences and sd in enhancers
#otherwise read file containing motif mean and sd
option = 1

#recreate file of tf motifs not found in coding or gap regions
filterMatch = 0 #read in filtered file
#filterMatch = 1 #create filtered motif file

#read fewer tf motifs from .match file if True
debug="False" #read in all lines of motif file
#debug="True" #read a few lines

#maximum number of tf-tf-length matched sets
max_bg_match = 10


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

#pick random location from roi that are of a particular length
def pickRandomRegion(regionsOfInterest3, length): 
 	randStart = int("".join(map(str, random.sample(regionsOfInterest3, 1))))
	randEnd = randStart + length
	if randEnd in regionsOfInterest3:
		return (randStart, randEnd)
	#print length, randStart, randEnd
	return pickRandomRegion(regionsOfInterest3, length)

#identify whether to continue searching for tf-tf-length matches regions
def tf_tf_pairs_length_matched ( tf_tf_pairs, tf, length, max_bg_match ):
	#done searching for matched length region if done remains equivalent to 1
	done=1
	for tf2 in tf_tf_pairs[tf].keys():
		if tf2 in tf_motifs_within_roi.keys():
			count = tf_tf_pairs[tf][tf2][length]
			if count < max_bg_match:
				done = 0
				return done

#returns list of tf that could not be pairs to another tf and size matched
def tf_tf_pairs_length_return_unmatched ( tf_tf_pairs, tf, length, max_bg_match ):
	#done searching for matched length region if done remains equivalent to 1
	done=1
	keep=[]
	for tf2 in tf_tf_pairs[tf].keys():
		if tf2 in tf_motifs_within_roi.keys():
			count = tf_tf_pairs[tf][tf2][length]
			if count < max_bg_match:
				done = 0
				keep.append(tf2)
				return keep

###########################################################################################################################################################
#enhancer set
###########################################################################################################################################################
#retrieve distributions of enhancer lengths
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

#create file if it does not exist (if option == 0)
#otherwise read a file containing mean motif occurrences from enhancer set
if option == 0: 
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
	with open(enh_motif_mean_sd, 'wb') as f:
		#np.savetxt(f, mean_motif_occur, delimiter="\t")
		f.write("\t".join(unique_all_seen_tfnames))
		f.write("\n")
		f.write("\t".join(map (str, mean_motif_occur.tolist() ) ))
		f.write("\n")
		f.write("\t".join(map (str, sd_motif_occur.tolist() ) ))
		f.write("\n")
	f.close()

elif option == 1:
	f = open(enh_motif_mean_sd, "r")
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
print chrom
chromSizeHandle = open(genome_size_file, "r")
( chrom, end) = chromSizeHandle.readline().split("\t")
chromSizeHandle.close()
name = chrom
start = 0
regionsOfInterest2 = [region(name, chrom, start, int(end), "+")]

#create .match file that filters out motifs within repeats/codign regions (definined by regionsOfInterest3/4 )
if filterMatch==1:
	#regions that exlcude coding and gaps
	regionsToKeep = open(keepRegions, "r")
	regionsOfInterest3 = []
	for line in regionsToKeep: 
		fields = line.rstrip("\n").split("\t")
		chrom1 = fields[0]
		if chrom1 == chrom:
			start = int(fields[1])
			end = int(fields[2])
			regionsOfInterest3.extend(range(start, end))
	regionsOfInterest4 = set(regionsOfInterest3)

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
	motifFile3 = folder4 + "/motifs_" + stage + "_" + chrom + "_filtered_regions.match"
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
			else: count =  count+1
	motifHandle3.close() 

outputResults = open(outputResultsFile, "w")
outputResults2 = open(outputResultsFile2, "w")
outputResults.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("inc", "tf", "tf2", "length", "chrom", "regStart", "regEnd", "motif1_count", "motif2_count", "tf1_mean", "tf2_mean", "tf1_sd", "tf2_sd" ) )
print "Looking for TF-TF background regions"
#pick random motif
#for tf in tf_tf_pairs.keys():
	#if len(tf_motifs_within_roi[tf]) == 0: continue
for length in bed_lengths: 
	for iter in range(1, bed_lengths[length] +1 ) : 
		done_with_length=0
		inc = 0
		#keep looking for regions containing tf with this length
		#while tf_tf_pairs[tf][tf2][length] < max_bg_match
		while done_with_length == 0:
			inc += 1
			print tf, length, done_with_length
	
			#pick random location for that motif stored in
			regionObject2 = regionsOfInterest2[0]
			randomMotifIndex = random.sample( tf_motifs_within_roi[tf], 1 )
			motifObject = regionObject2.motifs[randomMotifIndex[0]]
			motifStart = regionObject2.return_motif_start_position(motifObject)
			motifEnd = regionObject2.return_motif_end_position(motifObject)
			#randomPlacement of motif within length of region
			randStart = np.random.randint(0, (length - (motifEnd - motifStart)) )
			#assign genomic coordinates to roi
			regStart = motifStart - randStart 
			regEnd = regStart + length	
			print inc, motifStart, motifEnd, regStart, regEnd

			#look for other tf found within that region
			#count occurrences of each tf within this region and store 
			count_occurrences = {}
			for tf2 in tf_tf_pairs.keys():
				#don't bother if we've reach the max number of matched tf-tf-length regions
				if tf2 != tf and tf_tf_pairs[tf][tf2][length] == max_bg_match: continue
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
			if len(count_occurrences.keys()) < 10  and inc >= 1000: 
				#return tf-tf pairs that could not be size matched
				unmatched_tfs = tf_tf_pairs_length_return_unmatched ( tf_tf_pairs, tf, length, max_bg_match )
				print length, tf, ",".join(unmatched_tfs)
				outputResults2.write("%s\t%s\t%s\n" % ( length, tf, ",".join(unmatched_tfs) ) )
				break
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
					outputResults.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (inc, tf, tf2, length, chrom, regStart, regEnd, motif1_count, motif2_count, tf1_mean, tf2_mean, tf1_sd, tf2_sd))
			#check to see if matched length regions have been found for all tf-tf pairs
 			done_with_length = 	tf_tf_pairs_length_matched ( tf_tf_pairs, tf, length, max_bg_match )
	
outputResults.close()
outputResults2.close()
