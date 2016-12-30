import sys
import copy
sys.path.append('/home/tara/bin/mine')
#sys.path.append('/Users/ucsf/mine')
from maf_functions import *
from network_classes import * 
#import matplotlib.pyplot as plt
import numpy as np
import pickle
import os.path

############################################################################################
#functions
############################################################################################
def distances_between_motifs_within_regions_of_interest(regionsOfInterest):
	roi_cov_mat = []
	for regionObject in regionsOfInterest:
		#print regionObject.chrom, regionObject.start, regionObject.end, regionObject.strand
		#print regionObject.length
		#print regionObject.print_motifs()

		#co occuring motifs
		covariance_matrix = []
		for motifObject1 in regionObject.motifs:
			start1 = regionObject.return_motif_start_position(motifObject1)
			end1 = regionObject.return_motif_end_position(motifObject1)
			row=[]
			for motifObject2 in regionObject.motifs: 
				start2 = regionObject.return_motif_start_position(motifObject2)
				end2 = regionObject.return_motif_end_position(motifObject2)
				intersection = intersectF ( ("chr",start1, end1), ("chr", start2, end2))
				#motif positions along region of interest overlap if interseciton != 0
				if intersection == 0: 
					#first motif is to the right of the second motif
					if start1 > start2: 
						distance = start1 - end2
					#first motif is to the left of the second motif
					else: 
						distance = start2 - end1
				else: 
					#motifs overlap
					distance = 0
				row.append(distance)
			covariance_matrix.append(row) 
		roi_cov_mat.append( covariance_matrix )
	return roi_cov_mat
					
def distances_between_motifs_within_regions_of_interest_v2(regionsOfInterest):

	all_seen_tfnames = []
	for regionObject in regionsOfInterest:
		row_tfnames = [motifObject.tfname for motifObject in regionObject.motifs]
		all_seen_tfnames.extend(row_tfnames)

	unique_all_seen_tfnames = sort_and_deduplicate( all_seen_tfnames)

	#dictionary of all pairs of tf
	#that points to a list of distances between those two tf
	tf_pairs = {}
	for tf1 in unique_all_seen_tfnames: 
		for tf2 in unique_all_seen_tfnames: 
				key = tf1 + "_" + tf2
				key_rev = tf2 + "_" + tf1
				if not ( tf_pairs.has_key(key) or tf_pairs.has_key(key_rev) ) : 
					#create list
					tf_pairs[key] = []

	#roi_cov_mat = []
	for regionObject in regionsOfInterest:
		for motifObject1 in regionObject.motifs:
			tf1 = motifObject1.tfname
			start1 = regionObject.return_motif_start_position(motifObject1)
			end1 = regionObject.return_motif_end_position(motifObject1)
			row=[]
			for motifObject2 in regionObject.motifs: 
				tf2 = motifObject2.tfname
				start2 = regionObject.return_motif_start_position(motifObject2)
				end2 = regionObject.return_motif_end_position(motifObject2)
				intersection = intersectF ( ("chr",start1, end1), ("chr", start2, end2))
				#motif positions along region of interest overlap if interseciton != 0
				if intersection == 0: 
					#first motif is to the right of the second motif
					if start1 > start2: 
						distance = start1 - end2
					#first motif is to the left of the second motif
					else: 
						distance = start2 - end1
				else: 
					#motifs overlap
					distance = 0

				#append distances to list
				#each list is the set  of distances between two tfs across roi
				key = tf1 + "_" + tf2
				key_rev = tf2 + "_" + tf1
				if tf_pairs.has_key(key) and not tf_pairs.has_key(key_rev):
					tf_pairs[key].append(distance)
				elif tf_pairs.has_key(key_rev) and not tf_pairs.has_key(key): 
					tf_pairs[key_rev].append(distance)
				elif tf_pairs.has_key(key) and tf_pairs.has_key(key_rev): 
					if tf1 == tf2: 
						tf_pairs[key].append(distance)
					else: 
						print "tf pair has been duplicated"
						sys.exit()
				else: 
					print "tf pair key not found"
					sys.exit()
	return tf_pairs

def distances_between_motifs_within_regions_of_interest_v3(regionsOfInterest, listOfTFs):
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %("roi_ID", "roi_chrom", "roi_start", "roi_end", "tf1", "tf1_chrom", "tf1_start", "tf1_end", "tf1_position", "tf1_strand", "tf2", "tf2_chrom", "tf2_start", "tf2_end", "tf2_position", "tf2_strand", "dist")
	for regionObject in regionsOfInterest:
		for motifObject1 in regionObject.motifs:
			tf1 = motifObject1.tfname
			if tf1.replace("V$", "") in listOfTFs:
				start1 = regionObject.return_motif_start_position(motifObject1)
				end1 = regionObject.return_motif_end_position(motifObject1)
				for motifObject2 in regionObject.motifs:
					tf2 = motifObject2.tfname
					if tf2.replace("V$", "") in listOfTFs:
						start2 = regionObject.return_motif_start_position(motifObject2)
						end2 = regionObject.return_motif_end_position(motifObject2)
						intersection = intersectF ( ("chr",start1, end1), ("chr", start2, end2))
						#motif positions along region of interest overlap if interseciton != 0
						if intersection == 0:
							#first motif is to the right of the second motif
							if start1 > start2: 
								distance = start1 - end2
							#first motif is to the left of the second motif
							else: 
								distance = start2 - end1
						else: 
							#motifs overlap
							distance = 0
						#print regionObject.name, regionObject.chrom, regionObject.start, regionObject.end, tf1, chrom, start1, end1, motifObject1.position, motifObject1.strand, tf2, chrom,  start2, end2, motifObject2.position, motifObject2.strand, distance
						print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (regionObject.name, regionObject.chrom, regionObject.start, regionObject.end, tf1, chrom, start1, end1, motifObject1.position, motifObject1.strand, tf2, chrom,  start2, end2, motifObject2.position, motifObject2.strand, distance)


#return pairs of motifs that overlap 					
def distances_between_motifs_within_regions_of_interest_v4(regionsOfInterest):
	inc =0
	store_by_roi = {}
	for index0, regionObject in enumerate(regionsOfInterest):
		for index1, motifObject1 in enumerate(regionObject.motifs):
			inc = inc + 1
			if inc == 10000:
				return store_by_roi
			tf1 = motifObject1.tfname
			start1 = regionObject.return_motif_start_position(motifObject1)
			end1 = regionObject.return_motif_end_position(motifObject1)
			for index2, motifObject2 in enumerate(regionObject.motifs):
				if index1 != index2:
					tf2 = motifObject2.tfname
					start2 = regionObject.return_motif_start_position(motifObject2)
					end2 = regionObject.return_motif_end_position(motifObject2)
					intersection = intersectF ( ("chr",start1, end1), ("chr", start2, end2))
					#motif positions along region of interest overlap if interseciton != 0
					if intersection == 0:
						#first motif is to the right of the second motif
						if start1 > start2:
							distance = start1 - end2
						#first motif is to the left of the second motif
						else: 
							distance = start2 - end1
					else:
						#motifs overlap
						distance = 0
						#keep track of the fact that these two motifs interact
						store_pair = [index1, index2]
						if store_by_roi.has_key(index0): 
							store_by_roi[index0].append( store_pair )
						else:
							store_by_roi[index0] = [ store_pair ]
						#print tf1, tf2

	return store_by_roi
					
def locate_variants_within_motifs_within_hars(regionsOfInterest):
	print "\t".join( ["mismatch_position", "motif_start", "motif_end", "motif_strand", "motif_length", "motif_seq", "motif_tfname", "region_length", "region_name"] ) 
	for regionObject in regionsOfInterest: 
		for motifObject1 in regionObject.motifs: 
		 	for mismatchObject in regionObject.mismatches: 
				if  mismatchObject.position_along_region == regionObject.return_motif_start_position(motifObject1): 
		 			output =  [mismatchObject.position_along_region, regionObject.return_motif_start_position(motifObject1), regionObject.return_motif_end_position(motifObject1),  motifObject1.strand, motifObject1.motifLength, motifObject1.motifSeq, motifObject1.tfname, regionObject.length,  regionObject.name ]
					print "\t".join(map(str, output))  
				elif mismatchObject.position_along_region == regionObject.return_motif_end_position(motifObject1):
		 			output =  [mismatchObject.position_along_region, regionObject.return_motif_start_position(motifObject1), regionObject.return_motif_end_position(motifObject1),  motifObject1.strand, motifObject1.motifLength, motifObject1.motifSeq, motifObject1.tfname, regionObject.length,  regionObject.name ]
					print "\t".join(map(str, output))  
					
		 		else: 
					intersection = intersectF ( ("chr", mismatchObject.position_along_region, mismatchObject.position_along_region, "+"), ("chr", regionObject.return_motif_start_position(motifObject1), regionObject.return_motif_end_position(motifObject1), "+"))	
					if intersection != 0: 
		 				output =  [intersection[1], regionObject.return_motif_start_position(motifObject1), regionObject.return_motif_end_position(motifObject1), motifObject1.strand,  motifObject1.motifLength, motifObject1.motifSeq, motifObject1.tfname, regionObject.length,  regionObject.name ]
		 				#output =  [intersection[1], regionObject.return_motif_start_position(motifObject1), regionObject.return_motif_end_position(motifObject1),  motifObject1.motifLength, motifObject1.motifSeq, motifObject1.tfname, regionObject.length,  regionObject.name, motifObject1.position, motifObject1.strand, regionObject.start, regionObject.end ]
						print "\t".join(map(str, output))  

def print_matrix(matrix):
	print '\t', '\t'.join(map(str, range(0, len(matrix))))
	for i, row in enumerate(matrix):
		print "%s\t%s" % ( i, '\t'.join(map(str,row)))

def print_matrix2(matrix, row_col_names):
	print '\t', '\t'.join(row_col_names)
	for i, row in enumerate(matrix):
		print "%s\t%s" % ( row_col_names[i], '\t'.join(map(str,row)))

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

#old function, use merge_overlapping_motifs2
def merge_overlapping_motifs(overlapping_motifs):
	prev_len = len(overlapping_motifs)
	overlapping_motifs_new = {}
	for roi_index in overlapping_motifs.keys():
		overlapping_pairs = overlapping_motifs[roi_index]
		#print overlapping_pairs
		for pair1 in overlapping_pairs:
			if overlapping_motifs_new.has_key(roi_index):
				overlapping_motifs_new[roi_index].append(pair1)
			else:
				overlapping_motifs_new[roi_index] = [ pair1 ]
			for pair2 in overlapping_pairs:
				if pair1 == pair2: continue
				print roi_index, pair1, pair2
				for element1 in pair1: 
					if element1 in pair2:
						pair1_temp = set(pair1)
						pair2_temp = set(pair2)
						pair_merge = list(pair1_temp.union(pair2_temp))
						pair_merge2 = sort_and_deduplicate(pair_merge)
						if overlapping_motifs_new.has_key(roi_index):
							if overlapping_motifs_new[roi_index][-1] == pair1:
								check = overlapping_motifs_new[roi_index].pop()
								print "remove pair1: ", check
							overlapping_motifs_new[roi_index].append(pair_merge2)
							print "add merge: ", pair_merge2
					 	break
		print overlapping_motifs_new[roi_index]
		if prev_len == len( overlapping_motifs_new[roi_index]):
			return overlapping_motifs_new
		else: 
			merge_overlapping_motifs(overlapping_motifs_new)
		prev_len= len( overlapping_motifs_new[roi_index] )


def merge_overlapping_motifs2(overlapping_sets):
	print "++++++++++++++++++++++++++ ++++++++++++++++++++++++++++++++++++"
	prev_num_sets = len(overlapping_sets)
	print overlapping_sets
	#pair1 = A-C
	for pair1 in overlapping_sets:
		index_to_merge=[]  #store inde x of elements to mergem
		#A-C, A-D, A-E, B-F, C-D
		#[A-C, A-D, A-E, C-D]
		#delete A-C, A-D, A-E, C-D
		#delete 0,1,2,4 in list
		#add [A, C, D, E]
		for index2, pair2 in enumerate(overlapping_sets):
			for element1 in pair1: 
				if element1 in pair2:
					if index2 not in index_to_merge:
		 			 	index_to_merge.append(index2)
						print "overlap found" 
						print "p1, p2	 		" , pair1, pair2
					break #no need to l oop through pair1 elements if a single element is found in pair2, pair1 and pair2 should be merged 

		#alter copy of distionary of lists
		values_to_merge = [overlapping_sets[index3] for index3 in index_to_merge]
		print "values to merge		", values_to_merge
		print "index_to_merge		", index_to_merge
		values_to_merge2 = []
		
		#to delete and add new merged  element, the number of positions has to be greater than 1 
		#since pair1 overlaps with pair2 at least once
		if len(index_to_merge) > 1:
			for inc,  index3 in enumerate(index_to_merge):
				print "index3, inc		", index3, inc
				values_to_merge2.extend( overlapping_sets[index3-inc] )
				print "values to merge loop	", values_to_merge2
				check = overlapping_sets.pop(index3-inc)
				print "check:			", check
				if check not in values_to_merge:
				 	print "Error: ", check, " not in ", values_to_merge 
					sys.exit() 

			#remove duplicates in [A. C. A. D. A. E. C. D]
			values_to_merge3 = sort_and_deduplicate(values_to_merge2)
			print "after", values_to_merge3
			#add [A, C, D, E]
			overlapping_sets.append(values_to_merge3)
		
			#remove empty lists
			#not sure why there are empty lists though!!!!!
			values_to_merge5 = [ x for x in overlapping_sets if len(x) > 0]

			#update sets for roi
			overlapping_sets = values_to_merge5	

		#after looping through all pairs once, use recursion to loop again
		#within a single roi
		#number of sets for each roi should decrease until it plateaus
		#return sets for each roi when it is no longer decreasing
		if prev_num_sets > len(overlapping_sets):
			return merge_overlapping_motifs2(overlapping_sets)
		elif prev_num_sets == len(overlapping_sets):
			print "end of recurisve round	", overlapping_sets
			print "max num of merges"
			return overlapping_sets 
		else:
			print "Error: the number  of sets should not be increasing"
			sys.exit()


############################################################################################
#file input
############################################################################################
#regions of interest
#bedFile = "hars_merged_hg19_chr1.bed"
#bedFile = "brain_expression_hars.bed"
bedFile = sys.argv[1]

#relevant transcription factors bindign
#motifFile = "hars_merged_hg19_chr1.match"
#motifFile = "brain_expression_hars.match"
motifFile = sys.argv[2] 

#options covariation matrix or motifs falling har variants
option = int(sys.argv[3])

#distances between motifs
#option == 1: 

#positions of variants along regons of interest
#har_mismatch_position_file = "brain_expression_hars_blast_to_panTro4_mismatch_v2.txt"
if option == 2: 
	har_mismatch_position_file = sys.argv[4]

#number of motif co-occurrences across regions of interest
if option == 3:
	co_occurrence_file = sys.argv[4]

#motifs present or absent across regions of interest
if option == 5:
	motif_present_across_roi = sys.argv[4]

#co-occurrence counts
if option == 6:
	co_occurrence_file_2 = sys.argv[4]

if option == 7:
	pickled_distances_file = sys.argv[4]

if option == 8: 
	pickled_distances_file = sys.argv[4]

#co-occurrence counts
if option == 10:
	co_occurrence_file_3 = sys.argv[4]

	
############################################################################################
#main
############################################################################################
#bed file containing regions of interest
bed = open(bedFile, "r")
regionsOfInterest=[]
for line in bed: 
	fields = line.rstrip("\n").split("\t")
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	name = fields[3]
	regionsOfInterest.append(region(name, chrom, start, end, "+"))

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
	
	

#distances between motifs
if option == 1:
	pickled_distances_bw_motifs_file = "co/cardiac/version4/distance_between_motifs.p"
	#pickled_distances_bw_motifs_file_header = "co/cardiac/version4/distance_between_motifs_header.p"

	#function call to calculate the distances between all motifs
	roi_cov_mat = distances_between_motifs_within_regions_of_interest(regionsOfInterest)
	for index, covariance_matrix in enumerate(roi_cov_mat):
		print regionsOfInterest[index].name, regionsOfInterest[index].length 
		row_col_names = [motifObject.tfname for motifObject in regionsOfInterest[index].motifs]
		print_matrix2(covariance_matrix, row_col_names)
	
	#pickle matrices
	if not os.path.isfile(pickled_distances_bw_motifs_file):
		pickle.dump(roi_cov_mat,    open(pickled_distances_bw_motifs_file, "w")   )
		#pickle.dump(row_col_names,  open(pickled_distances_bw_motifs_file_header, "w") )

#variants
if option == 2:
	har_mismatch_position_file_handle = open(har_mismatch_position_file, "r")
	for line in har_mismatch_position_file_handle: 
		fields = line.rstrip("\n").split("\t")
		mismatch_position = int(fields[0])
		regionName = fields[1]
		for index, regionObject in enumerate(regionsOfInterest): 
			if regionObject.name == regionName: 
				index_roi = index
		mismatchObject = mismatch(regionName, mismatch_position)
		regionsOfInterest[index_roi].add_mismatch(mismatchObject)
		#print mismatchObject.position_along_region

	locate_variants_within_motifs_within_hars(regionsOfInterest)

#matrix of motifs that occur in same enhancer
if option == 3: 
	#list of all tfnames present in all reigons of interest
	all_seen_tfnames = []
	for regionObject in regionsOfInterest:
		row_tfnames = [motifObject.tfname for motifObject in regionObject.motifs]
		all_seen_tfnames.extend(row_tfnames)

	unique_all_seen_tfnames = sort_and_deduplicate( all_seen_tfnames)
	
	#print out tf name in order (header)
	co_occurrence_file_tf = open( co_occurrence_file.replace(".txt", "_tfname_header.txt"), "w")
	co_occurrence_file_tf.write( "\t".join([ tfname_motif.replace("V$", "")  for tfname_motif in unique_all_seen_tfnames]) )
	co_occurrence_file_tf.write( "\n" )
	co_occurrence_file_tf.close()

	#matrix of all unique tf seen
	matrix  = np.zeros(shape=(len(unique_all_seen_tfnames),len(unique_all_seen_tfnames)) , dtype=np.int)

	#now count number of times two tf are seen together within a region of interest
	for regionObject in regionsOfInterest:
		row_tfnames = [str(motifObject.tfname) for motifObject in regionObject.motifs]
		unique_row_tfnames = sort_and_deduplicate(row_tfnames)
		for tfname1 in unique_row_tfnames: 
			position1 = [index for index,name in enumerate(unique_all_seen_tfnames) if str(name) == str(tfname1)]
			if len(position1) > 1: 
				print "this tf name should be unique"
				sys.exit()
			elif len(position1) == 0: 
				print "tf not found"
				sys.exit()
			for tfname2 in unique_row_tfnames: 
				position2 = [index for index,name in enumerate(unique_all_seen_tfnames) if str(name) == str(tfname2)]
				if len(position2) > 1: 
					print "this tf name should be unique"
					sys.exit()
				elif len(position2) == 0: 
					print "tf not found"
					sys.exit()
				#print position1, position2
				count = matrix[np.array(map(int, position1)), np.array(map(int, position2))]
				count = count + 1
				matrix[np.array(map(int, position1)), np.array(map(int, position2))] = count
			#print matrix

	#print matrix
	np.savetxt(co_occurrence_file, matrix)

if option == 4:
	test = [5,4,3,5,1,5,6,2,7,5,7,2,8,0]
	print test
	test2 = sort_and_deduplicate(test)
	print test2
	test3 = ['m', 'o', 'a', 'm', 'c', 'r', 'd', 'm', 'o']
	print test3
	test4 = sort_and_deduplicate(test3)
	print test4
	
# binary motif occurrence across roi
if option == 5:
	all_seen_tfnames = []
	for regionObject in regionsOfInterest:
		row_tfnames = [motifObject.tfname for motifObject in regionObject.motifs]
		all_seen_tfnames.extend(row_tfnames)

	unique_all_seen_tfnames = sort_and_deduplicate( all_seen_tfnames)

	#print out tf name in order (header)
	motif_occurrence_tf = open( motif_present_across_roi.replace(".txt", "_tfname_header.txt"), "w")
	motif_occurrence_tf.write( "\t".join([ tfname_motif.replace("V$", "")  for tfname_motif in unique_all_seen_tfnames]) )
	motif_occurrence_tf.write( "\n" )
	motif_occurrence_tf.write( "\t".join ( [ roi.name for roi in regionsOfInterest ]) )
	motif_occurrence_tf.write( "\n" )
	motif_occurrence_tf.close()

	#number of TF motifs (unique) by number of roi
	matrix  = np.zeros(shape=(len(unique_all_seen_tfnames),len(regionsOfInterest)) , dtype=np.int)
	
	for index1,regionObject in enumerate(regionsOfInterest): 
		row_tfnames = [str(motifObject.tfname) for motifObject in regionObject.motifs]
		unique_row_tfnames = sort_and_deduplicate(row_tfnames)
		for tfname1 in unique_row_tfnames: 
			index2 = [index for index,name in enumerate(unique_all_seen_tfnames) if str(name) == str(tfname1)]
			if len(index2) > 1: 
				print "this tf name should be unique"
				sys.exit()
			elif len(index2) == 0:
				print "tf not found"
				sys.exit()
			count = matrix[np.array(map(int, index2)), np.array( index1 )]
			count = count + 1
			matrix[np.array(map(int, index2)), np.array(index1)] = count

	np.savetxt(motif_present_across_roi, matrix)

#matrix of motifs that occur in same enhancer
if option == 6: 
	#list of all tfnames present in all reigons of interest
	all_seen_tfnames = []
	for regionObject in regionsOfInterest:
		row_tfnames = [motifObject.tfname for motifObject in regionObject.motifs]
		all_seen_tfnames.extend(row_tfnames)

	unique_all_seen_tfnames = sort_and_deduplicate( all_seen_tfnames)
	
	#print out tf name in order (header)
	co_occurrence_file_tf_2 = open( co_occurrence_file_2.replace(".txt", "_tfname_header.txt"), "w")
	co_occurrence_file_tf_2.write( "\t".join([ tfname_motif.replace("V$", "")  for tfname_motif in unique_all_seen_tfnames]) )
	co_occurrence_file_tf_2.write( "\n" )
	co_occurrence_file_tf_2.write( "\t".join ( [ roi.name for roi in regionsOfInterest ]) )
	co_occurrence_file_tf_2.write( "\n" )
	co_occurrence_file_tf_2.write( "\t".join ( map( str, [ roi.length for roi in regionsOfInterest ]) ) )
	co_occurrence_file_tf_2.write( "\n" )
	co_occurrence_file_tf_2.close(      )

	#matrix of all unique tf seen
	matrix  = np.zeros(shape=(len(unique_all_seen_tfnames),len(regionsOfInterest)) , dtype=np.int)

	#now count number of times two tf are seen together within a region of interest
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

	np.savetxt(co_occurrence_file_2, matrix)

if option == 7:
	#roi_cov_mat = pickle.load(open("co/cardiac/version4/distance_between_motifs.p", "rb"))
	
	#for index, covariance_matrix in enumerate(roi_cov_mat):
		#print regionsOfInterest[index].name, regionsOfInterest[index].length 
		#row_col_names = [motifObject.tfname for motifObject in regionsOfInterest[index].motifs]
		#print_matrix2(covariance_matrix, row_col_names)

	tf_pairs = distances_between_motifs_within_regions_of_interest_v2(regionsOfInterest)	

	#pickled_distances_file = "co/cardiac/version4/distance_between_motifs.p"
	if not os.path.isfile(pickled_distances_file):
		pickle.dump(tf_pairs, open(pickled_distances_file, "wb" ) )

	for key in tf_pairs.keys():
		med = np.median(tf_pairs[key])
		print key, "\t", med, "\t", "\t".join( map( str, tf_pairs[key] ) )

if option == 8:
	#pickled_distances_file = "co/cardiac/version4/distance_between_motifs.p"
	if os.path.isfile(pickled_distances_file):
		tf_pairs = pickle.load( open(pickled_distances_file, "rb" ) )

	for key in tf_pairs.keys():
		med = np.median(tf_pairs[key])
		print key, "\t", med, "\t", "\t".join( map( str, tf_pairs[key] ) )

#distances between all tf,  
if option == 9:
	#listOfTFs = ["GLIS2_04", "FOXO3_01", "FOXO4_01",  "MEF2A_Q6", "MEF2C_Q4", "ARID3A_04", "GATA6_Q5", "GATA6_01", "GATA6_04", "GATA4_Q5_01"]
	listOfTFs = ["SOX2_01", "SOX2_Q6", "POU5F1_01"]
	distances_between_motifs_within_regions_of_interest_v3(regionsOfInterest, listOfTFs)

#merge motifs that overlap
if option == 10:
	output  = open("overlapping_motifs.txt", "w")

	#find pairs of overlapping motifs
	overlapping_motifs = distances_between_motifs_within_regions_of_interest_v4(regionsOfInterest)
	
	clusters_motifs = [] #list of clusters of overlapping motifs across roi
	clusters_and_singlets_by_roi = {}
	singlets_motifs2 = [] #list of singlet motifs across roi 

	for roi_index in overlapping_motifs.keys():
		regionObject = regionsOfInterest[roi_index] 
		
		clusters_index = [] #list of motif indices  

		print "++++++++++++++++++++++++++++++"
		print "roi_index", roi_index
	
		#store motif names for each roi
		clusters_and_singlets_by_roi[roi_index] = []

		#merge pairs of overlapping motifs
		roi_overlapping_motifs = merge_overlapping_motifs2(overlapping_motifs[roi_index])
		for element in roi_overlapping_motifs:

			merge_motifs = [ str(motifObject.tfname) for index, motifObject in enumerate(regionObject.motifs) if index in element ]
			merge_motifs_positions  = [ motifObject.position for index, motifObject in enumerate(regionObject.motifs) if index in element ]
			merge_motifs_strands  = [ motifObject.strand for index, motifObject in enumerate(regionObject.motifs) if index in element ]
			merge_motifs_positions2  = [ regionObject.return_motif_start_position(motifObject) for index, motifObject in enumerate(regionObject.motifs) if index in element ]
			merge_motifs_positions3  = [ regionObject.return_motif_end_position(motifObject) for index, motifObject in enumerate(regionObject.motifs) if index in element ]
			#convert clusters into string motif names
			merge_motifs2 = sort_and_deduplicate(merge_motifs)
			merge_motifs3 = map(str, merge_motifs2)
			merge_motifs4 = "|".join([ tfname_motif.replace("V$", "")  for tfname_motif in merge_motifs3 ])

			clusters_index.extend(element) #list of cluster motifs indices found for single roi 
			clusters_motifs.append(merge_motifs4) #list of clusters (also lists) of overlapping motifs across roi
			clusters_and_singlets_by_roi[roi_index].append(merge_motifs4) #retrieve motif cluster names for single roi

			output.write("%s\t%s\n" %( roi_index, "\t".join( map(str, element ) ) ) )  
			output.write( "mm\t%s\n" % ( "\t".join( map(str, merge_motifs) ) ) )
			if len(merge_motifs) == 1: output.write("hi\n")
			if len(element) == 1: output.write("hi2\n")
			output.write( "\t%s\n" % ( "\t".join( map(str, merge_motifs_positions) ) ) )
			output.write( "\t%s\n" % ( "\t".join( map(str, merge_motifs_strands) ) ) )
			output.write( "\t%s\n" % ( "\t".join( map(str, merge_motifs_positions2) ) ) )
			output.write( "\t%s\n" % ( "\t".join( map(str, merge_motifs_positions3) ) ) )
			output.write( "\n" )
		
		#list of motifs not overlapping with other motifs
		roi_indices = set(range( 0, (len(regionObject.motifs) + 1 )))
		singlets_index = list(roi_indices.difference(clusters_index))
		singlets_motifs = [ str(motifObject.tfname) for index, motifObject in enumerate(regionObject.motifs) if index in singlets_index ]
		clusters_and_singlets_by_roi[roi_index].extend(singlets_motifs) #retrieve motif singlet names for single roi
		singlets_motifs2.extend( singlets_motifs ) #store motif singlets across roi

		output.write( "single\t%s\n" % ( "\t".join( map(str, singlets_motifs) ) ) )


	#set of all motifs present (clusters of motifs + single motifs)
	clusters_motifs2 = sort_and_deduplicate(clusters_motifs)
	singlets_motifs3 = sort_and_deduplicate(singlets_motifs2)
	output.write( "singlea2\t%s\n" % ( "\t".join( map(str, singlets_motifs2) ) ) )
	output.write( "singlea3\t%s\n" % ( "\t".join( map(str, sorted(singlets_motifs2, reverse=True)) ) ) )
	output.write( "singlea4\t%s\n" % ( "\t".join( map(str, singlets_motifs3) ) ) )
	output.write( "cls1\t%s\n" % ( "\t".join( map(str, clusters_motifs) ) ) )
	output.write( "cls2\t%s\n" % ( "\t".join( map(str, sorted(clusters_motifs, reverse=True)) ) ) )
	output.write( "cls3\t%s\n" % ( "\t".join( map(str, clusters_motifs2) ) ) )
	output.close()
	all_sets = clusters_motifs2 + singlets_motifs3
	all_sets2 = sort_and_deduplicate(all_sets)

	#print out tf name in order (header)
	co_occurrence_file_tf_3 = open( co_occurrence_file_3.replace(".txt", "_tfname_header.txt"), "w")
	co_occurrence_file_tf_3.write( "\t".join([ tfname_motif.replace("V$", "")  for tfname_motif in all_sets2]) )
	co_occurrence_file_tf_3.write( "\n" )
	co_occurrence_file_tf_3.write( "\t".join ( [ roi.name for roi in regionsOfInterest ]) )
	co_occurrence_file_tf_3.write( "\n" )
	co_occurrence_file_tf_3.write( "\t".join ( map( str, [ roi.length for roi in regionsOfInterest ]) ) )
	co_occurrence_file_tf_3.write( "\n" )

	#matrix of all unique tf seen
	matrix  = np.zeros(shape=(len(all_sets2),len(regionsOfInterest)) , dtype=np.int)

	#now count number of times two tf are seen together within a region of interest
	for index1,regionObject in enumerate(regionsOfInterest): 
		if index1 > 18: break
		#clusters and singlets for each roi
		#row_tfnames = [str(motifObject.tfname) for motifObject in regionObject.motifs]
		row_tfnames = clusters_and_singlets_by_roi[index1]
		for tfname1 in row_tfnames: 
			index2 = [index for index,name in enumerate(all_sets2) if str(name) == str(tfname1)]
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

	np.savetxt(co_occurrence_file_3, matrix)

