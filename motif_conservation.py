########################################################################################
#  Calculate the average motif conservation for each instance of a motif               #
########################################################################################
import sys
sys.path.append('/home/tara/bin/mine')
from network_classes import *
from maf_functions import * 
import pickle
import os.path
import numpy as np

#option = int(sys.argv[1])

#files
#motifFile = "co/cardiac/version4/table_s4_union_enhancer_ESC_chr19.match"
##wigFile = "phyloP46way_vertebrate_basewise/chr10.phyloP46way.wigFix"
#wigFile = "chr19.phyloP60way.wigFix"
#bedFile = "co/cardiac/table_s4_union_enhancer_ESC_chr19.bed"
#pickled_phyloP_scores_file = "phyloP46way_vertebrate_basewise_intersected_with_regions_of_interest_chr19_v2.p"

motifFile = sys.argv[1]
wigFile = sys.argv[2]
bedFile = sys.argv[3]
pickled_phyloP_scores_file = sys.argv[4]


if os.path.isfile(pickled_phyloP_scores_file):
	option = 0
else: 
	option = 1

#intersect motif genomic coordinates with store_all_scores[0] coordinates
#sorted bed file
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

#retrieve basewise conservation scores 
#option == 0 : upload pickled store_all_scores
#option == 1 : create new store_all_scores
if option == 1:
	wig = open(wigFile, "r")
	prevChrom = ""
	prevStart = 1
	currPos = 0 
	store_all_scores = {}
	for line in wig: 
		fields = line.rstrip("\n").split()
		if line.startswith("fixed"):
			#print scores for regions of interest that overlap with this wig interval 	
			for regionObject in regionsOfInterest: 
				name = regionObject.name
				roi = ( regionObject.chrom, regionObject.start, regionObject.end, "+" )
				overlap = intersectF ( (prevChrom, prevStart, currPos), roi )
				if overlap != 0: 
					#print overlap, prevStart, currPos, roi
					all_scores = store_scores [ ( overlap[1] - prevStart ) : ( overlap[2] - prevStart ) ] 
					if store_all_scores.has_key(name):
						store_all_scores[name][0].append(overlap) 
						store_all_scores[name][1].append(store_scores)
					else: 
						store_all_scores[name] = ([overlap], [store_scores])
		
			#print line.rstrip("\n")
			store_scores= []
			currChrom = fields[1].split("=")[1]
			start = int(fields[2].split("=")[1])
			currPos = start
			prevChrom = currChrom
			prevStart = start
		else: 
			currPos = currPos + 1
			score = float(fields[0])
			store_scores.append(score)
		
	pickle.dump(store_all_scores, open(pickled_phyloP_scores_file, "wb" ) )
	wig.close()

#basewise conservation score by region of interest
if option == 3: 
	for regionObject in regionsOfInterest: 
		name = regionObject.name
		print "###########################"
		print name
		print regionObject.chrom, regionObject.start, regionObject.end
		if store_all_scores.has_key(name): 
			for index, coords in enumerate(store_all_scores[name][0]):
				print coords
				print store_all_scores[name][1][index]
		else: 
			print "region did not overlap with conservation wig"


#basewise conservation score for each instance of a motif within a region of interest
#loop through each motif instance within each region of interest
store_all_scores = pickle.load(  open(pickled_phyloP_scores_file, "rb") )
print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("region_ID", "region_chr", "region_start", "region_end", "score_block_chr", "score_block_start", "score_block_end", "overlap_chr", "overlap_start", "overlap_end", "motif_ID", "motif_chr", "motif_start", "motif_end", "mean_motif_conservation", "sd_motif_conservation")
for regionObject in regionsOfInterest:
	name = regionObject.name
	if store_all_scores.has_key(name): 
		regionObject_scores = store_all_scores[name]
		motifChrom = regionObject.chrom
		for motifObject in regionObject.motifs:
			motifStart = regionObject.return_motif_genomic_start_position(motifObject)
			motifEnd = regionObject.return_motif_genomic_end_position(motifObject)
			#loop through conservation wig file coordinates
			for index, regionObject_scores_coords in enumerate(regionObject_scores[0]):
				#intersect motif genomic coordinates with store_all_scores[0] coordinates
				overlap = intersectF ( (motifChrom, motifStart, motifEnd), regionObject_scores_coords ) 
				if overlap != 0:
					#pull out store_all_scores[1] that correspond to that intersection
					motif_scores = np.array(regionObject_scores[1][index][overlap[1] - regionObject_scores_coords[1]: overlap[2] - regionObject_scores_coords[1] ] )
					#take the average
					motif_average_score = np.mean(motif_scores, dtype=np.float64)
					motif_sd = np.std(motif_scores, dtype=np.float64)
					#print name, "|roi",  regionObject.start, regionObject.end, "|o", overlap[0], overlap[1], overlap[2], "|", motifObject.tfname,  motifChrom, motifStart, motifEnd, "|score_coords",  regionObject_scores_coords,  motif_average_score
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (name, regionObject.chrom,  regionObject.start, regionObject.end, regionObject_scores_coords[0], regionObject_scores_coords[1], regionObject_scores_coords[2], overlap[0], overlap[1], overlap[2],  motifObject.tfname,  motifChrom, motifStart, motifEnd, motif_average_score, motif_sd)
					
#issue: how do I deal with missing conservation scores

