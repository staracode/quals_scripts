import sys
#blastFile = "brain_expression_hars_blast_to_panTro4.blast"
blastFile = sys.argv[1] 

blastHandle = open(blastFile, "r")

keepNextTwoLines=-1
lastQueryStart=-200
global queryName
for line in blastHandle: 
	fields = line.rstrip("\n").split()
	if line.startswith ("Query="): 
		queryName = fields[1]
		#print queryName
		
	elif line.startswith("Query:"): 
		keepNextTwoLines = 2
		queryStart = int(fields[1])
		querySeq = fields[2]
		queryEnd = int(fields[3])
		#use this to find position of mismatch along whole alignment
		lastQueryStart = queryStart
		#print querySeq, queryStart, queryEnd, len(querySeq)
	else: 
		if keepNextTwoLines == 2: 
			sbjctSeq = line.lstrip("                ").rstrip("\n")
			#print sbjctSeq
			#mismatch_position = [i for (i,x) in enumerate(sbjctSeq) if x == " "]
			mismatch_position2 = [i + lastQueryStart for (i,x) in enumerate(sbjctSeq) if x == " "]
			#print mismatch_position, mismatch_position2
			#print mismatch_position2
			if mismatch_position2:
				print "".join([str(mismatch) + "\t" + queryName + "\n"  for mismatch  in mismatch_position2]).rstrip("\n")
		elif keepNextTwoLines == 1: 
			subjStart = int(fields[1])
			subjSeq = fields[2]
			subjEnd = int(fields[3])
			#print subjSeq, subjStart, subjEnd
			
		keepNextTwoLines = keepNextTwoLines -1 

		#else don't print the lines
	
