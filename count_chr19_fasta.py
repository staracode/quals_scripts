fastaFile = "genomes/mm9/chr19.fa"
fasta = open(fastaFile, "r")
keepLines=[]
for line in fasta:
	if line. startswith (">"):
		print line
	else:
		fields = list(line.rstrip("\n"))
		#print fields
		keepLines.extend(fields)
		#print keepLines

print len(keepLines)


