import sys
import os
#python find_spanning_pairs.py input.sam exons.bed outprefix

for line in open(sys.argv[1]):
	if line[0] == "@":
		print(line.rstrip())
		continue

	content=line.rstrip().split()
	start=int(content[3])
	end=int(content[7])

	if start == end:
		continue
	#print(line.strip())
	exons_start=set([])
	exons_end=set([])

	os.system( "tabix {} {}:{}-{} > {}.tmp.bed".format(sys.argv[2],content[2],start,start+150,sys.argv[3]) )
	for e in open("{}.tmp.bed".format(sys.argv[2])):
		c=e.rstrip().split("\t")
		if len(c) > 2:
			exons_start.add(c[-1])

	if not exons_start:
		continue

	os.system( "tabix {} {}:{}-{} > {}.tmp.bed".format(sys.argv[2],content[2],end,end+150,sys.argv[3]) )
	for e in open("{}.tmp.bed".format(sys.argv[2])):
		c=e.rstrip().split("\t")
		if len(c) > 2:
			exons_end.add(c[-1])

	if not exons_end:
		continue

	same_exon=False
	for exon in exons_end:
		if exon in exons_start:
			same_exon=True
	
	if same_exon:
		continue

	print(line.rstrip())
