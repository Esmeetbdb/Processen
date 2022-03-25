import sys

read_length=150
min_cov=3

first=True
for line in open(sys.argv[1]):
	if first:
		print line.strip()
		first=False
		continue
	
	content=line.strip().split()
	if read_length*float(content[-1])/float(content[2]) > min_cov and float(content[-1]) > 20:
		print(line.strip())
