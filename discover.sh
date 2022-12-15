#usage:
#bash discover.sh input.bam exons.bed.gz prefix salmon_index

#extract exonic reads having mate aligning on same chromosome
samtools view -h $1 -L $2 -F 2 | grep -E "@|	=	" | samtools view -Sh -q 20 -F 4 -F 8 - > $3.sam
#extract pairs 
python discover/find_spanning_pairs.py $3.sam $2 $3 | samtools view -Shb - > $3.filt.bam
#convert to fastq
samtools sort -n $3.filt.bam | samtools fastq - -F 2048 -F 256 -s $3_single.fq -1 $3_R1.fq -2 $3_R2.fq
#quantify transcripts using salmon
salmon quant -i $4 -l A -1 $3_R1.fq  -2 $3_R2.fq -o $3 --validateMapping
python filter_quant.py $3/quant.sf | grep -v "|MT-" > $3.pseudogenes.sf
