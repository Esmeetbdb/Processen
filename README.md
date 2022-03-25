# Processen
Pipeline for detection of non-reference processed pseudogenes

# Discover
Pseudogenes are discovered through the discovery workflow. 
The discovery workflow uses samtools to extract read pairs bridging distant exons, and aligns these pairs to the reference using Salmon quant.
The steps are detailed in discover.sh

The workflow requires bam file as input, as well as a bed file describing the position of the exons of protein coding genes.
The bed file needs to be in the following format:

	chromosome	start	end	exon_id

the bed file needs be indexed using tabix.

The discovery pipeline is run through the folloing command:

	bash discover.sh input.bam exons.bed.gz prefix salmon_index

# Position


# Dependencies

	samtools
	salmon
	tabix
