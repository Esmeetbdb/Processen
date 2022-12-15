def create_pseudogene_dict(summary_file):
	pseudogene_dict = {}
	for line in open(summary_file);
		if 'vidual' in line:
			continue
		content = line.strip().split(',')
		ind = contetn[0]
		gene = content[2]

		if ind not in pseudogene_dict and gene != 'none':
			pseudogene_dict[ind] = [gene]
		elif ind in pseudogene_dict and gene != 'none':
			pseudogene_dict[ind].append(gene)
	
	return pseudogene_dict

 
def create_gene_info_dict(gene_position_info):
	gene_info_dict = {}

	for line in open(gene_position_info):
		if 'name' in line:
			continue
		content = line.strip().split(',')
		gene = content[0]
		start = content[1]
		end = content[2]
		chr = content[3]

		gene_info_dict[gene] = (chr, start, end)
	return gene_info_dict

def pseudogene_localization(path_to_VCF, pseudogene_dict, gene_info_dict, output_folder, prefix, r):
	import readVCF 
	import sys
	for path in open(path_to_VCF):
		content = path.strip().split('\')
		ind = content[1] # Change this so that it is more universal
		output = []
		gene_list_ind = pseudogene_dict[ind]
		
		for gene in gene_list:
			chr = gene_info_dict[gene][0]
			start =  gene_info_dict[gene][1]
			end = gene_info_dict[gene][2]
			
			for line in open(path.strip()):
				if line[0] == '#':
					continue
				try:
					variant = readVCF.readVCFLine(line)
				except:
					continue
				# Figure out how this works again to annotate code and add more checks probably
				if variant[0] == variant[2]:
					if abs(int(variant[1]) - int(variant[3])) < abs(start - end):	
						continue
				if variant[0] == chr:
					if abs(start - int(variant[1])) < r:
						new_chr = variant[2]
						insert_pos = variant[3]
						output.append([gene, chr, start, end, new_chr, insert_pos])
						continue
					elif abs(end - int(variant[1])) < r:
						new_chr = variant[2]
						insert_pos = variant[3]
						output.append([gene, chr, start, end, new_chr, insert_pos])
						continue
				
				if variant[2] == chr:
					if abs(start - int(variant[3])) < r:
						new_chr = variant[0]
						insert_pos = variant[1]
						output.append([gene, chr, start, end, new_chr, insert_pos])
						continue
					elif abs(end - int(variant[3])) < r:
						new_chr = variant[0]
						insert_pos = variant[1]
						output.append([gene, chr, start, end, new_chr, insert_pos])
						continue

		if len(output) >= 1:
			file_name = output_folder + '/{}_{}.csv'.format(prefix,ind)
			with open(file_name, 'w') as new_file:
				new_file.writelines('Gene,Chr,Start,End,Insert_Chr,Insert_Position\n')
				for loc in output:
					new_file.writelines('{}, {}, {}, {}, {}, {}\n'.format(loc[0],str(loc[1],str(loc[2]),str(loc[3]),str(loc[4]),str(loc[5])))


