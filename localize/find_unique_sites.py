def unique_pos(position_directory, unique_length, prefix):
	import os
	insert_sites = {}

	for filename in os.listdir(position_directory):
		for line in open(os.path.join(position_directory, filename)):
			if 'nsert' in line:
				continue
			
			content = line.strip().split(',')
			gene = content[0]
			ins_chr = content[4]
			ins_pos = content[5]
			parent_chr = content[1]
			parent_start = content[2]
			parent_end = content[3]

			if parent_chr == ins_chr and parent_start <= ins_pos <= parent_end:
				continue
			if gene not in insert_sites:
				insert_sites[gene] = [[ins_chr, ins_pos, 1, [ind], parent_chr, parent_start, parent_end]]
			elif gene in insert_sites:
				counter = 0
				for s in range(len(insert_sites[gene])):
					temp_chr = insert_sites[gene][s][0]
					temp_pos = insert_sites[gene][s][1]
					if ins_chr == temp_chr and abs(int(temp_pos)-int(ins_pos)) <= unique_length:
						temp += 1
						insert_sites[gene][s][2] += 1
						insert_sites[gene][s][3] += [ind]
				if counter == 0:
					insert_sites[gene] += [[ins_chr, ins_pos, 1, [ind], parent_chr, parent_start, parent_end]]		

	with open('{}_unique_inserts.txt'.format(prefix), 'w') as file:
		file.writelines('gene,chromosome,position,count,parent_chr,parent_start,parent_end)
		for gene in insert_sites:
			for pos in insert_sites[gene]:
				file.writelines('{},{},{},{},{},{},{}\n'.format(gene, post[0], pos[1], str(pos[2]), pos[4], pos[5], pos[6]))
	return insert_sites

