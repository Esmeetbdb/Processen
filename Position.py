def Unique_Pos(dir_name, unique_length, output_file):
    '''This function checks whether pseudogene insertion sites between different individuals with the same pseudogene are unique.

    Input:
    dir_name: name of the directory where files with individuals with pseudogenes are stored. Files in this directory should be comma
    separated
    unique_length: the minimum number of nucleotides that should be between two insertion points to be considered unique.
    output_file: name of the desired output file

    Output:
    file (name given by user): structure is insert chromosome, insert position, name of pseudogene,
    number of individuals with this insert point, first 3 individuals with that insert site.
    '''

    import os
    insert_sites = {}

    for filename in os.listdir(dir_name):
        with open(os.path.join(dir_name, filename), 'r') as file:
            for line in file:
                if 'nsert' in line:
                    continue
                else:
                    content = line.strip().split(',')
                    ind = str(filename)+'_'
                    gene = content[0]
                    ins_chr = content[4]
                    ins_pos = content[5]
                    og_chr = content[1]
                    og_start = content[2]
                    og_end = content[3]

                    if og_chr == ins_chr and og_start <= ins_pos <= og_end:
                        continue
                    elif gene not in insert_sites:
                        insert_sites[gene] = [[ins_chr,ins_pos,1, ind, og_chr, og_start, og_end]]
                    elif gene in insert_sites:
                        temp = 0

                        for s in range(len(insert_sites[gene])):
                            temp_chr = insert_sites[gene][s][0]
                            temp_pos = insert_sites[gene][s][1]
 if ins_chr == temp_chr and abs(int(temp_pos)-int(ins_pos)) <= unique_length:
                                temp += 1
                                insert_sites[gene][s][2] += 1
                                if insert_sites[gene][s][2] <4:
                                    insert_sites[gene][s][3] += ind
                        if temp == 0:
                            insert_sites[gene] += [[ins_chr, ins_pos, 1, ind, og_chr, og_start, og_end]]

    with open(output_file, 'w') as nfile:
        nfile.write('chromosome,position,gene_name,count,individuals\n')
        for gen in insert_sites:
            for positions in insert_sites[gen]:
                chr = positions[0]
                pos = positions[1]
                count = str(positions[2])
                inds = positions[3]
                og_chr = positions[4]
                og_start = positions[5]
                og_end = positions[6]
                nfile.write(gen+','+chr+','+pos+','+count+','+inds+','+og_chr+','+og_start+','+og_end+'\n')
