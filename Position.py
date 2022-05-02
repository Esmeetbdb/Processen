def ind_gene_dict(Summary_file):
    '''Generate a dictionary that has the individual as a key and a lit of all genes found as pseudogenes in this ind
    as the value.
    Input:
    Summary_file: csv file that contains information about individuals and their pseudogenes.
    Output:
    Dictionary that contains information about individuals and what pseudogenes they have.'''

    dict = {}
    with open(Summary_file) as file:
       for line in file:
            if 'vidual' in line:
               continue

            content = line.strip().split(',')
            ind = content[0]
            gene = content[2]

            if ind not in dict and gene != 'none':
                dict[ind] = [gene]
            elif ind in dict and gene != 'none':
                dict[ind] += [gene]
    return dict

def biomart_dict(biomart_file):
    '''Turns a file containing genome information into a dictionary. The Biomart file must contain gene names and positional information
    of each gene.'''
    Biomart_dict = {}
    with open(biomart_file) as file:
        for line in file:
            if 'name' in line:
                continue
            content = line.strip().split(',')
            gene = content[0]
            start = content[1]
            end = content[2]
            chr = content[3]
            Biomart_dict[gene] = [chr, start, end]
    return Biomart_dict

def pseudogene_localization(VCFpaths, Summary_file, biomart_file, output_folder, r = 1000):
    '''Finds the insert position of pseudogenes using VCF files. The function checks if one of the VCF variant call positions is near the original gene position,
    which is found using the biomart file. If the other position is very far away it is assumed that this is the insert location of the pseudogene.
    Input:
    VCFpaths: file that contains paths to all VCF files that should be checked
    r: the max distance a point can be from the o.g. location of the gene to be considered further. If this number is made smaller less variant calls will be attributed to this gene.
    If it is made too big, variants that are outside of this gene may also be considered for the location of the pseudogene.
    Summary_file: file that contains information about which individuals have which pseudogenes.
    biomart_file: file with information about original gene locations
    output folder: the location where you want all files with pseudogene location information to be stored.

    Output:
    a file for each individual that the pseudogenes could be localized off. The file has one line for each pseudogene insert location, with information about which gene,
    it's original location and insert location'''
    import readVCF
    import sys
    Ind_gene_dict = ind_gene_dict(Summary_file)
    Biomart_Dict = biomart_dict(biomart_file)
    for path in open(VCFpaths):
        content = path.strip().split('/')
        ind = content[1]
        output = []
        gene_list = Ind_gene_dict[ind]
#        print(gene_list)
        for gene in gene_list:
#            print(gene)
            chr = Biomart_Dict[gene][0]
            start = int(Biomart_Dict[gene][1])
            end = int(Biomart_Dict[gene][2])

            with open(path.strip()) as VCF_file:
                for line in VCF_file:

                    if line[0] == '#':
                        continue

                    try:
                        variant = readVCF.readVCFLine(line)
                    except:
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
            file_name = output_folder + '/{}.csv'.format(ind)
            with open(file_name, 'w') as new_file:
                new_file.writelines('Gene,Chr,Start,End,Insert_Chr,Insert_Position\n')
                for l in output:
                    a = l[0]
                    b = str(l[1])
                    c = str(l[2])
                    d = str(l[3])
                    e = str(l[4])
                    f = str(l[5])
                    new_file.writelines(a + ',' + b + ',' + c + ',' + d + ',' + e + ',' + f + '\n')

                            

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
