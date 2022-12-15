def insert_generation(n, output_file, out_path, prefix):
	import random

	chr_list = ['hg1', 'hg2', 'hg3', 'hg4', 'hg5', 'hg6', 'hg7', 'hg8', 'hg9', 'hg10', 'hg11', 'hg12', 'hg13', 'hg14',
                'hg15','hg16', 'hg17', 'hg18', 'hg19', 'hg20', 'hg21', 'hg22', 'hgX', 'hgY']
	chr_len = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616,
                 64444167, 46709983, 50818468, 156040895, 57227415]  # Chromosome lengths in order hg1-22, hgX, hgY

	rand_chr = []
	temp = []
	total_length = [0, 248956422]
	MC_inserts = []
	for i in range(1, len(chrom_len)):
		total_length.append(sum(chrom_len[0:i + 1]))

	for i in range(n):
		temp.append(random.randint(1, total_length[-1]))

	for item in temp:
		for i in range(0,len(chr_list)):
			if item >= total_length[i] and item < total_length[i+1]:
				rand_chr.append(chr_list[i])
				break
			else:
				continue
	for chr in rand_chr:
		index = chr_list.index(chrom)
		start = str(random.randint(1, chr_len[index]))
		end = start +1
		MC_inserts.append((chr,start,end))

	if output_file != False:
		with open(out_path+prefix+'MC_insertlocations.txt', 'w') as file:
			for chr in rand_chr:
				index = chr_list.index(chrom)
				file.writelines('{},{}\n'.format(chr,str(random.randint(1, chr_len[index]))))
	return MC_inserts

def loc_file_to_list(insert_location_file):
	inserts = []
	for line in open(insert_location_file):
		if 'ene' in line or if 'art' in line:
			continue
		content = line.strip().split(',')
		chr = 'chr' + content[1].replace('hg', '').replace('chr', '').replace('CHR', '')
		start = int(content[2])
		end = start += 1
		inserts.append((chr,start,end))
	return inserts

def count_chrom(insert_list):
	count = 0
	chr_list = ['hg1', 'hg2', 'hg3', 'hg4', 'hg5', 'hg6', 'hg7', 'hg8', 'hg9', 'hg10', 'hg11', 'hg12', 'hg13', 'hg14',
                'hg15', 'hg16', 'hg17', 'hg18', 'hg19', 'hg20', 'hg21', 'hg22', 'hgX', 'hgY']

	chr_count = {}
	for loc in insert_list:
		chr = loc[0]
		if chr not in chr_list:
			continue
		count+=1
		if chr not in chr_count:
			chr_count[chr] = 1
		else:
			chr_count += 1

	for chr in chr_list:
		if chr not in chr_count:
			chr_count[chr] = 0

	return chr_count, count


def chromosome_binomial(MC_insert_list, insert_list, prefix, out_path):
	from scipy.stats import binom_test
	import math
	MC_chr_count, MC_count = count_chrom(MC_insert_list)
	Actual_chr_count, Actual_count = count_chrom(insert_list)


	actual_chr_count_proportion = {}
	for chr in Actual_chr_count:
		actual_chr_count_proportion[chr] = Actual_chr_count[chr]/Actual_count

	MC_chr_count_proportion = {}
	for chr in MC_chr_count:
		MC_chr_count_proportion[chr] = MC_chr_count[chr]/MC_count

	with open(out_path+prefix+'Pseudogene_chr_binomial.csv', 'w') as file:
		file.writelines('Chromosome, log value, subtracted value, p-value (not corrected for multiple testing\n')
		for chr in chr_count_loc:
			k = chr_count_loc[chr]
			p = MC_chr_count_proportion[chr]
			sig = binom_test(k, count_loc, p)
			
			val_div = actual_chr_count_proportion[chr]/MC_chr_count_proportion[chr]
			if val_div != 0:
				log_val = math.log(val_div, 2)
			else:
				log_val = 'N.A.'
			val_sub = actual_chr_count_proportion[chr]-MC_chr_count_proportion[chr]
			file.writelines('{},{},{},{}\n'.format(chr,str(log_val),str(val_sub),str(sig))

# Change so that making output file is optional
# Instead of specifying full output file name make to specify a prefix for naming
def check_gene_track(track_loc_file, insert_list, out_file):
	import os
	gen_track = []
	for loc in insert_list:
		chr = loc[0]
		start = loc[1]
		end = loc[2]

		query = '{}:{}-{}'.format(chr,start,end)
		os.system('tabix {} {} > {}'.format(track_loc_file, query, 'tmp.txt'))
	
		for info in open('tmp.txt'):
			gen_track.append(info.strip())



	if output_file != False:
		with open(out_file, 'w') as file:
			for item in gen_track:
				file.writelines('{}\n'.format(item))

	return gen_track


def bionomial_dist(n, k, x, y):
	from scipy.stats import binom_test

	p=y/x
	sig = binom_test(k,n,p)
	value_sub = (k/n)-p
	value_div = (k/n)/p
	return sig, value_sub, value_div

def check_all_track(path_to_tabix_files, Actual_insert_list, MC_insert_list, out_file, out_path, prefix):
	outfile = open(out_path+'pseudogene_track_overrepresentation.csv', 'w')
	outfile.writelines('Trackname, p-value, value_subtracted, value_divided\n')
	for file in open(path_to_tabix_files):
		path = file.strip()
		content = file.strip().split('/')
		track = content[-1]
		if out_file != False:
			output = out_path + '{}.{}_locations'.format(prefix,track)
		else:
			output = False

		if len(Actual_insert_list) ==0:
			outfile.writelines('{}, NA, NA, NA\n'.format(track))		
			continue
		track_actual = check_gene_track(path, Actual_insert_list, output + '_actual.txt')
		track_MC = check_gene_track(path, MC_insert_list, output + '_MC.txt')

		sig, value_sub, value_div = binomial_dist(len(Actual_insert_list), len(track_actual), len(MC_insert_list), len(track_MC))
		outfile.writelines('{}, {}, {}, {}\n'.format(track, sig, value_sub, value_div)
