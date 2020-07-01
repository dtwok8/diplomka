import pysam
import sys
import os
import pickle

import config

import subprocess
""" Because there can be more files, like nuc and another"""
END_REFERENCE_GEN_FILE = "_gen.fasta"




"""
	Make statistics like coverage of allels and sum of histogram covarage.
"""
def count_statistic(sort_bam_file: str):
	coverage_string = pysam.depth("-aa", sort_bam_file)
	coverage_list= coverage_string.splitlines()

	curent_alel =""
	# number of positions covered 
	number_nuc_coverage = 0
	sum_nuc_coverage = 0
	alel_size = 0
	max_histogram = 0
	alels_statistics = dict()
	cov_nucleotids = list()

	for item_nuc in coverage_list:
		item_nuc_list = item_nuc.split()

		if(item_nuc_list[0] != curent_alel):
			#skip the first change
			if(curent_alel != ""):
				alels_statistics[curent_alel] = {}
				alels_statistics[curent_alel]['number_nuc_coverage'] = number_nuc_coverage
				alels_statistics[curent_alel]['sum_nuc_coverage'] = sum_nuc_coverage
				alels_statistics[curent_alel]['alel_size'] = alel_size
				alels_statistics[curent_alel]['max_histogram'] = max_histogram
				alels_statistics[curent_alel]['cov_nucleotids'] = cov_nucleotids


			curent_alel = item_nuc_list[0]
			number_nuc_coverage = 0	
			sum_nuc_coverage = 0
			alel_size = 0
			max_histogram = 0
			cov_nucleotids = list()

		alel_size += 1
		cov_nucleotids.append(int(item_nuc_list[2]));

		if(int(item_nuc_list[2]) != 0):
			number_nuc_coverage += 1
			sum_nuc_coverage += int(item_nuc_list[2])

			if(max_histogram < int(item_nuc_list[2])):
				max_histogram = int(item_nuc_list[2])

	return alels_statistics


"""
	Evaluate statistics, find max coverage of each gene
	and write it into result file.
"""
def evaluate_statistics(haplotype_gen_statistics: dict, output_file_name: str):

	for alel, values in haplotype_gen_statistics.items():
		values["normalize_coverage"] = values['number_nuc_coverage'] / (values['alel_size'] / 100)
			
	# asi muzeme predpokladat za budou serazeny
	# takze by to asi jen chtelo seradit
	# {'KIR:KIR00001': {'number_nuc_coverage': 7407, 'sum_nuc_coverage': 19500, 'alel_size': 14738, 'max_histogram': 6, 'normalize_coverage': 50.25783688424481}, 'KIR:KIR00978': 
	#sorted_x = sorted(haplotype_gen_statistics.items(), key=operator.itemgetter(1["normalize_coverage"]))
	
	sorted_x = {k: v for k, v in sorted(haplotype_gen_statistics.items(), key=lambda item: item[1]["normalize_coverage"], reverse= True)}
	# takze nejdriv je seradit podle coverage
	output_file = open(output_file_name, "w")

	for alel, values in sorted_x.items():
		#if(values['number_nuc_coverage'] > 0):
		print(alel, " c: ", values['normalize_coverage'], " sum: ", values['sum_nuc_coverage'], "max histogram: ", values['max_histogram'] , file=output_file)

		#print(alel, " c: ", values['normalize_coverage'], "max_value_in_histogram", values['max_histogram'], file=output_file)

	output_file.close()
	print("create result file ", output_file_name)
	return sorted_x


""""
	Find all haplotypes in sam folder.
"""
def find_haplotypes(aligments_sam_files: list):
	gene_statistics = dict()
	haplotype_list = list()
	for file in aligments_sam_files:
		file_name_split = file.split('_')
		# find haplotype name 
		haplotype_name = ""	
		for item in file_name_split:
			if(item.startswith("KIR")):
				break
			else:
				haplotype_name += item+"_"
		haplotype_list.append(haplotype_name)

	# find unique
	haplotype_name_set = set(haplotype_list)
	haplotype_name_list_unique = list(haplotype_name_set)
	return haplotype_name_list_unique


def prepare_sam_bam_file(sam_file: str, haplotype: str,file_sufix: str):
	file_name_split = sam_file.split('_')
	
	# find KIR 
	kir_name = ""
	
	for item in file_name_split:
		if(item.startswith("KIR")):
			kir_name= item
			break

	print("preparing file for ", haplotype, kir_name, "...")
	# create bam file and  conver sam to bam and then sort 
	bam_file = os.path.join(config.BAM_FOLDER, haplotype+kir_name+".bam")
	fh = open(bam_file, 'w')
	fh.close()

	pysam.view('-o', bam_file, '-b', os.path.join(config.ALIGNMENT_FOLDER, sam_file), save_stdout=bam_file)

	sort_bam_file = os.path.join(config.BAM_FOLDER, haplotype+kir_name+"_sort_"+file_sufix+".bam")
	pysam.sort("-o", sort_bam_file,  bam_file)


	return kir_name, sort_bam_file


""""
	conver to BAM (binary version) because it is faster
	then run statistics,
	then eva
"""
def first_iteration(haplotype, aligments_file_haplotype):

	haplotype_gen_statistics = {}
	
	for sam_file in aligments_file_haplotype:
		kir_name, sort_bam_file =  prepare_sam_bam_file(sam_file, haplotype, "exp2_2")
		
		#count statistics
		print("counting ", haplotype, kir_name, "...")
		haplotype_gen_statistics = count_statistic(sort_bam_file)

	result_file = os.path.join(config.ALELS_STATISTICS_FOLDER, haplotype+"_exp2.txt")

	alels_statistics = evaluate_statistics(haplotype_gen_statistics, result_file)

	pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, haplotype+"_"+"exp2.pyc")
	with open(pyc_file, 'wb') as handle:
		pickle.dump(alels_statistics, handle, protocol=pickle.HIGHEST_PROTOCOL)	 

	return haplotype, alels_statistics


"""
	Prepare dictionary for create new dictionary.
"""
def create_alels_dictionary():
	# make dictionary KIR:KIR00978 => KIR2DL1*0010102
	ref_gen_files = [f for f in os.listdir(config.REFERENCE_KIR_GENS_FOLDER) if os.path.isfile(os.path.join(config.REFERENCE_KIR_GENS_FOLDER, f)) and f.endswith(END_REFERENCE_GEN_FILE) ]
	print("Preparing alels from files ", ref_gen_files)
	KIR_alels_dictionary = dict()
	file_content = ""

	for ref_file in ref_gen_files:
		with open(os.path.join(config.REFERENCE_KIR_GENS_FOLDER, ref_file)) as fileHandle:
			file_content = fileHandle.read()

		alels = file_content.split('>')

		for alel in alels[1:]: # alels[0] is empty
			# KIR:KIR00037 KIR2DS2*0010101 14577 bp
			alel_head, alel_body = alel.split('\n', 1) 
			alel_marker = alel_head.split()[0]
			alel_number = alel_marker.split(":")[1]

			alel_number = alel_number.replace('KIR', '')

			# have to put back >
			KIR_alels_dictionary[alel_number] = '>'+alel

	return KIR_alels_dictionary


"""
	Make a KIR dictionary. 
	Then find KIR gen in dictionary
	try to find which alels start with wanted alels number
"""
def create_new_reference(haplotype: list, haplotype_name: str):
	KIR_alels_dictionary = create_alels_dictionary()
	#print(KIR_alels_dictionary)

	result_haplotype = ""
	result_haplotype_legend = ""

	for item in haplotype:	
		found = False
		# KIR:KIR00432
		wanted_alel = item.split(":")[1]

		wanted_alel = wanted_alel.replace('KIR', '')
		wanted_alel = wanted_alel.strip()

		if(wanted_alel in KIR_alels_dictionary):
			result_haplotype += KIR_alels_dictionary[wanted_alel]
		else:
			print("Can not find item ", wanted_gen+"*"+wanted_alel)
							

	output_file = os.path.join("/home/kate/Dokumenty/FAV/Diplomka/software/data/reference/", haplotype_name+"reference.fasta")
	haplotype_result_file_rename = open(output_file, "w+")
	haplotype_result_file_rename.write(result_haplotype)
	haplotype_result_file_rename.close()

	#print("Haplotype legend ", haplotype_name, ":", result_haplotype_legend)
	print("Create file: ", output_file)
	return output_file


def run():
	# get all alignment .sam
	aligments_sam_files = [f for f in os.listdir(config.ALIGNMENT_FOLDER) if os.path.isfile(os.path.join(config.ALIGNMENT_FOLDER, f)) and os.path.splitext(f)[1] == '.sam']

	haplotype_name_list_unique = find_haplotypes(aligments_sam_files)
	print("Find haplotype: ", haplotype_name_list_unique)

	for haplotype in haplotype_name_list_unique:
		aligments_file_haplotype = [f for f in aligments_sam_files if f.startswith(haplotype+"KIR_gen")]
		print(aligments_file_haplotype)

		haplotype_name, alels_statistics = first_iteration(haplotype, aligments_file_haplotype)
		
		haplotype_gens = list()
		for key, statistic in alels_statistics.items():
			if statistic['normalize_coverage'] > config.CUT_COVERAGE_ALELS:
				haplotype_gens.append(key)	

		new_reference = create_new_reference(haplotype_gens, haplotype_name)
		print("chci vytv orit index")
		#process = subprocess.run([config.BOWTIE_HOME_DIRECTORY+"/bowtie2-build", "--threads", str(config.BOWTIE_THREADS), new_reference, haplotype_name], cwd="/home/kate/Dokumenty/FAV/Diplomka/software/data/temp")
		print("aligment build index", haplotype_name)
		if haplotype_name!="amala_":
			read1= os.path.join(config.READS_FOLDER, "bob1.fq")
			read2= os.path.join(config.READS_FOLDER, "bob2.fq")
			process = subprocess.run([config.BOWTIE_HOME_DIRECTORY+"/bowtie2", "--threads", str(config.BOWTIE_THREADS), "-x", haplotype_name, "-1", read1,"-2", read2, "-S", os.path.join("/home/kate/Dokumenty/FAV/Diplomka/software/data/temp", haplotype_name+"align.sam")], cwd="/home/kate/Dokumenty/FAV/Diplomka/software/data/temp") 	

	#print(haplotype_gens)
	# first cut coverage
	#


