import pysam
import sys
import os
import pickle
import copy

import config
import src.alleles_identification.make_new_reference_and_align as make_new_reference_and_align


""" Because there can be more files, like nuc and another"""
END_REFERENCE_GEN_FILE = "_gen.fasta"


"""
	Make statistics like coverage of allels and sum of histogram covarage - deep_coverage.
"""
def count_statistic(sort_bam_file: str):
	coverage_string = pysam.depth("-aa", sort_bam_file)
	coverage_list= coverage_string.splitlines()
	curent_alel =""
	# number of positions covered 
	width_coverage = 0
	deep_coverage = 0
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
				alels_statistics[curent_alel]['width_coverage'] = width_coverage
				alels_statistics[curent_alel]['deep_coverage'] = deep_coverage
				alels_statistics[curent_alel]['alel_size'] = alel_size
				alels_statistics[curent_alel]['max_histogram'] = max_histogram
				alels_statistics[curent_alel]['cov_nucleotids'] = cov_nucleotids


			curent_alel = item_nuc_list[0]
			width_coverage = 0	
			deep_coverage = 0
			alel_size = 0
			max_histogram = 0
			cov_nucleotids = list()

		alel_size += 1
		cov_nucleotids.append(int(item_nuc_list[2]));

		if(int(item_nuc_list[2]) != 0):
			width_coverage += 1
			deep_coverage += int(item_nuc_list[2])

			if(max_histogram < int(item_nuc_list[2])):
				max_histogram = int(item_nuc_list[2])

	# last item
	alels_statistics[curent_alel] = {}
	alels_statistics[curent_alel]['width_coverage'] = width_coverage
	alels_statistics[curent_alel]['deep_coverage'] = deep_coverage
	alels_statistics[curent_alel]['alel_size'] = alel_size
	alels_statistics[curent_alel]['max_histogram'] = max_histogram
	alels_statistics[curent_alel]['cov_nucleotids'] = cov_nucleotids

	print("Found alels in aligments: ", len(alels_statistics))
	return alels_statistics


"""
	Evaluate statistics, find max coverage of each gene
	and write it into result file.
"""
def evaluate_statistics(genotype_gen_statistics: dict, output_file_name: str):
	for alel, values in genotype_gen_statistics.items():
		values["normalize_coverage"] = values['width_coverage'] / (values['alel_size'] / 100)
			

	# {'KIR:KIR00001': {'number_nuc_coverage': 7407, 'sum_nuc_coverage': 19500, 'alel_size': 14738, 'max_histogram': 6, 'normalize_coverage': 50.25783688424481}, 'KIR:KIR00978': 
	sorted_x = {k: v for k, v in sorted(genotype_gen_statistics.items(), key=lambda item: item[1]["normalize_coverage"], reverse= True)}
	# takze nejdriv je seradit podle coverage
	output_file = open(output_file_name, "w")

	for alel, values in sorted_x.items():
		print(alel, " c: ", values['normalize_coverage'], " deep coverage: ", values['deep_coverage'], "max histogram: ", values['max_histogram'] , file=output_file)


	output_file.close()
	print("create result file ", output_file_name)
	return sorted_x


"""
	Make sort sam file and then bam file. Bam is binary sam. 
"""
def prepare_sam_bam_file(sam_file: str, basic_name:str, file_sufix: str):
	print("creating sam/bam file for ", sam_file, "...")
	# create bam file and  conver sam to bam and then sort 
	bam_file = os.path.join(config.BAM_FOLDER, basic_name+".bam")
	fh = open(bam_file, 'w')
	fh.close()

	print("bam_file ", bam_file)
	pysam.view('-o', bam_file, '-b', os.path.join(config.ALIGNMENT_FOLDER, sam_file), save_stdout=bam_file)

	sort_bam_file = os.path.join(config.BAM_FOLDER, basic_name+"_sort"+file_sufix+".bam")
	pysam.sort("-o", sort_bam_file,  bam_file)

	return sort_bam_file


"""
	Prepare file,
	then run evaluate statistics
	and put result into pyc file.
"""
def create_statistics(genotype_aligment_file, basic_name: str, file_sufix: str):
	sort_bam_file = prepare_sam_bam_file(genotype_aligment_file, basic_name,file_sufix)
	genotype_gen_statistics = count_statistic(sort_bam_file)

	result_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+file_sufix+".txt")

	alels_statistics = evaluate_statistics(genotype_gen_statistics, result_file)

	pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+file_sufix+".pyc")
	
	with open(pyc_file, 'wb') as handle:
		pickle.dump(alels_statistics, handle, protocol=pickle.HIGHEST_PROTOCOL)

	print("Create result file ", pyc_file)
	return alels_statistics
		

"""
	make dictionary KIR:KIR00978 => KIR2DL1*0010102
"""
def get_translate_dictionary():
	KIR_dictionary = dict()


	with open(config.REFERENCE_KIR_GENS_FILE, "r") as openfileobject:
		for line in openfileobject:
			if(line.startswith('>KIR')):
				split_line = line.split()
				# split_line[0][1:] - cut '>'
				KIR_dictionary[split_line[0][1:]] = split_line[1]

	print("make dictionary: ", config.REFERENCE_KIR_GENS_FILE)
	
	return KIR_dictionary


"""
	cut similar alels in one gen
	compare distance beetwen two alels in one gen,
	and where they are close enought then compare deep coverage on place where are not same.
	Alels with smaller deep coverage will be delete.
"""
def cut_similar_alels(alels_statistics):
	KIR_alels_dictionary = get_translate_dictionary()
	

	with open(config.ALELS_DISTANCE_FILE_PYC, 'rb') as handle:
		alels_distance = pickle.load(handle)

	# add translation na KIR3DL1*000120
	for key, statistics in alels_statistics.items():
		alels_statistics[key]['translation'] = KIR_alels_dictionary[key]

	cut_alel_statistics_deep_copy = copy.deepcopy(alels_statistics)

	for key1, statistics1 in alels_statistics.items():
		for key2, statistics2 in alels_statistics.items():
			if(key1==key2):
				continue
			
			#just same gen 
			gen1 = statistics1['translation'].split('*')[0]
			gen2 = statistics2['translation'].split('*')[0]

			if(gen1 == gen2):

				#get distance - because it same across diagonal, so it was save a count just one times
				if(key2 in alels_distance[key1]):
					distance = alels_distance[key1][key2]['distance']
					changes = alels_distance[key1][key2]['changes']
				elif(key1 in alels_distance[key2]):
					distance = alels_distance[key2][key1]['distance']
					changes = alels_distance[key2][key1]['changes']
				else:
					print("Not found distance", key1, key2)


				if(distance < config.CLOSE_DISTANCE):	
					#changes = Levenshtein.editops(KIR_alels_dictionary[key1], KIR_alels_dictionary[key2])
					key1_coverage_changes = 0
					key2_coverage_changes = 0

					for one_change in changes:
						# it sometimes touch behind string
						if(len(alels_statistics[key1]['cov_nucleotids']) > one_change[1]):
							key1_coverage_changes+=alels_statistics[key1]['cov_nucleotids'][one_change[1]]

						if(len(alels_statistics[key2]['cov_nucleotids']) > one_change[2]):
							key2_coverage_changes+=alels_statistics[key2]['cov_nucleotids'][one_change[2]]
						
					if(2*key1_coverage_changes < key2_coverage_changes):
						# delete key1
						if(key1 in cut_alel_statistics_deep_copy):
							del cut_alel_statistics_deep_copy[key1]
							print("delete", statistics1['translation'], "because ", statistics2['translation'])
					elif(2*key2_coverage_changes < key1_coverage_changes):
						if(key2 in cut_alel_statistics_deep_copy):
							del cut_alel_statistics_deep_copy[key2]
							print("delete", statistics2['translation'], "because ", statistics1['translation'])

	return cut_alel_statistics_deep_copy


"""
	Find peak in gene
		compare two alels max_histogram_alel1 > 2*max_histogram_alel2
		2*average_all_alles_in_gene < max_histrogram_alel1
		then delete alel2
"""
def cut_avg(alels_statistics):
	alels_statistics_copy = copy.deepcopy(alels_statistics)
	average_deep_cov = dict()
	for key1, statistics in alels_statistics.items():
		gen = statistics['translation'].split('*')[0]

		if gen not in average_deep_cov:
			average_deep_cov[gen] = {
										'sum': statistics['deep_coverage']/statistics['alel_size'],
										'count': 1 
									}

		else:
			average_deep_cov[gen]['sum'] +=	statistics['deep_coverage']/statistics['alel_size']
			average_deep_cov[gen]['count'] += 1 


	for key, item in average_deep_cov.items():
		average_deep_cov[key]['avg'] = item['sum'] /item['count'] 


	for key1, statistics1 in alels_statistics.items():
		for key2, statistics2 in alels_statistics.items():
			if(key1==key2):
				continue

			#just same gen 
			gen1 = statistics1['translation'].split('*')[0]
			gen2 = statistics2['translation'].split('*')[0]

			if(gen1 == gen2):
				if(statistics1['max_histogram'] > 2*statistics2['max_histogram'] ):
					if(statistics1['max_histogram'] > 2 * average_deep_cov[gen1]['avg']):
					#if(statistics1['deep_coverage']/statistics1['alel_size']> 2 * average_deep_cov[gen1]['avg']): # dont delete peak
						if(key2 in alels_statistics_copy):
							print("del avg ",  statistics2['translation'])
							del alels_statistics_copy[key2]		

				elif(statistics2['max_histogram'] > 2*statistics1['max_histogram']):
					if(statistics2['max_histogram'] > 2 * average_deep_cov[gen2]['avg']):
						if(key1 in alels_statistics_copy):
							print("del avg ",  statistics1['translation'])
							del alels_statistics_copy[key1]	
						
	return  alels_statistics_copy


"""
	make dictionary KIR:KIR00978 => KIR2DL1*0010102
"""
def get_translate_dictionary():
	# make dictionary KIR:KIR00978 => KIR2DL1*0010102

	KIR_dictionary = dict()

	with open(config.REFERENCE_KIR_GENS_FILE, "r") as openfileobject:
		for line in openfileobject:
			if(line.startswith('>KIR')):
				split_line = line.split()
				# split_line[0][1:] - cut '>'
				KIR_dictionary[split_line[0][1:]] = split_line[1]

	print("make dictionary: ", config.REFERENCE_KIR_GENS_FILE)
	
	return KIR_dictionary


def run():
	all_references_alels = make_new_reference_and_align.get_all_references_alels()
	# get all alignment .sam
	aligments_sam_files = [f for f in os.listdir(config.ALIGNMENT_FOLDER) if os.path.isfile(os.path.join(config.ALIGNMENT_FOLDER, f)) and os.path.splitext(f)[1] == '.sam']

	print("Aligment files: ", aligments_sam_files)
	translate_dictionary = get_translate_dictionary()

	for genotype_aligment_file in aligments_sam_files:
		# step 1 
		basic_name = genotype_aligment_file.split('.')[0] 
		alels_statistics = create_statistics(genotype_aligment_file, basic_name, "_exp2_step1")

		# step 2
		genotype_gens = list()
		for key, statistic in alels_statistics.items():
			if statistic['normalize_coverage'] > config.CUT_COVERAGE_ALLELES:
				genotype_gens.append(key)

		new_reference = make_new_reference_and_align.create_new_reference(genotype_gens, basic_name, all_references_alels, "_exp2_step2")
		new_align_file = make_new_reference_and_align.align(config.TEMP_FOLDER, new_reference, basic_name, "_exp2_step2")

		#new_align_file = os.path.join("/home/kate/Dokumenty/FAV/Diplomka/software/data/temp", basic_name+"_exp2_step2.sam")
		alels_statistics = create_statistics(new_align_file, basic_name, "_exp2_step2")
		
		# step 3 
		# check distance beetween alels
		res = cut_similar_alels(alels_statistics)
		pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+"_exp2_step3.pyc")
		with open(pyc_file, 'wb') as handle:
			pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)


		# step 4
		res = cut_avg(res)

		pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+"_exp2_step4.pyc")
		with open(pyc_file, 'wb') as handle:
			pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)

		# result to result file
		result_txt_file_name = os.path.join(config.RESULT_FOLDER, basic_name+"_exp2"+".txt")
		result_txt_file = open(result_txt_file_name, "w")

		for alel, values in res.items():
			print(translate_dictionary[alel], " c: ", values['normalize_coverage'], file=result_txt_file)

		result_txt_file.close()