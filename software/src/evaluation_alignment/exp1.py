import pysam
import sys
import os
import pickle
import copy

import config


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
			
	# asi muzeme predpokladat za budou serazeny
	# takze by to asi jen chtelo seradit
	# {'KIR:KIR00001': {'number_nuc_coverage': 7407, 'sum_nuc_coverage': 19500, 'alel_size': 14738, 'max_histogram': 6, 'normalize_coverage': 50.25783688424481}, 'KIR:KIR00978': 
	#sorted_x = sorted(genotype_gen_statistics.items(), key=operator.itemgetter(1["normalize_coverage"]))
	
	sorted_x = {k: v for k, v in sorted(genotype_gen_statistics.items(), key=lambda item: item[1]["normalize_coverage"], reverse= True)}
	# takze nejdriv je seradit podle coverage
	output_file = open(output_file_name, "w")

	for alel, values in sorted_x.items():
		#if(values['number_nuc_coverage'] > 0):
		print(alel, " c: ", values['normalize_coverage'], " deep coverage: ", values['deep_coverage'], "max histogram: ", values['max_histogram'] , file=output_file)

		#print(alel, " c: ", values['normalize_coverage'], "max_value_in_histogram", values['max_histogram'], file=output_file)

	output_file.close()
	print("create result file ", output_file_name)
	return sorted_x


def prepare_sam_bam_file(sam_file: str, file_sufix: str):
	print("creating sam/bam file for ", sam_file, "...")
	# create bam file and  conver sam to bam and then sort 
	bam_file = os.path.join(config.BAM_FOLDER, sam_file+".bam")
	fh = open(bam_file, 'w')
	fh.close()

	pysam.view('-o', bam_file, '-b', os.path.join(config.ALIGNMENT_FOLDER, sam_file), save_stdout=bam_file)

	sort_bam_file = os.path.join(config.BAM_FOLDER, sam_file+"_sort"+file_sufix+".bam")
	pysam.sort("-o", sort_bam_file,  bam_file)

	return sort_bam_file


def create_statistics(genotype_aligment_file, basic_name: str, file_sufix: str):
	sort_bam_file = prepare_sam_bam_file(genotype_aligment_file, file_sufix)
	genotype_gen_statistics = count_statistic(sort_bam_file)

	result_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+file_sufix+".txt")

	alels_statistics = evaluate_statistics(genotype_gen_statistics, result_file)

	pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+file_sufix+".pyc")
	
	with open(pyc_file, 'wb') as handle:
		pickle.dump(alels_statistics, handle, protocol=pickle.HIGHEST_PROTOCOL)

	print("Create result file ", pyc_file)
	return alels_statistics


"""
	cut similar alels in one gen
"""
def cut_similar_alels(alels_statistics):
	KIR_alels_dictionary = create_alels_dictionary()
	

	with open(config.ALELS_DISTANCE_FILE_PYC, 'rb') as handle:
		alels_distance = pickle.load(handle)

	# add translation na KIR3DL1*000120
	for key, statistics in alels_statistics.items():
		alels_statistics[key]['translation'] = KIR_alels_dictionary[key]

	cut_alel_statistics_deep_copy = copy.deepcopy(alels_statistics)

	#print(alels_distance['KIR2DL4*042']['KIR3DL2*0070102'])

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
							# delete key2
						#print(type(changes))

	return cut_alel_statistics_deep_copy

def create_alels_dictionary():
	print("Creating dictionary...")
	# make dictionary KIR:KIR00978 => KIR2DL1*0010102
	print("Preparing alels from files ", config.REFERENCE_KIR_GENS_FILE)
	KIR_alels_dictionary = dict()
	file_content = ""

	with open(config.REFERENCE_KIR_GENS_FILE) as fileHandle:
		file_content = fileHandle.read()

	alels = file_content.split('>')

	for alel in alels[1:]: # alels[0] is empty
		# KIR:KIR00037 KIR2DS2*0010101 14577 bp
		alel_head, alel_body = alel.split('\n', 1)
		alel_marker = alel_head.split()[0] # KIR:KIR00037

		KIR_alels_dictionary[alel_marker] = alel_marker = alel_head.split()[1] 
	return KIR_alels_dictionary


""""
	conver to BAM (binary version) because it is faster
	then run statistics,
	then eva
"""
def run():
	# get all alignment .sam
	aligments_sam_files = [f for f in os.listdir(config.ALIGNMENT_FOLDER) if os.path.isfile(os.path.join(config.ALIGNMENT_FOLDER, f)) and os.path.splitext(f)[1] == '.sam']

	print("Aligment files: ", aligments_sam_files)

	for genotype_aligment_file in aligments_sam_files:
		basic_name = genotype_aligment_file.split('.')[0] 
		alels_statistics = create_statistics(genotype_aligment_file, basic_name, "_exp1_step1")	

		statistic_cov_cut = dict()
		for alel, statistic in alels_statistics.items():
			if(statistic['normalize_coverage'] > config.CUT_COVERAGE_ALELS):
				print(alel, statistic['normalize_coverage'])
				statistic_cov_cut[alel] = statistic

		pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+"_exp1_step2.pyc")
		with open(pyc_file, 'wb') as handle:
			pickle.dump(statistic_cov_cut, handle, protocol=pickle.HIGHEST_PROTOCOL)


		res = cut_similar_alels(statistic_cov_cut)
		pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+"_exp1_step3.pyc")
		with open(pyc_file, 'wb') as handle:
			pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)