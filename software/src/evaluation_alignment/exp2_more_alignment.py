import pysam
import sys
import os
import pickle
import copy

import config
import src.evaluation_alignment.make_new_reference_and_align as make_new_reference_and_align


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

	return alels_statistics


"""
	Evaluate statistics, find max coverage of each gene
	and write it into result file.
"""
def evaluate_statistics(haplotype_gen_statistics: dict, output_file_name: str):

	for alel, values in haplotype_gen_statistics.items():
		values["normalize_coverage"] = values['width_coverage'] / (values['alel_size'] / 100)
			
	# asi muzeme predpokladat za budou serazeny
	# takze by to asi jen chtelo seradit
	# {'KIR:KIR00001': {'number_nuc_coverage': 7407, 'sum_nuc_coverage': 19500, 'alel_size': 14738, 'max_histogram': 6, 'normalize_coverage': 50.25783688424481}, 'KIR:KIR00978': 
	#sorted_x = sorted(haplotype_gen_statistics.items(), key=operator.itemgetter(1["normalize_coverage"]))
	
	sorted_x = {k: v for k, v in sorted(haplotype_gen_statistics.items(), key=lambda item: item[1]["normalize_coverage"], reverse= True)}
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


def create_statistics(haplotype_aligment_file, basic_name: str, file_sufix: str):
	sort_bam_file = prepare_sam_bam_file(haplotype_aligment_file, file_sufix)
	haplotype_gen_statistics = count_statistic(sort_bam_file)

	result_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+file_sufix+".txt")

	alels_statistics = evaluate_statistics(haplotype_gen_statistics, result_file)

	pyc_file = os.path.join(config.ALELS_STATISTICS_FOLDER, basic_name+file_sufix+".pyc")
	
	with open(pyc_file, 'wb') as handle:
		pickle.dump(alels_statistics, handle, protocol=pickle.HIGHEST_PROTOCOL)

	print("Create result file ", pyc_file)
	return alels_statistics
		

def create_alels_dictionary():
	print("Creating dictionary...")
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
			alel_marker = alel_head.split()[0] # KIR:KIR00037

			KIR_alels_dictionary[alel_marker] = alel_marker = alel_head.split()[1] 
	return KIR_alels_dictionary


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

def cut_avg(alels_statistics):
	# find piky
	# prumer v ramci genu a 2x prumer by mohl byt lepsi extrem
	# takze kdyz to tam najit tak to tam nechat a zbytek co? 
	# potrebuju to projit po genech
	#for key, alel in alels_statistics.items():
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
			# kdyz je v tom genu nekde vetsi pik nez je prumer
			# tak ty ktery to nemaji zahodit.. 
			#just same gen 
			gen1 = statistics1['translation'].split('*')[0]
			gen2 = statistics2['translation'].split('*')[0]

			if(gen1 == gen2):
				if(gen1 == "KIR2DL1"):
					print("------------------------------------------------------------")
				# print(gen1,  statistics1['translation'], statistics1['max_histogram'], average_deep_cov[gen1]['avg'],  statistics2['translation'], statistics2['max_histogram'] )
				# if(statistics1['max_histogram'] > 3*average_deep_cov[gen1]['avg']  and statistics2['max_histogram'] > 3*average_deep_cov[gen2]['avg'] ):
				# 	continue
				# elif(statistics1['max_histogram'] > 3 * average_deep_cov[gen1]['avg'] ):
				# 	if(key2 in alels_statistics_copy):
				# 		print("del avg ",  statistics2['translation'])
				# 		del alels_statistics_copy[key2]
				# elif(statistics2['max_histogram'] > 3 *average_deep_cov[gen2]['avg'] ):
				# 	if(key1 in alels_statistics_copy):
				# 		print("del avg",  statistics1['translation'])
				# 		del alels_statistics_copy[key1]

				if(statistics1['max_histogram'] > 2*statistics2['max_histogram'] ):
					if(statistics1['max_histogram'] > 2 * average_deep_cov[gen1]['avg']):
						if(key2 in alels_statistics_copy):
							print("del avg ",  statistics2['translation'])
							del alels_statistics_copy[key2]		

				elif(statistics2['max_histogram'] > 2*statistics1['max_histogram']):
					if(statistics2['max_histogram'] > 2 * average_deep_cov[gen2]['avg']):
						if(key1 in alels_statistics_copy):
							print("del avg ",  statistics1['translation'])
							del alels_statistics_copy[key1]	
						


				#del cut_alel_statistics_deep_copy[key1]

				#['deep_coverage'] = deep_coverage
			#alels_statistics[curent_alel]['alel_size'] = alel_size
	return  alels_statistics_copy


def run():
	all_references_alels = make_new_reference_and_align.get_all_references_alels()
	# get all alignment .sam
	aligments_sam_files = [f for f in os.listdir(config.ALIGNMENT_FOLDER) if os.path.isfile(os.path.join(config.ALIGNMENT_FOLDER, f)) and os.path.splitext(f)[1] == '.sam']

	print("Aligment files: ", aligments_sam_files)

	for haplotype_aligment_file in aligments_sam_files:
		# step 1 
		basic_name = haplotype_aligment_file.split('.')[0] 
		alels_statistics = create_statistics(haplotype_aligment_file, basic_name, "_exp2_step1")

		# step 2
		haplotype_gens = list()
		for key, statistic in alels_statistics.items():
			if statistic['normalize_coverage'] > config.CUT_COVERAGE_ALELS:
				haplotype_gens.append(key)

		new_reference = make_new_reference_and_align.create_new_reference(haplotype_gens, basic_name, all_references_alels, "_exp2_step2")
		new_align_file = make_new_reference_and_align.align("/home/kate/Dokumenty/FAV/Diplomka/software/data/temp", new_reference, basic_name, "_exp2_step2")

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