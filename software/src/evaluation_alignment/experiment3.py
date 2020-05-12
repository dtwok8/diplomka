#import pysam
import sys
import os
import pickle
import copy

#import Levenshtein

import config

END_REFERENCE_GEN_FILE = "_gen.fasta"


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

			KIR_alels_dictionary[alel_marker] = alel_body

	return KIR_alels_dictionary


def run():
	print("nahdow")
	KIR_alels_dictionary = create_alels_dictionary()

	with open(config.ALELS_STATISTICS_FILE_PYC, 'rb') as handle:
		alels_statistics = pickle.load(handle)

	with open(config.ALELS_DISTANCE_FILE_PYC, 'rb') as handle:
		alels_distance = pickle.load(handle)

	print(len(alels_statistics))
	exit(1)

	cut_alel_statistics = dict()
	for key, statistics in alels_statistics.items():
		if(statistics['normalize_coverage'] > config.CUT_COVERAGE_ALELS):
			cut_alel_statistics[key] = statistics

	print(len(cut_alel_statistics))

	cut_alel_statistics_deep_copy = copy.deepcopy(cut_alel_statistics)

	for key1, statistics1 in cut_alel_statistics.items():
		for key2, statistics2 in cut_alel_statistics.items():
			if(key1==key2):
				continue
				
			#get distance
			#distance = Levenshtein.distance(KIR_alels_dictionary[key1], KIR_alels_dictionary[key2])

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
					

				print("key1 ", key1_coverage_changes, "key2", key2_coverage_changes)
				
				if(abs(key1_coverage_changes - key2_coverage_changes) > 100):
					if(2*key1_coverage_changes < key2_coverage_changes):
						# delete key1
						if(key1 in cut_alel_statistics_deep_copy):
							del cut_alel_statistics_deep_copy[key1]
							print("delete", key1)
					elif(2*key2_coverage_changes < key1_coverage_changes):
						if(key2 in cut_alel_statistics_deep_copy):
							del cut_alel_statistics_deep_copy[key2]
							print("delete", key2)
						# delete key2
					#print(type(changes))

	print(len(cut_alel_statistics_deep_copy))

	for key1, statistics1 in cut_alel_statistics_deep_copy.items():
		print(key1)


	# musim vzit kazde dvě alelely a zase si urcit vzdalenost mezi nimi 
	# a od nejake vzdalenosti se mrknout jak jsou zarovanany tam kde je ten rozdil mezi nimi

