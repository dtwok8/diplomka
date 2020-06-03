#import pysam
import sys
import os
import pickle
import copy

import numpy as np
import Levenshtein
import matplotlib.pyplot as plt

GEN_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"
END_REFERENCE_GEN_FILE = "_gen.fasta"
ALELS_STATISTICS_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/alels_statistics.pyc"
ALELS_DISTANCE_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/alels_distance_full.pyc"
PLOT_OUTPUT_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/analysis_after_align_result"

# ANALYSIS IN ALIGMENT
ALELS_IN_ALIGMENT = ['3DL3: 0040201', '3DL3:00802', '2DS2: 0010101', '2DL2: 0030102', '2DL3: 0010109', '2DP1: 0020108', '2DL1: 0030201', '3DP1: 007', '3DP1: 0090101', '2DL4: 0010201', '2DL4: 00501', '3DL1: 0150201', '3DS1: 0130101', '2DL5A: 00102', '2DS5: 0020101', '2DS1: 0020106', '2DS4: 0010101', '3DL2: 0020105', '3DL2:0070102']

"""
	make dictionary KIR:KIR00978 => KIR2DL1*0010102
"""
def create_alels_translate_dictionary():
	KIR_dictionary = dict()

	with open(GEN_FILE, "r") as openfileobject:
		for line in openfileobject:
			if(line.startswith('>KIR')):
				split_line = line.split()
				# split_line[0][1:] - cut '>'
				KIR_dictionary[split_line[0][1:]] = split_line[1]

	print("make dictionary: ", GEN_FILE)
	
	return KIR_dictionary


def run():
	KIR_alels_dictionary = create_alels_translate_dictionary()

	with open(ALELS_STATISTICS_FILE_PYC, 'rb') as handle:
		alels_statistics = pickle.load(handle)

	with open(ALELS_DISTANCE_FILE_PYC, 'rb') as handle:
		alels_distance = pickle.load(handle)

	# convert to right format
	alels_in_aligmnet_format = []
	for item in ALELS_IN_ALIGMENT:
		item = item.replace(':', '*')
		gen, alel = item.split('*') 
		if('KIR' in item):
			alels_in_aligmnet_format.append(gen.strip()+"*"+alel.strip())
		else:
			alels_in_aligmnet_format.append("KIR"+gen.strip()+"*"+alel.strip())
	alels_statistics_sort_by_gens = dict()
	for key, statistics in alels_statistics.items():
		trasnlate_alel = KIR_alels_dictionary[key] # translate
		gen = trasnlate_alel.split('*')[0] 
		
		if(gen not in alels_statistics_sort_by_gens):
			alels_statistics_sort_by_gens[gen] = dict()
			
		alels_statistics_sort_by_gens[gen][trasnlate_alel] = statistics


	for gen, gen_alel_statistic in alels_statistics_sort_by_gens.items():
		fig, axs = plt.subplots(len(gen_alel_statistic))
		alel_number = 0

		max_lenght = 0
		max_height = 0
		for alel, statistics in gen_alel_statistic.items():
			if(max_lenght < len(statistics['cov_nucleotids'])):
				max_lenght = len(statistics['cov_nucleotids'])

			for i in range(0, len(statistics['cov_nucleotids'])):	
				if(max_height < statistics['cov_nucleotids'][i]):
					max_height = statistics['cov_nucleotids'][i]

				# because we want to see from plot that there is no coverage
				# plot not display nan
				if(statistics['cov_nucleotids'][i]==0):
					#statistics['cov_nucleotids'][i]=float('nan')
					statistics['cov_nucleotids'][i]=0

		for alel, statistics in gen_alel_statistic.items():	
			xf = np.linspace(0, len(statistics['cov_nucleotids']), len(statistics['cov_nucleotids']))

			#fig.suptitle(alel)
			# set x "label" just for last alel
			if(alel != list(gen_alel_statistic.keys())[-1]):
				axs[alel_number].set_xticks([])

			axs[alel_number].set_ylim(0, max_height+10)
			axs[alel_number].set_xlim(0, max_lenght)

			# if alel is in haplotype which was aligment
			make_red = False
			for item in alels_in_aligmnet_format:
				if(item in alel):
					make_red = True
					break

			if(make_red):
				axs[alel_number].plot(xf, statistics['cov_nucleotids'], '-r')
				axs[alel_number].fill_between(xf, 0, statistics['cov_nucleotids'], facecolor='red')
			else:	
				axs[alel_number].plot(xf, statistics['cov_nucleotids'], '-b')
				axs[alel_number].fill_between(xf, 0, statistics['cov_nucleotids'], facecolor='blue')
			
			axs[alel_number].set_aspect('auto', 'box')
			axs[alel_number].text(max_lenght+10, max_height*0.5, alel.split('*')[1]+", cov: "+str(round(statistics['normalize_coverage'], 2))+"%, len: "+str(len(statistics['cov_nucleotids'])) , fontsize=9)

			# so slow
			#axs[alel_number].bar(xf, statistics['cov_nucleotids'], color='#bd2a47', width=1, linewidth='0.1',  edgecolor='white', label="a")
			alel_number+=1
			print(alel)

		#plt.tight_layout()
		fig.suptitle(gen)

		plt.gcf().set_size_inches(max_lenght/1000, len(gen_alel_statistic))
		#plt.gcf().set_size_inches(20, len(alels_statistics_sort_by_gens))
		#plt.gcf().set_size_inches(len(alels_statistics_sort_by_gens), 9)
		#plt.savefig(os.path.join(PLOT_OUTPUT_FOLDER, gen+"compare_coverage.svg"), format='svg', bbox_inches='tight', aspect='auto', dpi=3000)
		plt.savefig(os.path.join(PLOT_OUTPUT_FOLDER, gen+"compare_coverage.svg"), format='svg', bbox_inches='tight', dpi=3000)
		plt.clf() # clear plot



	# musim vzit kazde dvě alelely a zase si urcit vzdalenost mezi nimi 
	# a od nejake vzdalenosti se mrknout jak jsou zarovanany tam kde je ten rozdil mezi nimi

run()