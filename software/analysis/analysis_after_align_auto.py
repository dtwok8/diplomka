#import pysam
import sys
import os
import pickle
import copy

import numpy as np
import Levenshtein
import matplotlib.pyplot as plt

""" Reference file """
GEN_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"

""" Make a file name from this and GENOMES_LIST """
# example ALELS_STATISTICS_PYC_FOLDER/GENOMES_LIST[1]_ALELS_STATISTICS_PYC_REFERENCE_NAME_ALELS_STATISTICS_PYC_EXPERIMENT
ALELS_STATISTICS_PYC_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/statistics/"
ALELS_STATISTICS_PYC_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/statistics/"
ALELS_STATISTICS_PYC_REFERENCE_NAME = "KIR_gen"
ALELS_STATISTICS_PYC_EXPERIMENT = "exp3"

""" Steps whitch you want to analyse """
STEPS = ["step1", "step2", "step3", "step4"]
""" Steps from which will be create a latex table """
STEPS_LATECH_TABLE = ["step2", "step3","step4"]

PLOT_OUTPUT_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/analysis_after_align_result"


""" which genome want to analyse """ 
# syntetic
GENOMES_LIST = ["amala", "bob", "cox", "ho301", "jvm", "kas011", "olga", "rsh", "wt51", "test1", "test2", "test3", "test4", "test5", "test6", "test7", "test8", "test9", "test10", "test11"]

# real 
#GENOMES_LIST = ["amala", "bob", "cox", "ho301", "jvm", "kas011", "olga", "rsh", "wt51"]

""" reference """
# real
"""
GENOMES_ALLELES = {
			"amala": [
 				"2DL1: 00302", "2DL2: 00301", "2DL3: 001", "2DL4: 00102", "2DL4: 00501", "2DL5A: 001", "2DP1: 00201", "2DS1: 00201", "2DS2: 00101", "2DS4: 001", "2DS5: 00201", "3DL1: 01502", "3DL2: 0020105", "3DL2: 0070102",
 				"3DL3: 00402", "3DL3: 00802", "3DP1: 007", "3DP1: 00901", "3DS1: 01301"
			],
			"bob": [
				"2DL1: 00302", "2DL2: 00301", "2DL3: 00201", "2DL4: 001", "2DL4: 005", "2DL5A: 00101", "2DP1: 00301", "2DS1: 00201", "2DS2: 00101", "2DS4: 001", "2DS5: 00201", "3DL1: 002", "3DL2: 0020101",  "3DL2: 0070102",
				"3DL3: 00101", "3DL3: 01303", "3DP1: 002", "3DP1: 00302", "3DS1: 01301"
			],
			"cox": [
				"2DL1: 00201", "2DL3: 00201", "2DL3: 007", "2DL4: 00501", "2DL4: 011", "2DL5A: 00101", "2DP1: 00301", "2DS1: 00201", "2DS4: 010", "2DS5: 00201", "3DL1: 005010", "3DL2: 00103",  "3DL2: 007", "3DL3: 00102", "3DL3: 00103",
				"3DP1: 005", "3DP1: 006", "3DS1: 055"
			],
			"ho301": [
				"2DL1: 004", "2DL1: 010", "2DL2: 00101", "2DL2: 00301", "2DL4: 00102", "2DL5B: 010", "2DP1: 00102", "2DS2: 00101", "2DS2: 002", "2DS3: 00103", "2DS3: 00201", "2DS3: 00103", "2DS3: 00201", "2DS4: 001", "3DL1: 002", 
				"3DL2: 00201", "3DL3: 014", "3DP1: 003010", "3DP1: 004"
			],
			"jvm": [
				"2DL1: 00302", "2DL2: 00301", "2DL4: 00103", "2DL4: 00801", "2DP1: 005", "2DS2: 00101", "2DS4: 003010", "3DL1: 00101", "3DL1: 008", "3DL2: 00101", "3DL2: 009", "3DL3: 007", "3DL3: 00801", "3DP1: 001", "3DP1: 00302"
			],
			"kas011": [
				"2DL1: 00201", "2DL1: 00302", "2DL4: 00103", "2DL4: 005", "2DL5A: 00101", "2DP1: 002", "2DP1: 00301", "2DS1: 00201", "2DS4: 00301", "2DS5: 00201", "3DL1: 008", "3DL2: 01001", "3DL2: 019", "3DL3: 00901", "3DL3: 01302",
				"3DP1: 00302", "3DP1: 006", "3DS1: 01301"	
			],
			"olga": [
				"2DL1: 00302", "2DL3: 00101", "2DL4: 005", "2DL4: 011", "2DL5A: 00103", "2DP1: 00201", "2DP1: 006", "2DS1: 002", "2DS4: 010", "2DS5: 002", "3DL1: 001", "3DL1: 00501", "3DL2: 00701", "3DL3: 00201", "3DL3: 00902", 
				"3DP1: 00302", "3DS1: 01301"
			],
			"rsh": [
				"2DL1: 00302", "2DL1: 01201", "2DL4: 0010307", "2DL4: 011", "2DL5B: 004", "2DP1: 00201", "2DP1: 009", "2DS2: 00101", "2DS5: 006", "3DL1: 00501", "3DL1: 017", "3DL3: 0040202", "3DL3: 00901", "3DP1: 00304", "3DP1: 008"
			],
			"wt51": [
				"2DL4: 00501", "2DL5A: 00101", "2DL5A: 00501", "2DL5B: 00201", "2DP1: 001", "2DP1: 004", "2DS1: 002", "2DS2: 00101", "2DS5: 002", "3DL3: 00103", "3DL3: 036", "3DS1: 01301"
			]

}

"""

# syntetic
GENOMES_ALLELES = {
			"bob": [ '3DL3: 00101', '3DL3: 019', '2DS2: 0010104', '2DL2: 0030101', '2DL3: 0020102', '2DP1: 0030101', '2DL1: 0030210', '3DP1: 002', '3DP1: 0030203', '2DL4: 0010202', '2DL4: 0050101', '3DL1: 002',
				'3DS1: 0130105', '2DL5A: 0010101', '2DS5: 0020104', '2DS1: 0020101', '2DS4: 0010105', '3DL2: 0020101' , '3DL2:0070102'
			],
			"amala": [
 				'3DL3: 0040201', '3DL3:00802', '2DS2: 0010101', '2DL2: 0030102', '2DL3: 0010109', '2DP1: 0020108', '2DL1: 0030201', '3DP1: 007', '3DP1: 0090101',
 				'2DL4: 0010201', '2DL4: 0050106', '3DL1: 0150201', '3DS1: 0130101', '2DL5A: 00102', '2DS5: 0020101', '2DS1: 0020106', '2DS4: 0010101', '3DL2: 0020105', '3DL2:0070102' 
			],
			"cox": [
				'3DL3: 00102', '3DL3: 0090101', '2DL3: 0020101', '2DL3: 006', '2DP1: 0030102', '2DP1: 0030102', '2DL1: 0020102', '2DL1: 0020102', '3DP1: 005', '3DP1: 006', '2DL4: 0050102', '2DL4: 00901', '3DL1: 0050103', 
				'3DS1: 055', '2DS5: 0020102', '2DS1: 0020105', '2DS4: 010', '3DL2: 0010301', '3DL2: 0070103'
			],
			"test1": [
				'3DL3: 0030101', '3DL3: 0140201', '2DS2: 0010104', '2DL2: 0030105', '2DL3: 0020101', '2DL5B: 01301', '2DS3: 0020102', '2DS3: 0020102', '2DP1: 0030101' , '2DP1:0010203', '2DL1: 0030203', 
				'2DL1:007', '3DP1: 004', '3DP1: 004', '2DL4: 0080104', '2DL4: 010', '3DL1: 0150101', '3DS1: 014', '2DL5A: 00102', '2DS1: 0020104', '2DS4: 0010103', '3DL2: 00501', '3DL2: 00501'
			],
			"test2": [
				'3DL3: 0090102', '3DL3: 019', '2DL2: 0010101', '2DL3: 0010111', '2DL5B: 0020101', '2DP1: 0020107', '2DP1: 0030102', '2DL1: 0020102', '2DL1: 0030210', '3DP1: 004', '3DP1: 01001', '2DL4: 0080104', '2DL4: 0080104',
				'3DL1: 0070101', '3DS1: 078', '2DL5A: 0050101', '2DS5: 010', '2DS1: 0020102', '2DS4:0010103', '3DL2: 0020101', '3DL2: 00501'	
			],
			"test3": [
				'3DL3: 005', '3DL3: 0140201', '2DL3: 0010101', '2DL3: 0020103', '2DP1: 0020108', '2DP1: 0020108', '2DL1: 0040101', '2DL1: 008', '3DP1: 0030102', '3DP1: 00902', '2DL4: 0010306', '2DL4: 0050104', '3DL1: 002', 
				'3DL1: 0040101', '3DL2: 0010301', '3DL2: 008'	
			],
			"test4": [
				'3DL3: 0030104', '3DL3: 007', '2DS2: 0010105', '2DL2: 0030101', '2DL3: 0010102', '2DP1: 008', '2DL1: 007', '3DP1: 007', '3DP1: 00902', '2DL4: 0010307', '2DL4: 0080104', '3DL1: 0150202', '3DS1: 055', '2DS5: 007', 
				'2DS1: 0020101', '2DS4: 0060101', '3DL2: 0020101', '3DL2: 00903'
			],
			"test5": [
				'3DL3: 0140202', '3DL3: 036', '2DL3: 0010109', '2DL3: 006', '2DS3: 0010301', '2DP1: 0030102', '2DP1: 009', '2DL1: 0030208', '2DL1: 00303', '3DP1: 001', '3DP1: 002', '2DL4: 0010202', '2DL4: 0010202', '3DL1: 0200101',
				'3DS1: 0130102', '2DL5A: 0010102', '2DS1: 0020105', '3DL2: 00202', '3DL2: 018'	
			],
			"test6": [
				'3DL3: 0090102', '3DL3: 0140203', '2DS2: 0010111', '2DL2: 0010105', '2DL3: 0010102', '2DL5B: 0080101', '2DS3: 0010302', '2DP1: 0020103', '2DP1: 010', '2DL1: 0030203', '2DL1: 0040102', '3DP1: 0030202', '3DP1: 0030402',
				'2DL4: 0010303', '2DL4: 00901', '3DL1: 0050102', '3DL1: 0250102', '2DS4: 0010104', '2DS4: 010', '3DL2: 0010302', '3DL2: 01001'	
			],
			"test7": [
				'3DL3: 00802', '3DL3: 0090103', '2DL3: 0010103', '2DL3: 0010108', '2DP1: 0020106', '2DP1: 004', '2DL1: 0030204', '2DL1: 0030205', '3DP1: 0030202', '3DP1: 0030202', '2DL4: 0010201', '2DL4: 0010305', '3DL1: 008',
				'3DL1: 0150203', '2DS4: 0010107', '2DS4: 0030104', '3DL2: 0020105', '3DL2: 00901'	
			],
			"test8": [
				'3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DL5B: 0070101', '2DS3: 0020101', '2DS3: 0020101', '2DS3: 0010302', '2DP1: 0030102', '2DL1: 00402', '3DP1: 0030101',
				'3DP1: 005', '2DL4: 00104', '2DL4: 0080104', '3DS1: 0130104', '3DS1: 055', '2DL5A: 0050102', '2DL5A: 0050102', '2DS1: 0020102', '2DS1: 0020105', '3DL2: 0010102', '3DL2: 0070102'
			],
			"test9": [
				'3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DL5B: 0070101', '2DS3: 0020101', '2DS3: 0020101', '2DP1: 0030102', '2DL1: 00402', '3DP1: 0030101', '3DP1: 005',
				'2DL4: 00104', '2DL4: 0080104', '3DL1: 0150208', '3DS1: 0130104', '2DL5A: 01201', '2DS1: 0020102', '2DS4: 0040101', '3DL2: 0010102', '3DL2: 0070102'	
			],
			"test10": [
				'3DL3: 0030101', '3DL3: 0140201', '2DS2: 0010104', '2DL2: 0030105', '2DL3: 0020101', '2DL5B: 01301', '2DS3: 0020102', '2DS3: 0020102', '2DP1: 0010203', '2DP1: 0030101', '2DL1: 0030203', '2DL1: 007', '3DP1: 004',
				'3DP1: 004', '2DL4: 0080104', '2DL4: 010', '3DL1: 0150101', '3DS1: 014', '2DS1: 0020104', '3DL2: 00501', '3DL2: 00501'
			],
			"test11": [
				'3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DP1: 0030102', '2DP1: 008', '2DL1: 00402', '2DL1: 00402', '3DP1: 0030101', '3DP1: 005', '2DL4: 00104', '2DL4: 0080104',
				'3DS1: 0130104', '3DS1: 055', '2DL5A: 0050102', '2DL5A: 0050102', '2DS5: 0020102', '2DS5: 0020103', '2DS1: 0020102', '2DS1: 0020105', '3DL2: 0010102', '3DL2: 0070102'	
			],
			"ho301": [
				'3DL3: 0140201', '3DL3: 0140201', '2DS2: 0010104', '2DS2: 0010106', '2DL2: 0010103', '2DL2: 0030107', '2DL5B: 010', '2DL5B: 010', '2DS3: 0010301', '2DS3: 0020103', '2DP1: 0010202', '2DP1: 0010202', '2DL1: 00402', 
				'2DL1: 010', '3DP1: 0030101', '3DP1: 004', '2DL4: 0010201', '2DL4: 0010201', '3DL1: 002', '3DL1: 002', '2DS4: 0010109', '2DS4: 0010109', '3DL2: 0020102', '3DL2: 0020106'
			],
			"jvm": [
				'3DL3: 00801', '3DL3: 0140201', '2DS2: 0010110', '2DL2: 0030102', '2DL3: 010', '2DP1: 004', '2DL1: 0030203', '3DP1: 001', '3DP1: 0030202', '2DL4: 0010304', '2DL4: 0080101', '3DL1: 0010104', '3DL1: 008', '2DS4: 0030103',
				'2DS4: 0030103', '3DL2: 0010101', '3DL2:018'	
			],
			"kas011": [
				'3DL3: 0090101', '3DL3: 0140203', '2DL3: 0020103', '2DL3: 0020103', '2DP1: 0020104', '2DP1: 0030101', '2DL1: 0020101', '2DL1: 0030209', '3DP1: 0030206', '3DP1: 009', '2DL4: 0010301', '2DL4: 0050107', '3DL1: 008', 
				'3DS1: 013011', '2DL5A: 0010102', '2DS5: 0020101', '2DS1: 0020101', '2DS4: 0030101', '3DL2: 01001', '3DL2: 018'	
			],
			"olga": [
				'3DL3: 00201', '3DL3: 00202' , '2DL3: 0010105', '2DL3: 0010105', '2DP1: 0020105', '2DP1: 006', '2DL1: 0030204', '2DL1: 0030204', '3DP1: 0030201', '3DP1: 0030201', '2DL4: 0050103', '2DL4: 00901', '3DL1: 0010102', 
				'3DL1: 0050101', '3DS1: 0130107', '2DL5A: 00103', '2DS5: 0020103', '2DS1: 0020101', '2DS4: 010', '3DL2: 0070101', '3DL2: 0070102'
			],
			"rsh": [
				'3DL3: 00202', '3DL3: 0040202', '2DS2: 0010108', '2DL2: 0030104', '2DL3: 0010107', '2DL5B: 004', '2DP1: 0020110', '2DP1: 009', '2DL1: 0030205', '2DL1: 01201', '3DP1: 0030401', '3DP1: 008', '2DL4: 0010307', 
				'2DL4: 00901', '3DL1: 0050101', '3DL1: 01701', '2DS5: 006', '2DS4: 0060102', '3DL2: 023', '3DL2: 056'
			],
			"wt51": [
				'3DL3: 0090101', '3DL3: 036', '2DS2: 0010103', '2DL2: 0010107', '2DL3: 006', '2DL5B: 0020103', '2DS3: 0020103', '2DS3: 0010302', '2DP1: 0010202', '2DP1: 004', '2DL1: 01201', '2DL1: 01201', '3DP1: 00303', '3DP1: 007', 
				'2DL4: 0050105', '2DL4: 0050103', '3DS1: 0130102', '3DS1: 0130102', '2DL5A: 0010103', '2DL5A: 0050104', '2DS5: 0020101', '2DS1: 0020103', '3DL2: 00202',  '3DL2: 00903'	
			]
}


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


"""
	Analysis step, make plot for each gen.
"""
def analysis_step(input_file_pyc: str, output_folder: str, alels_in_aligment: list):
	KIR_alels_dictionary = create_alels_translate_dictionary()

	
	with open(input_file_pyc, 'rb') as handle:
		alels_statistics = pickle.load(handle)


	# convert to right format
	alels_in_aligmnet_format = []
	for item in alels_in_aligment:
		item = item.replace(':', '*')
		gen, alel = item.split('*') 
		if('KIR' in item):
			alels_in_aligmnet_format.append(gen.strip()+"*"+alel.strip())
		else:
			alels_in_aligmnet_format.append("KIR"+gen.strip()+"*"+alel.strip())
	
	# lost have to be unique because there cen be some alels twice
	lost_set = set(alels_in_aligmnet_format.copy())
	lost = list(lost_set)

	alels_statistics_sort_by_gens = dict()
	for key, statistics in alels_statistics.items():
		trasnlate_alel = KIR_alels_dictionary[key] # translate
		gen = trasnlate_alel.split('*')[0] 
		
		if(gen not in alels_statistics_sort_by_gens):
			alels_statistics_sort_by_gens[gen] = dict()
			
		alels_statistics_sort_by_gens[gen][trasnlate_alel] = statistics

	moreover_gens = list()
	count_alels = 0
	for gen, gen_alel_statistic in alels_statistics_sort_by_gens.items():
		moreover_gens.append(gen)
		print(gen)
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
					statistics['cov_nucleotids'][i]=float('nan')
					#statistics['cov_nucleotids'][i]=0

		for alel, statistics in gen_alel_statistic.items():	
			count_alels += 1
			xf = np.linspace(0, len(statistics['cov_nucleotids']), len(statistics['cov_nucleotids']))

			if(len(gen_alel_statistic) == 1):
				axs.set_ylim(0, max_height+10)
				axs.set_xlim(0, max_lenght)

				# if alel is in haplotype which was aligment
				make_red = False
				for item in alels_in_aligmnet_format:
					# substring in string
					if(item in alel):
						make_red = True
						
						if(item in lost):
							lost.remove(item)

						if(gen in moreover_gens):
							moreover_gens.remove(gen)
						break

				if(make_red):
					axs.plot(xf, statistics['cov_nucleotids'], '-r')
					axs.fill_between(xf, 0, statistics['cov_nucleotids'], facecolor='red')
				else:	
					axs.plot(xf, statistics['cov_nucleotids'], '-b')
					axs.fill_between(xf, 0, statistics['cov_nucleotids'], facecolor='blue')
				
				# marker empty place of graph
				axs.axvspan(len(statistics['cov_nucleotids']), max_lenght, facecolor='grey')
				axs.set_aspect('auto', 'box')
				axs.text(max_lenght+10, max_height*0.5, alel.split('*')[1]+", cov: "+str(round(statistics['normalize_coverage'], 2))+"%, len: "+str(len(statistics['cov_nucleotids'])) , fontsize=9)
				plt.tight_layout()
			else:
				#fig.suptitle(alel)
				# set x "label" just for last alel
				if(alel != list(gen_alel_statistic.keys())[-1]):
					axs[alel_number].set_xticks([])

				axs[alel_number].set_ylim(0, max_height+10)
				axs[alel_number].set_xlim(0, max_lenght)

				# if alel is in haplotype which was aligment
				make_red = False
				for item in alels_in_aligmnet_format:
					# substring in string
					if(item in alel):
						make_red = True

						if(item in lost):
							lost.remove(item)
						if(gen in moreover_gens):
							moreover_gens.remove(gen)
						break

				if(make_red):
					axs[alel_number].plot(xf, statistics['cov_nucleotids'], '-r')
					axs[alel_number].fill_between(xf, 0, statistics['cov_nucleotids'], facecolor='red')
				else:	
					axs[alel_number].plot(xf, statistics['cov_nucleotids'], '-b')
					axs[alel_number].fill_between(xf, 0, statistics['cov_nucleotids'], facecolor='blue')
				
				# marker empty place of graph
				axs[alel_number].axvspan(len(statistics['cov_nucleotids']), max_lenght, facecolor='grey')
				axs[alel_number].set_aspect('auto', 'box')
				axs[alel_number].text(max_lenght+10, max_height*0.5, alel.split('*')[1]+", cov: "+str(round(statistics['normalize_coverage'], 2))+"%, len: "+str(len(statistics['cov_nucleotids'])) , fontsize=9)

				# so slow
				#axs[alel_number].bar(xf, statistics['cov_nucleotids'], color='#bd2a47', width=1, linewidth='0.1',  edgecolor='white', label="a")
				alel_number+=1

		#plt.tight_layout()
		if(len(gen_alel_statistic) == 1):
			axs.set_title(gen)
		else:
			fig.suptitle(gen)
		#fig.set_title(gen)

		plt.gcf().set_size_inches(max_lenght/1000, len(gen_alel_statistic))
		#plt.gcf().set_size_inches(20, len(alels_statistics_sort_by_gens))
		#plt.gcf().set_size_inches(len(alels_statistics_sort_by_gens), 9)
		#plt.savefig(os.path.join(PLOT_OUTPUT_FOLDER, gen+"compare_coverage.svg"), format='svg', bbox_inches='tight', aspect='auto', dpi=3000)
		plt.savefig(os.path.join(output_folder, gen+"compare_coverage.svg"), format='svg', bbox_inches='tight', dpi=3000)
		plt.clf() # clear plot
		plt.close()

	print("Lost: ", lost)
	print("Count alels: ", count_alels)

	return {"count_alels": count_alels, "lost": lost, "moreover_gens": moreover_gens}



"""
	Make result table to put into latex.
"""
def make_latech_table(genome: str, alels_in_aligment: list, result: dict):
	alels_unique = list(set(alels_in_aligment))
	output = genome+" & "+str(len(alels_in_aligment))+" ("+str(len(alels_in_aligment) - len(alels_unique))+")"

	more_over_step_before = []
	lost_step_before = []
	for step, data in result.items():
		if(step in STEPS_LATECH_TABLE):
			output += " & "+str(data['count_alels'])+" & "+str(len(data['lost']))+" & "

			new_lost = []
			for item in data['lost']:
				if(item not in lost_step_before):
					new_lost.append(item)
			
			lost_step_before = data['lost']

			if(len(new_lost) == 0):
				output+=" - "
			elif(len(new_lost) == 1):
				output+=str(new_lost[0].replace('KIR',''))

			else:
				output+="\\Gape[0pt][2pt]{\\makecell[l]{"

				for item in new_lost:
					if(item == new_lost[-1]):
						output+= item.replace('KIR','')
					else:	
						output+= item.replace('KIR','')+" \\\\ "
				output += "}}"


			output+=" & "+str(len(data['moreover_gens'])) + " & "

			new_more_over = []
			for item in data['moreover_gens']:
				if(item not in more_over_step_before):
					new_more_over.append(item)
			more_over_step_before = data['moreover_gens']

			if(len(new_more_over)  == 0):
				output += " - "
			elif(len(new_more_over) == 1):
				output += new_more_over[0].replace('KIR','')
			else:
				output+="\\Gape[0pt][2pt]{\\makecell[l]{"

				for item in new_more_over:
					if(item == new_more_over[-1]):
						output+= item.replace('KIR','')
					else:	
						output+= item.replace('KIR','')+" \\\\ "
				output += "}}"

			if(step == STEPS[-1]):
				true_positive = len(alels_in_aligment) - (len(alels_in_aligment) - len(alels_unique)) - len(data['lost'])
				false_positive = data['count_alels'] - true_positive
				false_negative = len(data['lost'])
				output += "& "+str(true_positive)+" & "+str(false_positive)+" & "+str(false_negative)

	output+= " \\\\ \n"
	return output



def run():
	result = dict()
	output = ""
	true_positive = 0
	false_positive = 0
	false_negative = 0

	for genome in GENOMES_LIST:
		print("analysis for genome: ", genome)
		alels_in_aligment = GENOMES_ALLELES[genome]
		
		for step in STEPS:
			file = os.path.join(ALELS_STATISTICS_PYC_FOLDER, genome+"_"+ALELS_STATISTICS_PYC_REFERENCE_NAME+"_"+ALELS_STATISTICS_PYC_EXPERIMENT+"_"+step+".pyc")
			output_folder = os.path.join(PLOT_OUTPUT_FOLDER,genome+"_"+step)
			os.mkdir(output_folder)
			result[step] = analysis_step(file, output_folder, alels_in_aligment)
		
		output += make_latech_table(genome, alels_in_aligment, result)

		# eval classifier
		alels_unique = list(set(alels_in_aligment))		
		true_positive_temp = len(alels_in_aligment) - (len(alels_in_aligment) - len(alels_unique)) - len(result[STEPS[-1]]['lost'])
		true_positive += true_positive_temp
		false_positive += result[STEPS[-1]]['count_alels'] - true_positive_temp
		false_negative += len(result[STEPS[-1]]['lost'])

	print(output)
	print("Eval classifier: ", "precision = ", true_positive/(true_positive+false_positive), "recall = ", true_positive/(true_positive+false_negative) )

run()
