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
ALELS_STATISTICS_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/data/statistics/test11_KIR_gen_exp2_step4.pyc"
PLOT_OUTPUT_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/analysis_after_align_result"

# ANALYSIS IN ALIGMENT
#AMALA 
#ALELS_IN_ALIGMENT = ['3DL3: 0040201', '3DL3:00802', '2DS2: 0010101', '2DL2: 0030102', '2DL3: 0010109', '2DP1: 0020108', '2DL1: 0030201', '3DP1: 007', '3DP1: 0090101', '2DL4: 0010201', '2DL4: 00501', '3DL1: 0150201', '3DS1: 0130101', '2DL5A: 00102', '2DS5: 0020101', '2DS1: 0020106', '2DS4: 0010101', '3DL2: 0020105', '3DL2:0070102']

#BOB 
#ALELS_IN_ALIGMENT = ['3DL3: 00101', '3DL3: 019', '2DS2: 0010104', '2DL2: 0030101', '2DL3: 0020102', '2DP1: 0030101', '2DL1: 0030210', '3DP1: 002', '3DP1: 0030203', '2DL4: 0010202', '2DL4: 00501', '3DL1: 002','3DS1: 0130105', '2DL5A: 0010101', '2DS5: 0020104', '2DS1: 0020101', '2DS4: 0010105', '3DL2: 0020101' , '3DL2:0070102']

#cox
#ALELS_IN_ALIGMENT = ['3DL3: 00102', '3DL3: 0090101', '2DL3: 0020101', '2DL3: 006', '2DP1: 0030102', '2DP1: 0030102', '2DL1: 0020102', '2DL1: 0020102', '3DP1: 005', '3DP1: 006', '2DL4: 00501', '2DL4: 00901', '3DL1: 0050103', '3DS1: 055', '2DS5: 0020102', '2DS1: 0020105', '2DS4: 010', '3DL2: 0010301', '3DL2: 0070103']

#ho301
#ALELS_IN_ALIGMENT = ['3DL3: 0140201', '3DL3: 0140201', '2DS2: 0010104', '2DS2: 0010106', '2DL2: 0010103', '2DL2: 0030107', '2DL5B: 010', '2DL5B: 010', '2DS3: 0010301', '2DS3: 0020103', '2DP1: 0010202', '2DP1: 0010202', '2DL1: 00402', '2DL1: 010', '3DP1: 0030101', '3DP1: 004', '2DL4: 0010201', '2DL4: 0010201', '3DL1: 002', '3DL1: 002', '2DS4: 0010109', '2DS4: 0010109', '3DL2: 0020102', '3DL2: 0020106']

#jvm
#ALELS_IN_ALIGMENT = ['3DL3: 00801', '3DL3: 0140201', '2DS2: 0010110', '2DL2: 0030102', '2DL3: 010', '2DP1: 004', '2DL1: 0030203', '3DP1: 001', '3DP1: 0030202', '2DL4: 0010304', '2DL4: 0080101', '3DL1: 0010104', '3DL1: 008', '2DS4: 0030103','2DS4: 0030103', '3DL2: 0010101', '3DL2:018']

#kas011
#ALELS_IN_ALIGMENT = ['3DL3: 0090101', '3DL3: 0140203', '2DL3: 0020103', '2DL3: 0020103', '2DP1: 0020104', '2DP1: 0030101', '2DL1: 0020101', '2DL1: 0030209', '3DP1: 0030206', '3DP1: 009', '2DL4: 0010301', '2DL4: 00501', '3DL1: 008', '3DS1: 013011', '2DL5A: 0010102', '2DS5: 0020101', '2DS1: 0020101', '2DS4: 0030101', '3DL2: 01001', '3DL2: 018'	]

#olga
#ALELS_IN_ALIGMENT = ['3DL3: 00201', '3DL3: 00202' , '2DL3: 0010105', '2DL3: 0010105', '2DP1: 0020105', '2DP1: 006', '2DL1: 0030204', '2DL1: 0030204', '3DP1: 0030201', '3DP1: 0030201', '2DL4: 00501', '2DL4: 00901', '3DL1: 0010102', '3DL1: 0050101', '3DS1: 0130107', '2DL5A: 00103', '2DS5: 0020103', '2DS1: 0020101', '2DS4: 010', '3DL2: 0070101', '3DL2: 0070102']

#rsh
#ALELS_IN_ALIGMENT = ['3DL3: 00202', '3DL3: 0040202', '2DS2: 0010108', '2DL2: 0030104', '2DL3: 0010107', '2DL5B: 004', '2DP1: 0020110', '2DP1: 009', '2DL1: 0030205', '2DL1: 01201', '3DP1: 0030401', '3DP1: 008', '2DL4: 0010307', '2DL4: 00901', '3DL1: 0050101', '3DL1: 01701', '2DS5: 006', '2DS4: 0060102', '3DL2: 023', '3DL2: 056']

#test1
#ALELS_IN_ALIGMENT = ['3DL3: 0030101', '3DL3: 0140201', '2DS2: 0010104', '2DL2: 0030105', '2DL3: 0020101', '2DL5B: 01301', '2DS3: 0020102', '2DS3: 0020102', '2DP1: 0030101' , '2DP1:0010203', '2DL1: 0030203', '2DL1:007', '3DP1: 004', '3DP1: 004', '2DL4: 0080104', '2DL4: 010', '3DL1: 0150101', '3DS1: 014', '2DL5A: 00102', '2DS1: 0020104', '2DS4: 0010103', '3DL2: 00501', '3DL2: 00501']		

#test2
#ALELS_IN_ALIGMENT = ['3DL3: 0090102', '3DL3: 019', '2DL2: 0010101', '2DL3: 0010111', '2DL5B: 0020101', '2DP1: 0020107', '2DP1: 0030102', '2DL1: 0020102', '2DL1: 0030210', '3DP1: 004', '3DP1: 01001', '2DL4: 0080104', '2DL4: 0080104','3DL1: 0070101', '3DS1: 078', '2DL5A: 0050101', '2DS5: 010', '2DS1: 0020102', '2DS4:0010103', '3DL2: 0020101', '3DL2: 00501']

#test3 
#ALELS_IN_ALIGMENT = ['3DL3: 005', '3DL3: 0140201', '2DL3: 0010101', '2DL3: 0020103', '2DP1: 0020108', '2DP1: 0020108', '2DL1: 0040101', '2DL1: 008', '3DP1: 0030102', '3DP1: 00902', '2DL4: 0010306', '2DL4: 00501', '3DL1: 002', '3DL1: 0040101', '3DL2: 0010301', '3DL2: 008']

#test4
#ALELS_IN_ALIGMENT = ['3DL3: 0030104', '3DL3: 007', '2DS2: 0010105', '2DL2: 0030101', '2DL3: 0010102', '2DP1: 008', '2DL1: 007', '3DP1: 007', '3DP1: 00902', '2DL4: 0010307', '2DL4: 0080104', '3DL1: 0150202', '3DS1: 055', '2DS5: 007', '2DS1: 0020101', '2DS4: 0060101', '3DL2: 0020101', '3DL2: 00903']

#test5
#ALELS_IN_ALIGMENT = ['3DL3: 0140202', '3DL3: 036', '2DL3: 0010109', '2DL3: 006', '2DS3: 0010301', '2DP1: 0030102', '2DP1: 009', '2DL1: 0030208', '2DL1: 00303', '3DP1: 001', '3DP1: 002', '2DL4: 0010202', '2DL4: 0010202', '3DL1: 0200101', '3DS1: 0130102', '2DL5A: 0010102', '2DS1: 0020105', '3DL2: 00202', '3DL2: 018']

#test6
#ALELS_IN_ALIGMENT = ['3DL3: 0090102', '3DL3: 0140203', '2DS2: 0010111', '2DL2: 0010105', '2DL3: 0010102', '2DL5B: 0080101', '2DS3: 0010302', '2DP1: 0020103', '2DP1: 010', '2DL1: 0030203', '2DL1: 0040102', '3DP1: 0030202', '3DP1: 0030402','2DL4: 0010303', '2DL4: 00901', '3DL1: 0050102', '3DL1: 0250102', '2DS4: 0010104', '2DS4: 010', '3DL2: 0010302', '3DL2: 01001']

#test7
#ALELS_IN_ALIGMENT = ['3DL3: 00802', '3DL3: 0090103', '2DL3: 0010103', '2DL3: 0010108', '2DP1: 0020106', '2DP1: 004', '2DL1: 0030204', '2DL1: 0030205', '3DP1: 0030202', '3DP1: 0030202', '2DL4: 0010201', '2DL4: 0010305', '3DL1: 008','3DL1: 0150203', '2DS4: 0010107', '2DS4: 0030104', '3DL2: 0020105', '3DL2: 00901']

#test8
#ALELS_IN_ALIGMENT = ['3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DL5B: 0070101', '2DS3: 0020101', '2DS3: 0020101', '2DS3: 0010302', '2DP1: 0030102', '2DL1: 00402', '3DP1: 0030101','3DP1: 005', '2DL4: 00104', '2DL4: 0080104', '3DS1: 0130104', '3DS1: 055', '2DL5A: 0050102', '2DL5A: 0050102', '2DS1: 0020102', '2DS1: 0020105', '3DL2: 0010102', '3DL2: 0070102']

#test9
#ALELS_IN_ALIGMENT = ['3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DL5B: 0070101', '2DS3: 0020101', '2DS3: 0020101', '2DP1: 0030102', '2DL1: 00402', '3DP1: 0030101', '3DP1: 005','2DL4: 00104', '2DL4: 0080104', '3DL1: 0150208', '3DS1: 0130104', '2DL5A: 01201', '2DS1: 0020102', '2DS4: 0040101', '3DL2: 0010102', '3DL2: 0070102']

#test10
#ALELS_IN_ALIGMENT = ['3DL3: 0030101', '3DL3: 0140201', '2DS2: 0010104', '2DL2: 0030105', '2DL3: 0020101', '2DL5B: 01301', '2DS3: 0020102', '2DS3: 0020102', '2DP1: 0010203', '2DP1: 0030101', '2DL1: 0030203', '2DL1: 007', '3DP1: 004','3DP1: 004', '2DL4: 0080104', '2DL4: 010', '3DL1: 0150101', '3DS1: 014', '2DS1: 0020104', '3DL2: 00501', '3DL2: 00501']

#test11
ALELS_IN_ALIGMENT = ['3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DP1: 0030102', '2DP1: 008', '2DL1: 00402', '2DL1: 00402', '3DP1: 0030101', '3DP1: 005', '2DL4: 00104', '2DL4: 0080104','3DS1: 0130104', '3DS1: 055', '2DL5A: 0050102', '2DL5A: 0050102', '2DS5: 0020102', '2DS5: 0020103', '2DS1: 0020102', '2DS1: 0020105', '3DL2: 0010102', '3DL2: 0070102']


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

	print(ALELS_STATISTICS_FILE_PYC)
	
	with open(ALELS_STATISTICS_FILE_PYC, 'rb') as handle:
		alels_statistics = pickle.load(handle)

	# convert to right format
	alels_in_aligmnet_format = []
	for item in ALELS_IN_ALIGMENT:
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

	count_alels = 0
	for gen, gen_alel_statistic in alels_statistics_sort_by_gens.items():
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
		plt.savefig(os.path.join(PLOT_OUTPUT_FOLDER, gen+"compare_coverage.svg"), format='svg', bbox_inches='tight', dpi=3000)
		plt.clf() # clear plot

	print("Lost: ", lost)
	print("Count alels: ", count_alels)

	# musim vzit kazde dvě alelely a zase si urcit vzdalenost mezi nimi 
	# a od nejake vzdalenosti se mrknout jak jsou zarovanany tam kde je ten rozdil mezi nimi

run()
