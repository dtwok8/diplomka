import pickle
import sys
import os
import statistics
import re #regular expresion

import matplotlib.pyplot as plt
import numpy as np
import Levenshtein

NUC_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_nuc.fasta"
GEN_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"

DISTANCE_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/alels_distance_new.pyc"
PLOT_OUTPUT_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/analysis_alels_result"



def compare_gen_and_nuc_first():
	not_found = 0
	something_wrong = 0
	match = 0
	count_nuc = 0
	count_gen = 0

	KIR_alels_dictionary = dict()
	
	with open(GEN_FILE, "r") as openfileobject:
		for line in openfileobject:			
			if(line.startswith('>KIR')):
				count_gen+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>KIR'
				KIR_alels_dictionary[split_line[0][5:]] = split_line[1]


	with open(NUC_FILE, "r") as openfileobject:
		for line in openfileobject:
			# IPD:KIR00874 KIR3DS1*108
			if(line.startswith('>IPD')):
				count_nuc+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>'
				key = split_line[0][5:]
				value = split_line[1]

				if(key in KIR_alels_dictionary):
					if(KIR_alels_dictionary[key] != value):
						#print("something wrong ", key, " : ", KIR_alels_dictionary[key], "!=", value )
						something_wrong+=1
					else:
						#print("match ", key)
						match+=1
				else:
					#print("key not found: ", key)
					not_found+=1

	print("count_gen: ", count_gen, ", count_nuc: ", count_nuc, ", match: ", match, ", key not found: ", not_found, ", something wrong: ", something_wrong)


def compare_gen_and_nuc_second():
	not_found = 0
	something_wrong = 0
	match = 0
	count_nuc = 0
	count_gen = 0

	KIR_alels_dictionary = dict()
	
	with open(GEN_FILE, "r") as openfileobject:
		for line in openfileobject:			
			if(line.startswith('>KIR')):
				count_gen+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>KIR'
				KIR_alels_dictionary[split_line[1]] = split_line[0][5:]


	with open(NUC_FILE, "r") as openfileobject:
		for line in openfileobject:
			# IPD:KIR00874 KIR3DS1*108
			if(line.startswith('>IPD')):
				count_nuc+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>IPD:'
				value = split_line[0][5:]
				key = split_line[1]

				if(key in KIR_alels_dictionary):
					if(KIR_alels_dictionary[key] != value):
						#print("something wrong ", key, " : ", KIR_alels_dictionary[key], "!=", value )
						something_wrong+=1
					else:
						#print("match ", key)
						match+=1
				else:
					#print("key not found: ", key)
					not_found+=1

	print("count_alel: ", count_gen, ", count_nuc: ", count_nuc, ", match: ", match, ", key not found: ", not_found, ", something wrong: ", something_wrong)


def get_min_max_distance():
	with open(DISTANCE_FILE_PYC, 'rb') as handle:
		alels_distance = pickle.load(handle)

	max_distance = 0
	min_distance = sys.maxsize 
	sum_distance = 0
	count_distance = 0
	list_distance = list()
	for alel1, dic_distance in alels_distance.items():
		for alel2, distance in dic_distance.items():
			if(alel1 != alel2):
				if(max_distance < distance):
					max_distance = distance
				if(min_distance > distance):
					min_distance = distance

				sum_distance += distance
				count_distance += 1
				list_distance.append(distance)

	print("max: ", max_distance, ", min: ", min_distance, ", average_distance: ", sum_distance/count_distance, ", median: ", statistics.median(list_distance)) 


def levenhstein_analyze():
	print("Levenhstein analyze")
	print(Levenshtein.distance('SPAM', 'PARK'))
	print(Levenshtein.editops('SPAM', 'PARK'))

	print(Levenshtein.distance('PARK', 'SPAM'))
	print(Levenshtein.editops('PARK', 'SPAM'))	
	print("----------------------------------")


def count_levenhstein_distance():
	print("counting distance ... ")
	# get reference file
	ref_gen_file = open(GEN_FILE, 'r')
	ref_gens = ref_gen_file.read()
	ref_gen_file.close()
	
	alels = ref_gens.split(">")[1:]

	alels_dict = dict()

	for alel in alels:
		head_alel, body_alel = alel.split('\n', 1)
		
		key = head_alel.split(' ', 1)[0]
		#get away all white space
		body_alel=re.sub(r"\s+", "", body_alel, flags=re.UNICODE)
		alels_dict[key] = body_alel


	alels_distance = dict()
	current_key = ""
	for key1, alel1 in alels_dict.items():
		print("counting distance for: ", key1)
		for key2, alel2 in alels_dict.items():
			
			if(current_key != key1):
				alels_distance[key1] = dict()
				current_key = key1


			alels_distance[key1][key2] = dict()
			alels_distance[key1][key2]['distance'] = Levenshtein.distance(alel1, alel2)
		

	with open(DISTANCE_FILE_PYC, 'wb') as handle:
   		pickle.dump(alels_distance, handle, protocol=pickle.HIGHEST_PROTOCOL)


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
	Translate from  KIR:KIR00978 => KIR2DL1*0010102
"""
def translate(alels_distance, alels_translate_dictionary):
	translate_alels_distance = dict()
	
	for key1, distance_dict in alels_distance.items():
		translate_key1 = alels_translate_dictionary[key1]
		translate_alels_distance[translate_key1] = dict()

		for key2, distance in distance_dict.items():
			translate_key2 = alels_translate_dictionary[key2]
			translate_alels_distance[translate_key1][translate_key2] = distance	

	return translate_alels_distance

"""
	In basic choose one alel from each gen. It can be easily change which alel analyse.
"""
def get_alels_for_distance_analyse(alels_distance):
	alels_for_distance_analyse = list()
	gen_complete = list()

	for alel in alels_distance:
		gen = alel.split('*')[0]
		
		if(gen not in gen_complete):
			alels_for_distance_analyse.append(alel)
			gen_complete.append(gen)

	return alels_for_distance_analyse

"""
	Make a plot for one alel from each gen to compare distance betweens alels.
"""
def analyze_distance_gene(alels_translate_dictionary):
	with open(DISTANCE_FILE_PYC, 'rb') as handle:
		alels_distance = pickle.load(handle)	

	# everything translate because need to recognize if it same gene or not
	translate_alels_distance = translate(alels_distance, alels_translate_dictionary)
	alels_distance = translate(alels_distance, alels_translate_dictionary)
	#find first alel for category, because you dont want plot all alels (460) I think that one alel from each gen is enought
	alels_for_distance_analyse = get_alels_for_distance_analyse(alels_distance) 

	gens_category = list()
	bars_by_gens_category = dict()

	# sort alels by gen to bars_by_gens_category
	# create  KIR2DL1*0010101[KIRD2L1] : list all distance between this alel and alles current gen
	# KIR2DL1*0010101[KIRD2L1] : [8, 10, 50, .. ]
	for alel1, distance_dict in translate_alels_distance.items():
		gen1 = alel1.split('*')[0]

		if(alel1 in alels_for_distance_analyse):
			if(gen1 not in gens_category):
				gens_category.append(gen1)

			if(alel1 not in bars_by_gens_category):
				bars_by_gens_category[alel1] = dict()

			for alel2, distance in distance_dict.items(): 
				gen2 = alel2.split('*')[0]

				if(gen2 not in bars_by_gens_category[alel1]):
					bars_by_gens_category[alel1][gen2] = dict()

				bars_by_gens_category[alel1][gen2][alel2] = distance	
	
	plot_one_alels_from_gen = dict()
	for alel, values in bars_by_gens_category.items():
		gen = alel.split('*')[0]

		if(gen not in plot_one_alels_from_gen):
			plot_one_alels_from_gen[gen] = alel
	
	#  width of bar
	barWidth = 0.5

	# draw bar plot for alel from bars_by_gens_category
	for ploting_alel, all_distance in bars_by_gens_category.items():
		last_positions = 0
		x_tick_position = list()

		for ploting_gen, gens_distances in all_distance.items():
			curent_gen_distance_plot = list()
			for key, distance in gens_distances.items():
				curent_gen_distance_plot.append(distance)

			x_tick_po = last_positions + (len(curent_gen_distance_plot) / 2)*barWidth
			x_tick_position.append(x_tick_po)
			current_positions = list()	

			#create positions for bars for alels distance current gen
			for i in range(len(curent_gen_distance_plot)):
				current_positions.append(last_positions)
				last_positions+=barWidth

			if(ploting_gen == ploting_alel.split('*')[0]): # make same gen red, other blue
				plt.bar(current_positions, curent_gen_distance_plot, color='#bd2a47', width=barWidth, linewidth='0.1',  edgecolor='white', label=gens_category)
			else:
				plt.bar(current_positions, curent_gen_distance_plot, color='#0d82ff', width=barWidth, linewidth='0.1',  edgecolor='white', label=gens_category)

			last_positions+=(5*barWidth) # space betweens gens	
		# Add xticks on the middle of the group bars
		plt.xlabel('Gen', fontweight='bold')
		plt.ylabel('Vzdálenost', fontweight='bold')
		plt.xticks(x_tick_position, gens_category, rotation='vertical')
		 
		# Create legend & Show graphic
		plt.tight_layout()
		plt.autoscale()

		print("ploting ", ploting_alel)
		plt.title(ploting_alel)
		plt.ylim(0, 16000)
		plt.gcf().set_size_inches(13, 9)
		plt.savefig(os.path.join(PLOT_OUTPUT_FOLDER, ploting_alel+"cross_gens.svg"), format='svg', bbox_inches='tight', aspect='auto', dpi=3000)
		plt.clf() # clear plot


		# ploting just one gen -----------------------------------------------------------------------
		alels_for_x_axis = list()
		last_positions = 0
		x_tick_position = list()
		ploting_gen = ploting_alel.split('*')[0]
		gens_distances = all_distance[ploting_gen]
		
		curent_gen_distance_plot = list()
		for key, distance in gens_distances.items():
			curent_gen_distance_plot.append(distance)
			alels_for_x_axis.append(key)

		current_positions = list()	

		#create positions for bars for alels distance current gen
		for i in range(len(curent_gen_distance_plot)):
			current_positions.append(last_positions)
			x_tick_position.append(last_positions)
			last_positions+=barWidth

		plt.bar(current_positions, curent_gen_distance_plot, color='#bd2a47', width=barWidth, linewidth='0.1',  edgecolor='white', label=alels_for_x_axis)

		last_positions+=(5*barWidth) # space betweens gens	
		# Add xticks on the middle of the group bars
		plt.xlabel('Gen', fontweight='bold')
		plt.ylabel('Vzdálenost', fontweight='bold')
		plt.xticks(x_tick_position, alels_for_x_axis, rotation='vertical')
		 
		# Create legend & Show graphic
		plt.tight_layout()
		plt.autoscale()
		print("ploting ", ploting_alel)
		plt.title(ploting_alel)
		plt.gcf().set_size_inches(13, 9)
		#plt.gcf().set_size_inches(13, 5)
		#plt.ylim(0, 100)
		plt.savefig(os.path.join(PLOT_OUTPUT_FOLDER, ploting_alel+"gen.svg"), format='svg', bbox_inches='tight', aspect='auto', dpi=3000)
		plt.clf() # clear plot


def main():
	compare_gen_and_nuc_first()
	count_levenhstein_distance()

	#compare_gen_and_nuc_second()
	#exit(1)
	# run just if you realy need, it takes long time more then 12 hours
	#
	#levenhstein_analyze()

	get_min_max_distance()
	alels_translate_dictionary = create_alels_translate_dictionary() 
	analyze_distance_gene(alels_translate_dictionary)

	

	


if __name__ == "__main__":
	main()