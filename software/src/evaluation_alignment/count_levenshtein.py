import pysam
import sys
import os
# regular expression for removing white space
import re
#save structure into file
import pickle
import json

import Levenshtein

import config

END_REFERENCE_FILE = ".fa"

"""
	half is enought - mirror matrix by diagonal
"""

def run():
	# get reference file
	ref_gen_file = open(config.REFERENCE_KIR_GENS_PSEUDOGENS_FILE, 'r')
	ref_gens = ref_gen_file.read()
	ref_gen_file.close()
	
	alels = ref_gens.split(">")[1:]

	alels_dict = dict()

	for alel in alels:
		head_alel, body_alel = alel.split('\n', 1)
		
		key = head_alel.split(' ', 1)[0]
		#print(key)
		#get away all white space
		body_alel=re.sub(r"\s+", "", body_alel, flags=re.UNICODE)
		alels_dict[key] = body_alel


	distance_txt_file = open(config.ALELS_DISTANCE_FILE_TXT, "a+")

	alels_distance = dict()
	current_key = ""
	for key1, alel1 in alels_dict.items():
		print("counting distance for: ", key1)
		for key2, alel2 in alels_dict.items():
			
			if(current_key != key1):
				alels_distance[key1] = dict()
				current_key = key1

			if(key1 == key2): #half is enought - mirror matrix by diagonal
				break

			alels_distance[key1][key2] = dict()
			alels_distance[key1][key2]['distance'] = Levenshtein.distance(alel1, alel2)
			alels_distance[key1][key2]['changes'] = Levenshtein.editops(alel1, alel2)

		print(key1+" - "+ json.dumps(alels_distance[key1]), file=distance_txt_file)
		

	with open(config.ALELS_DISTANCE_FILE_PYC, 'wb') as handle:
   		pickle.dump(alels_distance, handle, protocol=pickle.HIGHEST_PROTOCOL)

   	# podobnost alel danyho genu
   	# udelat na to pak nejakej graf jak jsou 
   	# pro jednu alelu mozna udělat graf kde budou puntiky s tim že tam bude počet chyb a bude tam jestli to je stejnej gen nebo jinej gen