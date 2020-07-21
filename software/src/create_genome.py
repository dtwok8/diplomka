import os
from datetime import datetime

import config

"""
	Get all alels in config.REFERENCE_KIR_GENS_FILE and make one dictionary 
	KIR_alels_dictionary[2DS2][0010101]: >KIR:KIR00037 KIR2DS2*0010101 14577 bp ATCGGGAT...
	Then all genomes in config and make new genome file from KIR_alels_dictionary
"""

""" Because there can be more files, like nuc and another"""
END_REFERENCE_GEN_FILE = "_gen.fasta"
END_READ_FILE = ".fa"

"""
	Prepare alels into dictionary. 
	KIR_alels_dictionary[2DS2][0010101]: >KIR:KIR00037 KIR2DS2*0010101 14577 bp ATCGGGAT...
"""

def create_alels_dictionary():
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
		alel_marker = alel_head.split()[1]
		gen, alel_number = alel_marker.split("*")

		gen = gen.replace('KIR', '')
		if(gen not in KIR_alels_dictionary):
			KIR_alels_dictionary[gen] = dict()
		
		# need cut KIR
		# have to put back >
		KIR_alels_dictionary[gen][alel_number] = '>'+alel

	return KIR_alels_dictionary


"""
	Make a KIR dictionary. 
	Then find KIR gen in dictionary
	try to find which alels start with wanted alels number
"""
def run():
	KIR_alels_dictionary = create_alels_dictionary()

	result_genome = ""
	result_genome_legend = ""

	for genome_name, genome in config.GENOMES.items():
		print("Creating genome: ", genome_name, " ...")
		result_genome_legend = ""
		result_genome = ""
		
		for item in genome:	
			found = False
			# can be 3DL3:00402 or 2DS2*00101
			wanted_gen, wanted_alel = item.replace(":", "*").split('*')

			# can be KIR3DL3 or 3DL3 
			wanted_gen = wanted_gen.strip() # 3DL3
			wanted_alel = wanted_alel.strip() #004002

			if(wanted_gen.startswith('KIR')):
				wanted_gen = wanted_gen.replace('KIR', '')

			if(wanted_gen in KIR_alels_dictionary):
				
				# alel does not have to be full, can be just the start in genome
				for gen_alel_key, gen_alel_body in KIR_alels_dictionary[wanted_gen].items():
					if(gen_alel_key.startswith(wanted_alel)):
						result_genome += gen_alel_body
						result_genome_legend = result_genome_legend + wanted_gen+"*"+ gen_alel_key+ ", "
						found = True
						break;
				
				if(found == False):
					print("Can not find item ", wanted_gen+"*"+wanted_alel)
							


		output_file = os.path.join(config.GENOME_FOLDER, genome_name+END_READ_FILE )
		genome_result_file_rename = open(output_file, "w+")
		genome_result_file_rename.write(result_genome)
		genome_result_file_rename.close()

		print("genome legend ", genome_name, ":", result_genome_legend)
		print("Create file: ", output_file)
