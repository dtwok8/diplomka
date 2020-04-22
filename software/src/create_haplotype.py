import os
from datetime import datetime

import config

"""
	Split content of file into list by alels and
	try to find the right alel
	compare KIR2DS2 and 0010101, separetly because blank space around it
	and compare KIR2DS3 because pseudogens are usually in another file (KIR_gen.fasta)
"""
def search_alel_in_file(filename: str, gen: list , result_haplotype, result_haplotype_legend):
	found_alel = False
	gen_file = open(filename, "r") 
	gen_file_str = gen_file.read()
	gen_file.close()

	alels = gen_file_str.split(">")
	
	for alel in alels[1:]: # alels[0] is empty
		# KIR:KIR00037 KIR2DS2*0010101 14577 bp
		head_alel = alel.split('\n')[0]
		alel_marker = head_alel.split()[1]
		compare_alel = alel_marker.split("*")

		if(compare_alel[0].strip() == gen[0].strip() and compare_alel[1].strip().startswith(gen[1].strip())):
			found_alel = True
			result_haplotype_legend = result_haplotype_legend + gen[0] + "*" + compare_alel[1] + ", "	
			# because remove > when split so need to add into 	
			result_haplotype = result_haplotype +">"+alel
			break

	if(found_alel == False):
		print("Alel not found: ", gen[0]+"*"+gen[1])

	return result_haplotype, result_haplotype_legend	


"""
	Try to find right file by 3DL3, 
	but pseudogens dont have specific file by this prefix
	so if it dont find file by prefix, try to search in REFERENCE_KIR_GENS_PSEUDOGENS_FILE
"""
def run():
	result_haplotype = ""
	result_haplotype_legend = ""

	for key, haplotype in config.HAPLOTYPES.items():
		print("Creating haplotype: ", key, " ...")
		for item in haplotype:	
			# can be 3DL3:00402 or 2DS2*00101
			gen = item.replace(":", "*").split('*')

			# can be KIR3DL3 or 3DL3
			gen[0] = gen[0].strip()
			gen[1] = gen[1].strip()

			if( not gen[0].startswith('KIR')):
				gen[0] = 'KIR'+gen[0]

			gen_file_name = os.path.join(config.REFERENCE_KIR_GENS_FOLDER, gen[0]+"_gen.fasta")
			if(os.path.isfile(gen_file_name)):
				result_haplotype, result_haplotype_legend = search_alel_in_file(gen_file_name, gen,  result_haplotype, result_haplotype_legend)	

			else:
				# pseudogen usually in file KIR_gen.fasta
				gen_file_name2 = config.REFERENCE_KIR_GENS_PSEUDOGENS_FILE

				if(os.path.isfile(gen_file_name2)):
					result_haplotype, result_haplotype_legend = search_alel_in_file(gen_file_name2, gen, result_haplotype, result_haplotype_legend)	
					
					if(gen[0] not in result_haplotype_legend):
						print("Can not find", item)
				else:
					print("Can not find file ", gen_file_name, "or", gen_file_name2)


		output_file = os.path.join(config.HAPLOTYPE_OUTPUT_FOLDER, key+".fa" )
		aligment_result_file_rename = open(output_file, "w+")
		aligment_result_file_rename.write(result_haplotype)
		aligment_result_file_rename.close()

		print("Haplotype legend: ", result_haplotype_legend)
		print("Create file: ", output_file)
