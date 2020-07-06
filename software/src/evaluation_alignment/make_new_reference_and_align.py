import pysam
import sys
import os
import pickle

import config

import subprocess
""" Because there can be more files, like nuc and another"""
END_REFERENCE_GEN_FILE = "_gen.fasta"


"""
	Prepare dictionary for create new dictionary.
"""
def get_all_references_alels():
	# make dictionary KIR:KIR00978 => KIR2DL1*0010102
	print("Preparing alels from files ", config.REFERENCE_KIR_GENS)
	KIR_alels_dictionary = dict()
	file_content = ""

	with open(config.REFERENCE_KIR_GENS) as fileHandle:
		file_content = fileHandle.read()

	alels = file_content.split('>')

	for alel in alels[1:]: # alels[0] is empty
		# KIR:KIR00037 KIR2DS2*0010101 14577 bp
		alel_head, alel_body = alel.split('\n', 1) 
		alel_marker = alel_head.split()[0]
		alel_number = alel_marker.split(":")[1]

		alel_number = alel_number.replace('KIR', '')

		# have to put back >
		KIR_alels_dictionary[alel_number] = '>'+alel

	return KIR_alels_dictionary


"""
	Make a KIR dictionary. 
	Then find KIR gen in dictionary
	try to find which alels start with wanted alels number
"""
def create_new_reference(haplotype: list, haplotype_name: str, all_references_alels: dict, file_sufix: str):
	result_haplotype = ""
	result_haplotype_legend = ""

	for item in haplotype:	
		found = False
		# KIR:KIR00432
		wanted_alel = item.split(":")[1]

		wanted_alel = wanted_alel.replace('KIR', '')
		wanted_alel = wanted_alel.strip()

		if(wanted_alel in all_references_alels):
			result_haplotype += all_references_alels[wanted_alel]
		else:
			print("Can not find item ", wanted_gen+"*"+wanted_alel)
							

	output_file = os.path.join("/home/kate/Dokumenty/FAV/Diplomka/software/data/reference/", haplotype_name+"_"+file_sufix+".fasta")
	haplotype_result_file_rename = open(output_file, "w+")
	haplotype_result_file_rename.write(result_haplotype)
	haplotype_result_file_rename.close()

	#print("Haplotype legend ", haplotype_name, ":", result_haplotype_legend)
	print("Create reference file: ", output_file)
	return output_file


"""
	Align by bowtie.
	First have to make bowtie index, then align.
	cwd - executable path - where will be generate bowtie index, and then have to run bowtie align
"""
def align(executable_path: str, new_reference, basic_name: str, file_sufix: str):
	# create bowtie index
	namefile_index = basic_name+file_sufix
	process = subprocess.run([config.BOWTIE_HOME_DIRECTORY+"/bowtie2-build", "--threads", str(config.BOWTIE_THREADS), new_reference, namefile_index], cwd=executable_path)

	# find reads
	
	reads_files = [f for f in os.listdir(config.READS_FOLDER) if os.path.isfile(os.path.join(config.READS_FOLDER, f)) and os.path.splitext(f)[1] == '.fq']
	print("Find reads: ", reads_files)

	find = False
	haplotype=""
	for read_file in reads_files:
		haplotype = read_file.split('.')[0]
		haplotype = haplotype[:-1]

		if basic_name.startswith(haplotype):
			find = True
			break

	if(find == False):
		print("Reads for", basic_name, "not found.")
		return

	read1= os.path.join(config.READS_FOLDER, haplotype+"1.fq")
	read2= os.path.join(config.READS_FOLDER, haplotype+"2.fq")
	
	print("read1: ", read1)
	print("read2: ", read2)
	# run bowtie
	namefile_align = os.path.join("/home/kate/Dokumenty/FAV/Diplomka/software/data/temp", basic_name+file_sufix+".sam")
	process = subprocess.run([config.BOWTIE_HOME_DIRECTORY+"/bowtie2", "--threads", str(config.BOWTIE_THREADS), "-x", namefile_index, "-1", read1,"-2", read2, "-S", namefile_align], cwd=executable_path) 	
	print("create file", namefile_align)

	return namefile_align