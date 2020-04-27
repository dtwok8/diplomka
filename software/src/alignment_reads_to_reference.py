import subprocess
import os

import config

""" Because there can be more files, like nuc and another"""
END_REFERENCE_GEN_FILE = "_gen.fasta"
"""
	Get all references gens and make index.
	Prepare for aligmnets.
	bowtie-build outputs a set of 6 files with suffixes .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt. 

"""
def bowtie_build_index():
	ref_kir_gen_files = [f for f in os.listdir(config.REFERENCE_KIR_GENS_FOLDER) if os.path.isfile(os.path.join(config.REFERENCE_KIR_GENS_FOLDER, f)) and f.endswith(END_REFERENCE_GEN_FILE) and (f!=".gitignore")]
	
	for i in range(0,len(ref_kir_gen_files)):
		name = os.path.splitext(os.path.basename(ref_kir_gen_files[i]))[0]
		# cw= index_foler - change folder where command will be execute, because there wasnt any parameter for output folder
		process = subprocess.run([config.BOWTIE_HOME_DIRECTORY+"/bowtie2-build", "--threads", str(config.BOWTIE_THREADS), config.REFERENCE_KIR_GENS_FOLDER+"/"+ref_kir_gen_files[i], name], cwd=config.BOWTIE_INDEX_FOLDER)
		print("aligment build index", name)


"""
	Compare all KIR gens with reads.
	Because Bowtie create 6 files for one reference gens, then we have to find unique files and the correct name for index.
"""
def alignment_read(read1: str, read2: str, basic_read_name: str):
	b_index = [f for f in os.listdir(config.BOWTIE_INDEX_FOLDER) if os.path.isfile(os.path.join(config.BOWTIE_INDEX_FOLDER, f)) and (f!=".gitignore")]

	if(len(b_index)==0):
		print("Error no index")

	b_index_basic_list = list()
	for index_file in b_index:
		b_index_basic = index_file.split('.')[0]

		b_index_basic_list.append(b_index_basic)	

	unique_set = set(b_index_basic_list)
	unique_b_index_list = list(unique_set)

	for i in range(0,len(unique_b_index_list)):
		index_name = os.path.basename(unique_b_index_list[i]).split('.')[0]
		print("align: ", basic_read_name," ", index_name)
		process = subprocess.run([config.BOWTIE_HOME_DIRECTORY+"/bowtie2", "--threads", str(config.BOWTIE_THREADS), "-x", index_name, "-1", read1,"-2", read2, "-S", os.path.join(config.ALIGNMENT_FOLDER, basic_read_name+"_"+index_name+"_align.sam")], cwd=config.BOWTIE_INDEX_FOLDER) 	



"""
	Run make index
	For all reads in folder reads run align for references gens
		- reads names has to be parser
"""
def run():
	if(config.BOWTIE_BUILD_INDEX):
		bowtie_build_index()
	
	all_reads = [f for f in os.listdir(config.READS_FOLDER) if os.path.isfile(os.path.join(config.READS_FOLDER, f)) and f.endswith(".fq") and (not f.endswith(".gitignore"))]

	print(all_reads)
	# we know that it fill end 1.fq a 2.fq
	# get basic name, unique withnout 1.fq and 2.fq
	read_basic_list = list()
	for read in all_reads:
		read_basic_name = read[:-4]

		read_basic_list.append(read_basic_name)	

	unique_set = set(read_basic_list)
	unique_read_list = list(unique_set)
	print("reads list: ", unique_read_list)

	for basic_read_name2 in unique_read_list:
		read1 = os.path.join(config.READS_FOLDER, basic_read_name2+"1.fq")
		read2 = os.path.join(config.READS_FOLDER, basic_read_name2+"2.fq")

		if(os.path.isfile(read1) and os.path.isfile(read2)):
			print("Aligment read: ", read1, " ", read2)
			alignment_read(read1, read2, basic_read_name2)
		else:
			print("Error not a files: ", read1, read2)			