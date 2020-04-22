import subprocess
import os

import config


"""
	Get all references gens and make index.
	Prepare for aligmnets.

"""
def bowtie_build_index():
	onlyfiles = [f for f in os.listdir(config.REFERENCE_KIR_GENS) if os.path.isfile(os.path.join(config.REFERENCE_KIR_GENS, f))]
	
	for i in range(0,len(onlyfiles)):
		name = os.path.splitext(os.path.basename(onlyfiles[i]))[0]
		# cw= index_foler - change folder where command will be execute, because there wasnt any parameter for output folder
		process = subprocess.run([config.BOWTIE2_HOME_DIRECTORY+"/bowtie2-build", config.REFERENCE_KIR_GENS+"/"+onlyfiles[i], name], cwd=config.BOWTIE_INDEX_FOLDER)
		print("aligment build index", name)


"""
	Compare all KIR gens with reads.
"""
def alignment_read(read1: str, read2: str, basic_read_name: str):
	onlyfiles = [f for f in os.listdir(config.BOWTIE_INDEX_FOLDER) if os.path.isfile(os.path.join(config.BOWTIE_INDEX_FOLDER, f))]
	#print(onlyfiles)

	for i in range(0,len(onlyfiles)):
		name = os.path.basename(onlyfiles[i]).split('.')[0]
		print("align: ", basic_read_name," ", name)
		process = subprocess.run([config.BOWTIE2_HOME_DIRECTORY+"/bowtie2","-x", name, "-1", read1,"-2", read2, "-S", os.path.join(config.ALIGNMENT_FOLDER, basic_read_name+"_"+name+"_align.sam")], cwd=config.BOWTIE_INDEX_FOLDER) 	



"""
	Run make index
	For all reads in folder reads run align for references gens
		- reads names has to be parser
"""
def run():
	if(config.BOWTIE_BUILD_INDEX):
		bowtie_build_index()
	
	all_reads = [f for f in os.listdir(config.READS_FOLDER) if os.path.isfile(os.path.join(config.READS_FOLDER, f)) and f.endswith(".fq") ]

	# we know that it fill end 1.fq a 2.fq
	# get basic name, unique withnout 1.fq and 2.fq
	read_basic_list = list()
	for read in all_reads:
		print(read)
		read_basic_name = read[:-4]

		read_basic_list.append(read_basic_name)	

	unique_set = set(read_basic_list)
	unique_read_list = list(unique_set)

	for basic_read_name2 in unique_read_list:
		read1 = os.path.join(config.READS_FOLDER, basic_read_name2+"1.fq")
		read2 = os.path.join(config.READS_FOLDER, basic_read_name2+"2.fq")

		if(os.path.isfile(read1) and os.path.isfile(read2)):
			print("Aligment read: ", read1, " ", read2)
			alignment_read(read1, read2, basic_read_name2)
		else:
			print("Error not a files: ", read1, read2)			