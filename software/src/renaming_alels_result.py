import os

# my file
import config

"""
	Renaming KIR:KIR00978  to KIR2DL1*0010102 in isert string. 

	Expected format of references gens
	TAGTTTTCCATCCTTCAAATAAACATGTCTGCCCCCAT
	>KIR:KIR00978 KIR2DL1*0010102 14738 bp
	GTTCGGGAGGTTGGATCTCAGACGTGTTTTGAGTTGGTCATAGTGAAGGACACTAGGTGT

	First create dictionary like KIR:KIR00978 => KIR2DL1*0010102
	Then replace KIR:KIR00978 by KIR2DL1*0010102
""" 


def make_gen_dictionary(reference_kir_gens_folder: str):
	# make dictionary KIR:KIR00978 => KIR2DL1*0010102
	onlyfiles = [f for f in os.listdir(reference_kir_gens_folder) if os.path.isfile(os.path.join(reference_kir_gens_folder, f)) and f.endswith("_gen.fasta") ]

	KIR_dictionary = dict()

	for i in range(0,len(onlyfiles)):
		with open(reference_kir_gens_folder+"/"+onlyfiles[i], "r") as openfileobject:
			for line in openfileobject:
				if(line.startswith('>KIR')):
					split_line = line.split()
					# split_line[0][1:] - cut '>'
					KIR_dictionary[split_line[0][1:]] = split_line[1]

		print("make dictionary: ", onlyfiles[i])
	
	return KIR_dictionary
		

def run(aligment_result_file_rename = None, aligment_result_file = None):
	KIR_dictionary = make_gen_dictionary(config.REFERENCE_KIR_GENS_FOLDER)

	result_files = [f for f in os.listdir(config.RESULT_FOLDER) if os.path.isfile(os.path.join(config.RESULT_FOLDER, f))]
	
	for file in result_files:
		aligment_result_file = open(os.path.join(config.RESULT_FOLDER, file), "r")
		aligment_result_str = aligment_result_file.read()
		aligment_result_file.close()
		print("replacing...", file)

		for key, item in KIR_dictionary.items():
			aligment_result_str = aligment_result_str.replace(key.strip(), item)

		aligment_result_file_rename = open(os.path.join(config.RESULT_FOLDER, file), "w")
		aligment_result_file_rename.write(aligment_result_str)
		aligment_result_file_rename.close()