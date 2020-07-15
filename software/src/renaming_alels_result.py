import os

import config


END_RESULT_FILE = '.txt'
"""
	Renaming KIR:KIR00978  to KIR2DL1*0010102 in isert string. 

	Expected format of references gens
	TAGTTTTCCATCCTTCAAATAAACATGTCTGCCCCCAT
	>KIR:KIR00978 KIR2DL1*0010102 14738 bp
	GTTCGGGAGGTTGGATCTCAGACGTGTTTTGAGTTGGTCATAGTGAAGGACACTAGGTGT

	First create dictionary like KIR:KIR00978 => KIR2DL1*0010102
	Then replace KIR:KIR00978 by KIR2DL1*0010102
""" 

""" Because there can be more files, like nuc and another"""
END_REFERENCE_GEN_FILE = "_gen.fasta"

"""
	Create dictionary like KIR:KIR00978 => KIR2DL1*0010102
	Then replace KIR:KIR00978 by KIR2DL1*0010102	
"""
def make_gen_dictionary():
	# make dictionary KIR:KIR00978 => KIR2DL1*0010102

	KIR_dictionary = dict()


	with open(config.REFERENCE_KIR_GENS_FILE, "r") as openfileobject:
		for line in openfileobject:
			if(line.startswith('>KIR')):
				split_line = line.split()
				# split_line[0][1:] - cut '>'
				KIR_dictionary[split_line[0][1:]] = split_line[1]

	print("make dictionary: ", config.REFERENCE_KIR_GENS_FILE)
	
	return KIR_dictionary
	
		
"""
	Get KIR dictionary like KIR:KIR00978 => KIR2DL1*0010102
	Get all result file
	Find KIR***** (example KIR:KIR00978)
	Find it in dictionary and replace for KIR2DL1*0010102
	Then save modified string into same file.
"""
def run(aligment_result_file_rename = None, aligment_result_file = None):
	KIR_dictionary = make_gen_dictionary()

	result_files = [f for f in os.listdir(config.RESULT_FOLDER) if os.path.isfile(os.path.join(config.RESULT_FOLDER, f)) and f.endswith(END_RESULT_FILE)]
	
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
