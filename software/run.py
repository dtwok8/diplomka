import config 

import src.create_syntetic_reads as create_syntetic_reads 
import src.alignment_reads_to_reference as alignment_reads_to_reference

#import src.evaluation_alignment.count_levenshtein as count_levenshtein
import src.evaluation_alignment.experiment1 as experiment1
import src.evaluation_alignment.experiment2 as experiment2
import src.evaluation_alignment.experiment3 as experiment3


import src.renaming_alels_result as renaming_alels_result

"""
	RUN all scripts
	1. create syntetics_reads - create haplotype + ART
	2. align - bowtie
	3. evaluate - make result from evaluation of align file
"""
def main():  
	if(config.CREATE_READS):
		create_syntetic_reads.run()

	if(config.ALIGN):
		alignment_reads_to_reference.run()

	if(config.EVALUATE):
		print("eval")
		#count_levenshtein.run()
		experiment1.run()
		#experiment3.run()
		#analysis.run()
		#renaming_alels_result.run()


if __name__ == "__main__":
	main()