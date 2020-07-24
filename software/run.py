import config 

import src.create_syntetic_reads as create_syntetic_reads 
import src.alignment_reads_to_reference as alignment_reads_to_reference

import src.alleles_identification.count_levenshtein as count_levenshtein
import src.alleles_identification.exp1 as exp1
import src.alleles_identification.exp2_more_alignment as exp2_more_alignment
import src.alleles_identification.exp3_clusters as exp3_clusters

import src.renaming_alels_result as renaming_alels_result

"""
	RUN all scripts
	1. create syntetics_reads - create haplotype + ART
	2. align - bowtie
	3. indetify - identify alels by align and another by steps
"""
def main():  
	if(config.CREATE_READS):
		create_syntetic_reads.run()

	if(config.ALIGN):
		alignment_reads_to_reference.run()

	if(config.IDENTIFY):
		print("identify")

		if(config.PRECOMPUTATION_DISTANCE):
			count_levenshtein.run()
		if(config.EXP1):
			exp1.run()
		if(config.EXP2):
			exp2_more_alignment.run()
		if(config.EXP3):
			exp3_clusters.run()


if __name__ == "__main__":
	main()