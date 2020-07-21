import pysam
import sys
import os

import config

"""
	Make statistics like coverage of allels and sum of histogram covarage.
"""
def count_statistic(sort_bam_file: str):
	coverage_string = pysam.depth("-aa", sort_bam_file)
	coverage_list= coverage_string.splitlines()

	curent_alel =""
	number_nuc_coverage = 0
	sum_nuc_coverage = 0
	alel_size = 0
	max_histogram = 0
	alels_statistics = dict()


	for item_nuc in coverage_list:
		item_nuc_list = item_nuc.split()


		if(item_nuc_list[0] != curent_alel):
			#skip the first change
			if(curent_alel != ""):
				alels_statistics[curent_alel] = {}
				alels_statistics[curent_alel]['number_nuc_coverage'] = number_nuc_coverage
				alels_statistics[curent_alel]['sum_nuc_coverage'] = sum_nuc_coverage
				alels_statistics[curent_alel]['alel_size'] = alel_size
				alels_statistics[curent_alel]['max_histogram'] = max_histogram

			curent_alel = item_nuc_list[0]
			number_nuc_coverage = 0	
			sum_nuc_coverage = 0
			alel_size = 0
			max_histogram = 0

		alel_size += 1

		if(int(item_nuc_list[2]) != 0):
			number_nuc_coverage += 1
			sum_nuc_coverage += int(item_nuc_list[2])

			if(max_histogram < int(item_nuc_list[2])):
				max_histogram = int(item_nuc_list[2])

	return alels_statistics


"""
	Evaluate statistics, find max coverage of each gene
	and write it into result file.
"""
def evaluate_statistics(haplotype_gen_statistics: dict, output_file_name: str):
	output_file = open(output_file_name, "w")
	max_coverage1 = 0
	max_coverage1_string = ""

	max_coverage2 = 0
	max_coverage2_string = ""

	max_coverage_sum1 = 0
	max_coverage_sum1_string = ""

	max_coverage_sum2 = 0
	max_coverage_sum2_string = ""

	max_histogram1 = 0
	max_histogram1_string = "" 

	max_histogram2 = 0
	max_histogram2_string = "" 


	for gene, alels_statistics in haplotype_gen_statistics.items():
		print(gene, file=output_file)

		max_coverage1 = 0
		max_coverage1_string = ""

		max_coverage2 = 0
		max_coverage2_string = ""

		max_coverage_sum1 = 0
		max_coverage_sum1_string = ""

		max_coverage_sum2 = 0
		max_coverage_sum2_string = ""

		for alel, values in alels_statistics.items():
			if(values['number_nuc_coverage'] > 0):
				coverage =  values['number_nuc_coverage'] / (values['alel_size'] / 100)
				if(coverage > max_coverage1):
					max_coverage2 = max_coverage1
					max_coverage2_string = max_coverage1_string	

					max_coverage1 = coverage
					max_coverage1_string = alel	
				elif (coverage > max_coverage2):
					max_coverage2 = coverage
					max_coverage2_string = alel	

				sum_nuc_coverage_normalized = values['sum_nuc_coverage']
				if(sum_nuc_coverage_normalized > max_coverage_sum1):
					max_coverage_sum2 = max_coverage_sum1
					max_coverage_sum2_string = max_coverage_sum1_string

					max_coverage_sum1 = sum_nuc_coverage_normalized
					max_coverage_sum1_string = alel

				elif (sum_nuc_coverage_normalized > max_coverage_sum2):
					max_coverage_sum2 = sum_nuc_coverage_normalized
					max_coverage_sum2_string = alel	

				if(max_histogram1 < values['max_histogram']):
					max_histogram2 = max_histogram1
					max_histogram2_string = max_histogram1_string

					max_histogram1 = values['max_histogram']
					max_histogram1_string = alel
				
				elif (max_histogram2 < values['max_histogram']):
					max_histogram2 = values['max_histogram']
					max_histogram2_string = alel

				print(alel, " c: ", coverage, " sum: ", sum_nuc_coverage_normalized, "max histogram: ", values['max_histogram'] , file=output_file)

		print("max_coverage1 ", max_coverage1_string, max_coverage1 , file=output_file)
		print("max_coverage2 ",  max_coverage2_string, max_coverage2 , file=output_file)
		print("max_sum1 ", max_coverage_sum1_string, max_coverage_sum1 , file=output_file)
		print("max_sum2 ", max_coverage_sum2_string, max_coverage_sum2 , file=output_file)
		print("max_histogram1 ", max_histogram1_string, max_histogram1 , file=output_file)
		print("max_histogram2 ", max_histogram2_string, max_histogram2 , file=output_file)
	output_file.close()
	print("create result file ", output_file_name)


""""
	Find all haplotypes in sam folder.
"""
def find_haplotypes(aligments_sam_files: list):
	gene_statistics = dict()
	haplotype_list = list()
	for file in aligments_sam_files:
		file_name_split = file.split('_')
		# find haplotype name 
		haplotype_name = ""	
		for item in file_name_split:
			if(item.startswith("KIR")):
				break
			else:
				haplotype_name += item+"_"
		haplotype_list.append(haplotype_name)

	# find unique
	haplotype_name_set = set(haplotype_list)
	haplotype_name_list_unique = list(haplotype_name_set)
	return haplotype_name_list_unique


""""
	conver to BAM (binary version) because it is faster
	then run statistics,
	then eva
"""
def run():
	# get all alignment .sam
	aligments_sam_files = [f for f in os.listdir(config.ALIGNMENT_FOLDER) if os.path.isfile(os.path.join(config.ALIGNMENT_FOLDER, f)) and os.path.splitext(f)[1] == '.sam']

	haplotype_name_list_unique = find_haplotypes(aligments_sam_files)
	print("Find haplotype: ", haplotype_name_list_unique)

	for haplotype in haplotype_name_list_unique:
		print(haplotype)
		algiments_file_haplotype = [f for f in aligments_sam_files if f.startswith(haplotype)]
		haplotype_gen_statistics = {}
		
		for sam_file in algiments_file_haplotype:
			file_name_split = sam_file.split('_')
			
			# find KIR 
			kir_name = ""
			
			for item in file_name_split:
				if(item.startswith("KIR")):
					kir_name= item
					break

			print("preparing file for ", haplotype, kir_name, "...")
			# create bam file and  conver sam to bam and then sort 
			bam_file = os.path.join(config.BAM_FOLDER, haplotype+kir_name+".bam")
			fh = open(bam_file, 'w')
			fh.close()

			pysam.view('-o', bam_file, '-b', os.path.join(config.ALIGNMENT_FOLDER, sam_file), save_stdout=bam_file)

			sort_bam_file = os.path.join(config.BAM_FOLDER, haplotype+kir_name+"_sort"+".bam")
			pysam.sort("-o", sort_bam_file,  bam_file)

			#count statistics
			print("counting ", haplotype, kir_name, "...")

			haplotype_gen_statistics[kir_name] = count_statistic(sort_bam_file)

		result_file = os.path.join(config.RESULT_FOLDER, haplotype+"exp2.txt")
		evaluate_statistics(haplotype_gen_statistics, result_file)	