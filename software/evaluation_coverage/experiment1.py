import Bio
import pysam
import sys
import os

from Bio import AlignIO
import matplotlib.pyplot as plt


def count_statistic(sort_bam_file):

	coverage_string = pysam.depth("-aa", sort_bam_file)
	coverage_list= coverage_string.splitlines()

	curent_alel =""
	number_nuc_coverage = 0
	sum_nuc_coverage = 0
	alel_size = 0
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

			curent_alel = item_nuc_list[0]
			number_nuc_coverage = 0	
			sum_nuc_coverage = 0
			alel_size = 0
			#print(curent_alel, end =" ")


		alel_size += 1

		if(int(item_nuc_list[2]) != 0):
			number_nuc_coverage += 1
			sum_nuc_coverage += int(item_nuc_list[2])


	print(alels_statistics)
	return alels_statistics


#soubor SAM, který vytvoří bowtie
#doporučuji převést do formátu BAM (binární verze) a následně seřadit (rychleji se s tím pracuje) - to lze provést pomocí samtools
# result prejmenovat na aligments

# get all aligments sam
# tohle asi bude pak chtit prehodit do nejakeho configu
aligment_folder = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata/aligments"
bam_folder = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata/bam"

aligments_sam_files = [f for f in os.listdir(aligment_folder) if os.path.isfile(os.path.join(aligment_folder, f)) and os.path.splitext(f)[1] == '.sam']
gene_statistics = dict()

for a_sam_file in aligments_sam_files:
	#create file
	gene_name = a_sam_file.split('_')[0]
	bam_file = os.path.join(bam_folder, gene_name+".bam")
	fh = open(bam_file, 'w')
	fh.close()	

	print(a_sam_file)
	pysam.view('-o', bam_file, '-b', os.path.join(aligment_folder, a_sam_file), save_stdout=bam_file)

	sort_bam_file = os.path.join(bam_folder, gene_name+"_sort"+".bam")
	pysam.sort("-o", sort_bam_file,  bam_file)

	gene_statistics[gene_name] = count_statistic(sort_bam_file)


print(gene_statistics)

# Evaluate statistics
max_coverage1 = 0
max_coverage1_string = ""

max_coverage2 = 0
max_coverage2_string = ""

max_coverage_sum1 = 0
max_coverage_sum1_string = ""

max_coverage_sum2 = 0
max_coverage_sum2_string = ""

for gene, alels_statistics in gene_statistics.items():
	print(gene)

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

			print(alel, "c:" ,coverage, "sum: ", sum_nuc_coverage_normalized)

		#print(values['number_nuc_coverage'], values['sum_nuc_coverage'], values['alel_size'])	
		#print(alel_statistics)
	print("max_coverage1", max_coverage1_string ,max_coverage1)
	print("max_coverage2", max_coverage2_string ,max_coverage2)
	print("max_sum1" , max_coverage_sum1_string, max_coverage_sum1)
	print("max_sum2" , max_coverage_sum2_string, max_coverage_sum2)

exit(1)


# ---------------------------BORDEL------------------------------------------------------------------
sam_file = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata/result/KIR2DL1_gen_result2.sam"
output_bam_file = "output.bam"
align_bam_file = "output_alig.bam"

fh = open(output_bam_file, 'w')
fh.close()

#samtools view -Sb -o alignment.bam alignment.sam // Převedení do BAM
#rows = pysam.view("-Sb", "-o",output_bam_file, file)
pysam.view('-o', output_bam_file, '-b', sam_file, save_stdout=output_bam_file)
#print(rows)

#samtools sort -o alignment.sorted.bam alignment.bam // Seřazení

pysam.sort("-o", align_bam_file,  output_bam_file)

#samtools index SAMPLE.bam musi se vygenrovat indexy .bam.bai
#pysam.index(align_bam_file, "output_alig.bam.bai")
pysam.index(align_bam_file)

# doporučuji se podívat na samtools depth, který dokáže vytvořit "CSV" s pokrytím.

# presunout jinam a pojemnovat experiment a zkusit to na vsechny testovaci nekam si to ukladat a pak budu muset je nejak vybirat,
# pak to mozna jeste hodnotit jak to je nejvic pokryty, a jeste delat mozna procentne vzhledem k delce ty dany alely.. a
# akorat je takovy zvlastni ze alely mohou mit ruzny delky

coverage_string = pysam.depth("-aa", align_bam_file)
coverage_list= coverage_string.splitlines()

curent_alel =""
number_nuc_coverage = 0
sum_nuc_coverage = 0
alel_size = 0
alels_statistics = dict()

#sum_nuc_coverage2 = 0


for item_nuc in coverage_list:
	item_nuc_list = item_nuc.split()


	if(item_nuc_list[0] != curent_alel):
		#skip the first change
		if(curent_alel != ""):
			alels_statistics[curent_alel] = {}
			alels_statistics[curent_alel]['number_nuc_coverage'] = number_nuc_coverage
			alels_statistics[curent_alel]['sum_nuc_coverage'] = sum_nuc_coverage
			alels_statistics[curent_alel]['alel_size'] = alel_size

		curent_alel = item_nuc_list[0]
		number_nuc_coverage = 0	
		sum_nuc_coverage = 0
		alel_size = 0
		#print(curent_alel, end =" ")


	alel_size += 1

	if(int(item_nuc_list[2]) != 0):
		number_nuc_coverage += 1
		sum_nuc_coverage += int(item_nuc_list[2])


print(alels_statistics)
#print(a)
#print(de)
#for item in de:
#	print(type(item))