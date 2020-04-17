import os
from datetime import datetime

# for easy to change
REFERENCE_KIR_GENS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta"
REFERENCE_KIR_GENS_PSEUDOGENS_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"

now = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
OUTPUT_FILE = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata/haplotype/haplotyp_"+ dt_string+".sam"

HAPLOTYP = ['3DL3: 00402', '3DL3:00802', '2DS2: 00101', '2DL2: 00301', '2DL3: 001', '2DP1: 00201', '2DL1: 00302', '3DP1: 007', 
'3DP1: 00901', '2DL4: 00102', '2DL4: 00501', '3DL1: 01502', '3DS1: 01301', '2DL5A: 001', '2DS5: 00201', '2DS1: 00201', '2DS4: 001', '3DL2: 0020105', '3DL2:0070102' ]

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
			result_haplotype = result_haplotype + alel.split('\n', 1)[1]
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


	for item in HAPLOTYP:	
		# can be 3DL3:00402 or 2DS2*00101
		gen = item.replace(":", "*").split('*')

		# can be KIR3DL3 or 3DL3
		gen[0] = gen[0].strip()
		gen[1] = gen[1].strip()

		if( not gen[0].startswith('KIR')):
			gen[0] = 'KIR'+gen[0]

		gen_file_name = os.path.join(REFERENCE_KIR_GENS_FOLDER, gen[0]+"_gen.fasta")
		if(os.path.isfile(gen_file_name)):
			result_haplotype, result_haplotype_legend = search_alel_in_file(gen_file_name, gen,  result_haplotype, result_haplotype_legend)	

		else:
			# pseudogen usually in file KIR_gen.fasta
			gen_file_name2 = REFERENCE_KIR_GENS_PSEUDOGENS_FILE

			if(os.path.isfile(gen_file_name2)):
				print(gen_file_name2)
				result_haplotype, result_haplotype_legend = search_alel_in_file(gen_file_name2, gen, result_haplotype, result_haplotype_legend)	
				
				if(gen[0] not in result_haplotype_legend):
					print("Can not find", item)
			else:
				print("Can not find file ", gen_file_name, "or", gen_file_name2)


	aligment_result_file_rename = open(OUTPUT_FILE, "w+")
	aligment_result_file_rename.write(result_haplotype)
	aligment_result_file_rename.close()

	print("Haplotype legend: ", result_haplotype_legend)
	print("Create file: ", OUTPUT_FILE)


run()