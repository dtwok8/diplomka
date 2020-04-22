""" create_haplotype """
REFERENCE_KIR_GENS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta"
REFERENCE_KIR_GENS_PSEUDOGENS_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"

HAPLOTYPE_OUTPUT_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/haplotype/"

""" haplotype in python dictionary example: 
HAPLOTYPES = {
			"haplotype1" : [
				'3DL3: 00402', '3DL3:00802', '2DS2: 00101', '2DL2: 00301', '2DL3: 001', '2DP1: 00201', '2DL1: 00302', '3DP1: 007', '3DP1: 00901',
				'2DL4: 00102', '2DL4: 00501', '3DL1: 01502', '3DS1: 01301', '2DL5A: 001', '2DS5: 00201', '2DS1: 00201', '2DS4: 001', '3DL2: 0020105', '3DL2:0070102' 
							],
			"haplotype2" : [
				'3DL3: 00402', '3DL3:00802', '2DS2: 00101', '2DL2: 00301', '2DL3: 001', '2DP1: 00201', '2DL1: 00302', '3DP1: 007', '3DP1: 00901', '2DL4: 00102', '2DL4: 00501',
				'3DL1: 01502', '3DS1: 01301', '2DL5A: 001', '2DS5: 00201', '2DS1: 00201', '2DS4: 001', '3DL2: 0020105', '3DL2:0070102' ]
}
"""

HAPLOTYPES = {
			"haplotype1" : [
				'3DL3: 00402', '3DL3:00802', '2DS2: 00101', '2DL2: 00301', '2DL3: 001', '2DP1: 00201', '2DL1: 00302', '3DP1: 007', '3DP1: 00901',
				'2DL4: 00102', '2DL4: 00501', '3DL1: 01502', '3DS1: 01301', '2DL5A: 001', '2DS5: 00201', '2DS1: 00201', '2DS4: 001', '3DL2: 0020105', '3DL2:0070102' 
							],
			"haplotype2" : [
				'3DL3: 00402', '3DL3:00802', '2DS2: 00101', '2DL2: 00301', '2DL3: 001', '2DP1: 00201', '2DL1: 00302', '3DP1: 007', '3DP1: 00901', '2DL4: 00102', '2DL4: 00501',
				'3DL1: 01502', '3DS1: 01301', '2DL5A: 001', '2DS5: 00201', '2DS1: 00201', '2DS4: 001', '3DL2: 0020105', '3DL2:0070102' ]
}


# renaming_alels_result.py
REFERENCE_KIR_GENS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta"


# run_bowtie.py
BOWTIE2_HOME_DIRECTORY = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/bowtie2-2.4.1-linux-x86_64"

# v tehle slozce by mel byt read vytvoreny artem

# reads create FROM ART, or from hospital
# run_art.py run_bowtie.py

READS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata"
#referencni kir geny z https://github.com/ANHIG/IPDKIR/tree/Latest/fasta
# vybiram si jen ty soubory ktere konci na _gen.fasta

REFERENCE_KIR_GENS = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/gen"
reference_kir_gens = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/gen"


# index create by bowtie for reference kir gens
BOWTIE_INDEX_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/index" 


#vysledne soubory z bowtie
results = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata/result" 



# renaming_alels_result.py
ALIGMENT_RESULT_FILE = "/home/kate/Dokumenty/FAV/Diplomka/software/data/result/vysledky.txt"

# renaming_alels_result.py
ALIGMENT_RESULT_FILE_RENAME = "/home/kate/Dokumenty/FAV/Diplomka/software/data/result/vysledky2.txt"
