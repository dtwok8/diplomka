""" Steps of program """
# values True/False
CREATE_READS = False
ALIGN = False
EVALUATE = True

""" create_haplotype """
REFERENCE_KIR_GENS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta"
# have to be same folder like REFERENCE_KIR_GENS_FOLDER!
REFERENCE_KIR_GENS_PSEUDOGENS_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"



""" haplotype in python dictionary
	haplotypes file will be have name like key of dictionary example "haplotype1"
	created reads also have name like key of dictionary example "haplotype1"

	example dictionary: 

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
			"haplotype_jedna" : [
				'3DL3: 00402', '3DL3:00802', '2DS2: 00101', '2DL2: 00301', '2DL3: 001', '2DP1: 00201', '2DL1: 00302', '3DP1: 007', '3DP1: 00901',
				'2DL4: 00102', '2DL4: 00501', '3DL1: 01502', '3DS1: 01301', '2DL5A: 001', '2DS5: 00201', '2DS1: 00201', '2DS4: 001', '3DL2: 0020105', '3DL2:0070102' 
							],
			"haplotype_dva" : [
				'3DL3: 00402', '3DL3:00802', '2DS2: 00101', '2DL2: 00301', '2DL3: 001', '2DP1: 00201', '2DL1: 00302', '3DP1: 007', '3DP1: 00901', '2DL4: 00102', '2DL4: 00501',
				'3DL1: 01502', '3DS1: 01301', '2DL5A: 001', '2DS5: 00201', '2DS1: 00201', '2DS4: 001', '3DL2: 0020105', '3DL2:0070102' ]
}

""" Folder for save new haplotype file """
HAPLOTYPE_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/haplotype/"

""" Folder to bowtie tools run_bowtie.py """
BOWTIE_HOME_DIRECTORY = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/bowtie2-2.4.1-linux-x86_64"


""" 
	reads create FROM ART, or from hospital
	run_art.py run_bowtie.py
"""
READS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/reads"


""" 
	index create by bowtie for reference kir gen 
"""
BOWTIE_INDEX_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/bowtie_index"
""" If bouwtie should create index TRUE/FALSE""" 
BOWTIE_BUILD_INDEX = False
BOWTIE_THREADS = 4

""" result from Bowtie """
ALIGNMENT_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/alignments"


""" bam folder for experiments """
BAM_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/bam"

""" Folder for result from evaluation alignment """
RESULT_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/result"

ALELS_STATISTICS_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/data/alels_statistics.pyc"

# Because not enought memory, use a trick and save just changes witch are lower then (average 4769)
LEVENSHTEIN_DISTANCE_CUT = 2500 #4769 
CUT_COVERAGE_ALELS = 60
CLOSE_DISTANCE = 1000

# k hovnu ne tak uplne
ALELS_DISTANCE_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/data/alels_distance.pyc"
ALELS_DISTANCE_FILE_TXT = "/home/kate/Dokumenty/FAV/Diplomka/software/data/alels_distance.txt"