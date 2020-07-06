""" Steps of program """
# values True/False
CREATE_READS = 	False
ALIGN = False
EVALUATE = True

""" create_haplotype """
REFERENCE_KIR_GENS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta"
# have to be same folder like REFERENCE_KIR_GENS_FOLDER!
REFERENCE_KIR_GENS_PSEUDOGENS_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"

# predelat na tohle
REFERENCE_KIR_GENS = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"



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

# HAPLOTYPES = {
# 			"amala" : [
# 				'3DL3: 0040201', '3DL3:00802', '2DS2: 0010101', '2DL2: 0030102', '2DL3: 0010109', '2DP1: 0020108', '2DL1: 0030201', '3DP1: 007', '3DP1: 0090101',
# 				'2DL4: 0010201', '2DL4: 00501', '3DL1: 0150201', '3DS1: 0130101', '2DL5A: 00102', '2DS5: 0020101', '2DS1: 0020106', '2DS4: 0010101', '3DL2: 0020105', '3DL2:0070102' 
# 							]
# }

HAPLOTYPES = {
			"bob": [ '3DL3: 00101', '3DL3: 019', '2DS2: 0010104', '2DL2: 0030101', '2DL3: 0020102', '2DP1: 0030101', '2DL1: 0030210', '3DP1: 002', '3DP1: 0030203', '2DL4: 0010202', '2DL4: 00501', '3DL1: 002',
				'3DS1: 0130105', '2DL5A: 0010101', '2DS5: 0020104', '2DS1: 0020101', '2DS4: 0010105', '3DL2: 0020101' , '3DL2:0070102'
			],
			"amala": [
 				'3DL3: 0040201', '3DL3:00802', '2DS2: 0010101', '2DL2: 0030102', '2DL3: 0010109', '2DP1: 0020108', '2DL1: 0030201', '3DP1: 007', '3DP1: 0090101',
 				'2DL4: 0010201', '2DL4: 00501', '3DL1: 0150201', '3DS1: 0130101', '2DL5A: 00102', '2DS5: 0020101', '2DS1: 0020106', '2DS4: 0010101', '3DL2: 0020105', '3DL2:0070102' 
			],
			"cox": [
				'3DL3: 00102', '3DL3: 0090101', '2DL3: 0020101', '2DL3: 006', '2DP1: 0030102', '2DP1: 0030102', '2DL1: 0020102', '2DL1: 0020102', '3DP1: 005', '3DP1: 006', '2DL4: 00501', '2DL4: 00901', '3DL1: 0050103', 
				'3DS1: 055', '2DS5: 0020102', '2DS1: 0020105', '2DS4: 010', '3DL2: 0010301', '3DL2: 0070103'
			],
			"test1": [
				'3DL3: 0030101', '3DL3: 0140201', '2DS2: 0010104', '2DL2: 0030105', '2DL3: 0020101', '2DL5B: 01301', '2DS3: 0020102', '2DS3: 0020102', '2DP1: 0030101' , '2DP1:0010203', '2DL1: 0030203', 
				'2DL1:007', '3DP1: 004', '3DP1: 004', '2DL4: 0080104', '2DL4: 010', '3DL1: 0150101', '3DS1: 014', '2DL5A: 00102', '2DS1: 0020104', '2DS4: 0010103', '3DL2: 00501', '3DL2: 00501'
			],
			"test2": [
				'3DL3: 0090102', '3DL3: 019', '2DL2: 0010101', '2DL3: 0010111', '2DL5B: 0020101', '2DP1: 0020107', '2DP1: 0030102', '2DL1: 0020102', '2DL1: 0030210', '3DP1: 004', '3DP1: 01001', '2DL4: 0080104', '2DL4: 0080104',
				'3DL1: 0070101', '3DS1: 078', '2DL5A: 0050101', '2DS5: 010', '2DS1: 0020102', '2DS4:0010103', '3DL2: 0020101', '3DL2: 00501'	
			],
			"test3": [
				'3DL3: 005', '3DL3: 0140201', '2DL3: 0010101', '2DL3: 0020103', '2DP1: 0020108', '2DP1: 0020108', '2DL1: 0040101', '2DL1: 008', '3DP1: 0030102', '3DP1: 00902', '2DL4: 0010306', '2DL4: 00501', '3DL1: 002', 
				'3DL1: 0040101', '3DL2: 0010301', '3DL2: 008'	
			],
			"test4": [
				'3DL3: 0030104', '3DL3: 007', '2DS2: 0010105', '2DL2: 0030101', '2DL3: 0010102', '2DP1: 008', '2DL1: 007', '3DP1: 007', '3DP1: 00902', '2DL4: 0010307', '2DL4: 0080104', '3DL1: 0150202', '3DS1: 055', '2DS5: 007', 
				'2DS1: 0020101', '2DS4: 0060101', '3DL2: 0020101', '3DL2: 00903'
			],
			"test5": [
				'3DL3: 0140202', '3DL3: 036', '2DL3: 0010109', '2DL3: 006', '2DS3: 0010301', '2DP1: 0030102', '2DP1: 009', '2DL1: 0030208', '2DL1: 00303', '3DP1: 001', '3DP1: 002', '2DL4: 0010202', '2DL4: 0010202', '3DL1: 0200101',
				'3DS1: 0130102', '2DL5A: 0010102', '2DS1: 0020105', '3DL2: 00202', '3DL2: 018'	
			],
			"test6": [
				'3DL3: 0090102', '3DL3: 0140203', '2DS2: 0010111', '2DL2: 0010105', '2DL3: 0010102', '2DL5B: 0080101', '2DS3: 0010302', '2DP1: 0020103', '2DP1: 010', '2DL1: 0030203', '2DL1: 0040102', '3DP1: 0030202', '3DP1: 0030402',
				'2DL4: 0010303', '2DL4: 00901', '3DL1: 0050102', '3DL1: 0250102', '2DS4: 0010104', '2DS4: 010', '3DL2: 0010302', '3DL2: 01001'	
			],
			"test7": [
				'3DL3: 00802', '3DL3: 0090103', '2DL3: 0010103', '2DL3: 0010108', '2DP1: 0020106', '2DP1: 004', '2DL1: 0030204', '2DL1: 0030205', '3DP1: 0030202', '3DP1: 0030202', '2DL4: 0010201', '2DL4: 0010305', '3DL1: 008',
				'3DL1: 0150203', '2DS4: 0010107', '2DS4: 0030104', '3DL2: 0020105', '3DL2: 00901'	
			],
			"test8": [
				'3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DL5B: 0070101', '2DS3: 0020101', '2DS3: 0020101', '2DS3: 0010302', '2DP1: 0030102', '2DL1: 00402', '3DP1: 0030101',
				'3DP1: 005', '2DL4: 00104', '2DL4: 0080104', '3DS1: 0130104', '3DS1: 055', '2DL5A: 0050102', '2DL5A: 0050102', '2DS1: 0020102', '2DS1: 0020105', '3DL2: 0010102', '3DL2: 0070102'
			],
			"test9": [
				'3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DL5B: 0070101', '2DS3: 0020101', '2DS3: 0020101', '2DP1: 0030102', '2DL1: 00402', '3DP1: 0030101', '3DP1: 005',
				'2DL4: 00104', '2DL4: 0080104', '3DL1: 0150208', '3DS1: 0130104', '2DL5A: 01201', '2DS1: 0020102', '2DS4: 0040101', '3DL2: 0010102', '3DL2: 0070102'	
			],
			"test10": [
				'3DL3: 0030101', '3DL3: 0140201', '2DS2: 0010104', '2DL2: 0030105', '2DL3: 0020101', '2DL5B: 01301', '2DS3: 0020102', '2DS3: 0020102', '2DP1: 0010203', '2DP1: 0030101', '2DL1: 0030203', '2DL1: 007', '3DP1: 004',
				'3DP1: 004', '2DL4: 0080104', '2DL4: 010', '3DL1: 0150101', '3DS1: 014', '2DS1: 0020104', '3DL2: 00501', '3DL2: 00501'
			],
			"test11": [
				'3DL3: 0030103', '3DL3: 00601', '2DS2: 0010103', '2DS2: 0010112', '2DL2: 0010102', '2DL2: 0030101', '2DP1: 0030102', '2DP1: 008', '2DL1: 00402', '2DL1: 00402', '3DP1: 0030101', '3DP1: 005', '2DL4: 00104', '2DL4: 0080104',
				'3DS1: 0130104', '3DS1: 055', '2DL5A: 0050102', '2DL5A: 0050102', '2DS5: 0020102', '2DS5: 0020103', '2DS1: 0020102', '2DS1: 0020105', '3DL2: 0010102', '3DL2: 0070102'	
			],
			"ho301": [
				'3DL3: 0140201', '3DL3: 0140201', '2DS2: 0010104', '2DS2: 0010106', '2DL2: 0010103', '2DL2: 0030107', '2DL5B: 010', '2DL5B: 010', '2DS3: 0010301', '2DS3: 0020103', '2DP1: 0010202', '2DP1: 0010202', '2DL1: 00402', 
				'2DL1: 010', '3DP1: 0030101', '3DP1: 004', '2DL4: 0010201', '2DL4: 0010201', '3DL1: 002', '3DL1: 002', '2DS4: 0010109', '2DS4: 0010109', '3DL2: 0020102', '3DL2: 0020106'
			],
			"jvm": [
				'3DL3: 00801', '3DL3: 0140201', '2DS2: 0010110', '2DL2: 0030102', '2DL3: 010', '2DP1: 004', '2DL1: 0030203', '3DP1: 001', '3DP1: 0030202', '2DL4: 0010304', '2DL4: 0080101', '3DL1: 0010104', '3DL1: 008', '2DS4: 0030103',
				'2DS4: 0030103', '3DL2: 0010101', '3DL2:018'	
			],
			"kas011": [
				'3DL3: 0090101', '3DL3: 0140203', '2DL3: 0020103', '2DL3: 0020103', '2DP1: 0020104', '2DP1: 0030101', '2DL1: 0020101', '2DL1: 0030209', '3DP1: 0030206', '3DP1: 009', '2DL4: 0010301', '2DL4: 00501', '3DL1: 008', 
				'3DS1: 013011', '2DL5A: 0010102', '2DS5: 0020101', '2DS1: 0020101', '2DS4: 0030101', '3DL2: 01001', '3DL2: 018'	
			],
			"olga": [
				'3DL3: 00201', '3DL3: 00202' , '2DL3: 0010105', '2DL3: 0010105', '2DP1: 0020105', '2DP1: 006', '2DL1: 0030204', '2DL1: 0030204', '3DP1: 0030201', '3DP1: 0030201', '2DL4: 00501', '2DL4: 00901', '3DL1: 0010102', 
				'3DL1: 0050101', '3DS1: 0130107', '2DL5A: 00103', '2DS5: 0020103', '2DS1: 0020101', '2DS4: 010', '3DL2: 0070101', '3DL2: 0070102'
			],
			"rsh": [
				'3DL3: 00202', '3DL3: 0040202', '2DS2: 0010108', '2DL2: 0030104', '2DL3: 0010107', '2DL5B: 004', '2DP1: 0020110', '2DP1: 009', '2DL1: 0030205', '2DL1: 01201', '3DP1: 0030401', '3DP1: 008', '2DL4: 0010307', 
				'2DL4: 00901', '3DL1: 0050101', '3DL1: 01701', '2DS5: 006', '2DS4: 0060102', '3DL2: 023', '3DL2: 056'
			]
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

ALELS_STATISTICS_FOLDER = "/home/kate/Dokumenty/FAV/Diplomka/software/data/statistics"
ALELS_STATISTICS_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/data/alels_statistics.pyc"

# Because not enought memory, use a trick and save just changes witch are lower then (average 4769)
LEVENSHTEIN_DISTANCE_CUT = 2500 #4769 
CUT_COVERAGE_ALELS = 75
CLOSE_DISTANCE = 100

# k hovnu ne tak uplne
ALELS_DISTANCE_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/data/alels_distance.pyc"
ALELS_DISTANCE_FILE_TXT = "/home/kate/Dokumenty/FAV/Diplomka/software/data/alels_distance.txt"