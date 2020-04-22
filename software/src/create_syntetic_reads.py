import subprocess
import os

import src.create_haplotype as create_haplotype
import config

def run():
	print("creating new haplotypes ...")
	create_haplotype.run()

	haplotype_files = [f for f in os.listdir(config.HAPLOTYPE_OUTPUT_FOLDER) if os.path.isfile(os.path.join(config.HAPLOTYPE_OUTPUT_FOLDER, f))]
	print("run ART for: ", haplotype_files)

	for ha_file in haplotype_files:
		# know that file will be .fa
		output_file = os.path.join(config.READS_FOLDER, ha_file[:-3])
		process = subprocess.run(["art_illumina" ,"-ss", "MSv1", "-sam", "-i", os.path.join(config.HAPLOTYPE_OUTPUT_FOLDER, ha_file), "-p", "-l", "250", "-f", "100", "-m", "300", "-s", "10", "-o", output_file])	

	#art_illumina -ss MSv1 -sam -i amplicon_reference.fa -p -l 250 -f 100 -m 300 -s 10 -o moje_art_data$
	#haplotype = "/home/kate/Dokumenty/FAV/Diplomka/software/data/haplotype/haplotype1.fa"
	#read_output = "/home/kate/Dokumenty/FAV/Diplomka/software/data/reads/haplotype"
	#process = subprocess.run(["art_illumina" ,"-ss", "MSv1", "-sam", "-i", haplotype, "-p", "-l", "250", "-f", "100", "-m", "300", "-s", "10", "-o", read_output])
	
	#print("hello")
#exit(1)
# pair end
# 250 dlouhy ready
# misto MSv3 pouzit MSv1 protoze tak budou i data co dostanu
# -f 100 pokryti 100
# -na značí že nemá vytvořit soubor zarovnání


# ART
# art_illumina -ss MSv3 -sam -i amplicon_reference.fa -p -l 250 -f 10 -m 300 -s 10 -o moje_art_data