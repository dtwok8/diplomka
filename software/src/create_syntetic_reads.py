import subprocess
import os

import src.create_haplotype as create_haplotype
import config

"""
	Create haplotypes.
	Run ART for all files in haplotype folder.

	SET ART:
		-p paired-end read
		-l 250 bp long reads
		-MSv1 because hospital data from ilumina have this setting
		-f 100 - coverage 100%
		-na dont create alignment files
	
"""
def run():
	print("creating new haplotypes ...")
	create_haplotype.run()

	haplotype_files = [f for f in os.listdir(config.HAPLOTYPE_FOLDER) if os.path.isfile(os.path.join(config.HAPLOTYPE_FOLDER, f)) and (f!=".gitignore")]
	print("run ART for: ", haplotype_files)

	for ha_file in haplotype_files:
		# know that file will be .fa, ART want it withnout extension
		output_file = os.path.join(config.READS_FOLDER, ha_file[:-3])
		#art_illumina -ss MSv1 -sam -i amplicon_reference.fa -p -l 250 -f 100 -m 300 -s 10 -o moje_art_data$
		process = subprocess.run(["art_illumina" ,"-ss", "MSv1", "-sam", "-i", os.path.join(config.HAPLOTYPE_FOLDER, ha_file), "-p", "-l", "250", "-f", "100", "-m", "300", "-s", "10", "-o", output_file])	