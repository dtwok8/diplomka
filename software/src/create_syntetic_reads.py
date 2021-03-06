import subprocess
import os

import src.create_genome as create_genome
import config


GENOME_END_FILE = '.fa'
"""
	Create genome.
	Run ART for all files in genome folder.

	SET ART:
		-p paired-end read
		-l 250 bp long reads
		-MSv1 because hospital data from ilumina have this setting
		-f 100 - coverage 100%
		-na dont create alignment files
	
"""
def run():
	print("creating new genomes ...")
	create_genome.run()

	genome_files = [f for f in os.listdir(config.GENOME_FOLDER) if os.path.isfile(os.path.join(config.GENOME_FOLDER, f)) and (f.endswith(GENOME_END_FILE))]
	print("run ART for: ", genome_files)

	for ge_file in genome_files:
		# know that file will be .fa, ART want it withnout extension
		output_file = os.path.join(config.READS_FOLDER, ge_file[:-3])
		#art_illumina -ss MSv1 -sam -i amplicon_reference.fa -p -l 250 -f 100 -m 300 -s 10 -o moje_art_data$
		process = subprocess.run(["art_illumina" ,"-ss", "MSv1", "-sam", "-i", os.path.join(config.GENOME_FOLDER, ge_file), "-p", "-l", "250", "-f", "100", "-m", "300", "-s", "10", "-o", output_file])	