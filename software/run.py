import sys, getopt

import src.create_syntetic_reads as create_syntetic_reads 
import src.alignment_reads_to_reference as alignment_reads_to_reference
import src.evaluation_alignment.experiment1 as experiment1
import src.evaluation_alignment.experiment2 as experiment2
import src.evaluation_alignment.experiment3 as experiment3
import src.renaming_alels_result as renaming_alels_result


def main(argv):
	create_reads = False
	align = False
	evaluate = True
   
	try:
		# : when is required
		opts, args = getopt.getopt(argv,"cae",["create_reads=","align=", "evaluate"])
	except getopt.GetoptError:
		print('run.py -c <create_reads-True/False> -a <align - True/False> -e <evaluate - True/False> ' )
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print('run.py - run with default seting, all settings are true and it is on all file in folders ')
			print('run.py -c <create_reads-True/False> -a <align - True/False> -e <evaluate - True/False>' )
			sys.exit()
		elif opt in ("-c", "--create_reads"):
			create_reads = arg
		elif opt in ("-a", "--align"):
			align = arg
		elif opt in ("-e", "--evaluate"):
			evaluate = arg

	if(create_reads):
		create_syntetic_reads.run()
		print("create_reads")

	if(align):
		alignment_reads_to_reference.run()

	if(evaluate):
		#experiment1.run()
		experiment3.run()
		renaming_alels_result.run()

if __name__ == "__main__":
	main(sys.argv[1:])