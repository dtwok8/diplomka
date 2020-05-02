import pickle
import sys


NUC_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_nuc.fasta"
GEN_FILE = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/KIR_gen.fasta"
DISTANCE_FILE_PYC = "/home/kate/Dokumenty/FAV/Diplomka/software/analysis/alels_distance.pyc"


def compare_gen_and_nuc_first():
	not_found = 0
	something_wrong = 0
	match = 0
	count_nuc = 0
	count_gen = 0

	KIR_alels_dictionary = dict()
	
	with open(GEN_FILE, "r") as openfileobject:
		for line in openfileobject:			
			if(line.startswith('>KIR')):
				count_gen+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>KIR'
				KIR_alels_dictionary[split_line[0][5:]] = split_line[1]


	with open(NUC_FILE, "r") as openfileobject:
		for line in openfileobject:
			# IPD:KIR00874 KIR3DS1*108
			if(line.startswith('>IPD')):
				count_nuc+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>'
				key = split_line[0][5:]
				value = split_line[1]

				if(key in KIR_alels_dictionary):
					if(KIR_alels_dictionary[key] != value):
						#print("something wrong ", key, " : ", KIR_alels_dictionary[key], "!=", value )
						something_wrong+=1
					else:
						#print("match ", key)
						match+=1
				else:
					#print("key not found: ", key)
					not_found+=1

	print("count_gen: ", count_gen, ", count_nuc: ", count_nuc, ", match: ", match, ", key not found: ", not_found, ", something wrong: ", something_wrong)


def compare_gen_and_nuc_second():
	not_found = 0
	something_wrong = 0
	match = 0
	count_nuc = 0
	count_gen = 0

	KIR_alels_dictionary = dict()
	
	with open(GEN_FILE, "r") as openfileobject:
		for line in openfileobject:			
			if(line.startswith('>KIR')):
				count_gen+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>KIR'
				KIR_alels_dictionary[split_line[1]] = split_line[0][5:]


	with open(NUC_FILE, "r") as openfileobject:
		for line in openfileobject:
			# IPD:KIR00874 KIR3DS1*108
			if(line.startswith('>IPD')):
				count_nuc+=1
				split_line = line.split()
				# split_line[0][1:] - cut '>IPD:'
				value = split_line[0][5:]
				key = split_line[1]

				if(key in KIR_alels_dictionary):
					if(KIR_alels_dictionary[key] != value):
						#print("something wrong ", key, " : ", KIR_alels_dictionary[key], "!=", value )
						something_wrong+=1
					else:
						#print("match ", key)
						match+=1
				else:
					#print("key not found: ", key)
					not_found+=1

	print("count_gen: ", count_gen, ", count_nuc: ", count_nuc, ", match: ", match, ", key not found: ", not_found, ", something wrong: ", something_wrong)


def get_min_max_distance():
	with open(DISTANCE_FILE_PYC, 'rb') as handle:
		alels_distance = pickle.load(handle)

	max_distance = 0
	min_distance = sys.maxsize 
	sum_distance = 0
	count_distance = 0

	for alel1, dic_distance in alels_distance.items():
		for alel2, distance in dic_distance.items():
			if(alel1 != alel2):
				if(max_distance < distance):
					max_distance = distance
				if(min_distance > distance):
					min_distance = distance

				sum_distance += distance
				count_distance += 1

	print("max: ", max_distance, ", min: ", min_distance, ", average_distance: ", sum_distance/count_distance)


def main(): 
	#compare_gen_and_nuc_first()
	#compare_gen_and_nuc_second()
	get_min_max_distance()



if __name__ == "__main__":
	main()