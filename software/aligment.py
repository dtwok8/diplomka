import Bio
import pysam
import sys

from Bio import AlignIO
import matplotlib.pyplot as plt

#soubor SAM, který vytvoří bowtie
#doporučuji převést do formátu BAM (binární verze) a následně seřadit (rychleji se s tím pracuje) - to lze provést pomocí samtools

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
bam = pysam.AlignmentFile(align_bam_file, 'rb')

#for line in pysam.depth(align_bam_file): ??
#    print(line)

#bam.count_coverage(contig, start=None, stop=None, region=None, quality_threshold=15, read_callback='all', reference=None, end=None)
asd = bam.count_coverage("KIR:KIR00005") # >KIR:KIR00005 KIR2DL1*00303 14740 bp
# four array.arrays of the same length in order A C G T
print(type(asd)) # tuple
print(len(asd)) # 4
print(len(asd[1])) # 14 740
print(type(asd[1])) # array. array
print(asd[0][50:100]) # 'L', [cisla ]
print(asd[1][50:100])
print(asd[2][50:100])
print(asd[3][50:100])


asd = bam.count_coverage("KIR:KIR00001") # KIR:KIR00001 KIR2DL1*0010101 14738 bp
print(type(asd)) # tuple
print(len(asd)) # 4
print(len(asd[1])) # 14 738
print(type(asd[1])) # array. array


#bamfile = sys.argv[1] # Input bam file
samfile = pysam.Samfile(align_bam_file, "rb" )
samfile.pileup(max_depth = 1000000)
X = []
#for pileupcolumn in samfile.pileup(max_depth = 1000000): 
for pileupcolumn in samfile.pileup("KIR:KIR00005", 0, 14740):      
    X.append(pileupcolumn.n)

Y = range(len(X))
print(len(X))
print(len(Y))
plt.title(align_bam_file)
plt.xlabel('Bases')
plt.ylabel('Coverage')
plt.plot(Y, X,c = 'g')
plt.grid(True, lw = 2, ls = '-', c = '.75')
plt.xlim(0, 14740)
plt.axhline(y=230, color='r', linestyle='--')
plt.savefig('Coverage_plot11c2.png', c = 'k')
#print(asd)
# doporučuji se podívat na samtools depth, který dokáže vytvořit "CSV" s pokrytím.

# presunout jinam a pojemnovat experiment a zkusit to na vsechny testovaci nekam si to ukladat a pak budu muset je nejak vybirat,
# pak to mozna jeste hodnotit jak to je nejvic pokryty, a jeste delat mozna procentne vzhledem k delce ty dany alely.. a
# akorat je takovy zvlastni ze alely mohou mit ruzny delky

coverage_string = pysam.depth("-aa", align_bam_file)
coverage_list= coverage_string.splitlines()

curent_alel =""
sum_nuc_coverage = 0
#sum_nuc_coverage2 = 0

for item_nuc in coverage_list:
	item_nuc_list = item_nuc.split()


	if(item_nuc_list[0] != curent_alel):
		print(sum_nuc_coverage)
		curent_alel = item_nuc_list[0]
		sum_nuc_coverage = 0	
		sum_nuc_coverage2 = 0
		print(curent_alel, end =" ")


	#print(item_nuc_list[2], end =" ")
	if(int(item_nuc_list[2]) != 0):
		sum_nuc_coverage += 1

		#sum_nuc_coverage2 += item_nuc_list[2]

#print(a)
#print(de)
#for item in de:
#	print(type(item))

exit(1)

samfile = pysam.AlignmentFile(align_bam_file, "rb" )
for pileupcolumn in samfile.pileup("KIR:KIR00001", 100, 14740):
    print ("\ncoverage at base %s = %s" %
           (pileupcolumn.pos, pileupcolumn.n))
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
            print ('\tbase in read %s = %s' %
                  (pileupread.alignment.query_name,
                   pileupread.alignment.query_sequence[pileupread.query_position]))

samfile.close()


#samfile = pysam.AlignmentFile(file, "rb")


#alignment = list(AlignIO.parse("/home/kate/Dokumenty/FAV/Diplomka/software/mydata/result/KIR2DL1_gen_result2.sam", "emboss"))
#print(alignment)


#print("Alignment length %i" % alignment.get_alignment_length())
#for record in alignment :
#    print(record.seq + " " + record.id)