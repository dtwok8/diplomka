from Bio import Seq
import Bio

my_seq = Seq.Seq("CCCCGGGAAA") # rikame mu jen ze je to nejaka sekvence
my_dna = Seq.Seq("CCCAAAGGG", Bio.Alphabet.generic_dna)
#my_rna = Bio.Seq("CCCAAAGGG", Bio.Alphabet.generic_rna)
#my_protein = Bio.Seq("AKKKGGUULL", Bio.Alphabet.generic_protein)

attributes = [a for a in dir(my_seq) if not a.startswith("_")]
print(attributes)
print(dir(my_seq))


my_gene = Seq.Seq("ACTGAC", Bio.Alphabet.generic_dna)

#get mRNA
my_transcript = my_gene.transcribe()
print(my_transcript)
print(my_transcript.alphabet)

my_protein = my_transcript.translate()
print(my_protein)
print(my_protein.alphabet)

