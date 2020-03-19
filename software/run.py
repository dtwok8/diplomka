import subprocess
import os




# !!!! ještě mě tak napadá ty máš identifikovat KIR alely takže možná budeš muset ještě roztříhat ty soubory na alely

# pair end
# 250 dlouhy ready
# misto MSv3 pouzit MSv1 protoze tak budou i data co dostanu
# -f 100 pokryti 100
# -na značí že nemá vytvořit soubor zarovnání


# ART
# art_illumina -ss MSv3 -sam -i amplicon_reference.fa -p -l 250 -f 10 -m 300 -s 10 -o moje_art_data
#tak já teoreticky můžu zkusit vytvořit haplotypy a porovnávat to s nima ne? 


# nejdřív to musíš sbuldit všechny ty referenční geny
bowtie2_home_directory = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/bowtie2-2.4.1-linux-x86_64"
# v tehle slozce by mel byt read vytvoreny artem
ready = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata"
#referencni kir geny z https://github.com/ANHIG/IPDKIR/tree/Latest/fasta
# vybiram si jen ty soubory ktere konci na _gen.fasta
reference_kir_gens = "/home/kate/Dokumenty/FAV/Diplomka/existujicisw/referencni/IPDKIR-Latest/fasta/gen"
#vysledne soubory z bowtie
results = "/home/kate/Dokumenty/FAV/Diplomka/software/mydata/result" 
#indexi vytvorene bowtie pro refenrenci kir geny
index_folder = "/home/kate/Dokumenty/FAV/Diplomka/software/index"

#build index
# cw= index_foler - change folder where command will be execute, because there wasnt any parameter for output folder
# nejdriv zbuldeni vsech referencnich genu
onlyfiles = [f for f in os.listdir(reference_kir_gens) if os.path.isfile(os.path.join(reference_kir_gens, f))]
for i in range(0,len(onlyfiles)):
	name = os.path.splitext(os.path.basename(onlyfiles[i]))[0]
	#process = subprocess.run([bowtie2_home_directory+"/bowtie2-build", reference_kir_gens+"/"+onlyfiles[i], name], cwd=index_folder)
	print("build", name)	

#porovnavat vsechny kir geny z danym readem
onlyfiles = [f for f in os.listdir(index_folder) if os.path.isfile(os.path.join(index_folder, f))]
print(onlyfiles)

for i in range(0,len(onlyfiles)):
	name = os.path.basename(onlyfiles[i]).split('.')[0]
	print("align", name)
	process = subprocess.run([bowtie2_home_directory+"/bowtie2","-x", name, "-1", ready+"/kir-genotyp1.fq","-2",ready+"/kir-genotyp2.fq", "-S", results+"/"+name+"_result2.sam"], cwd=index_folder) 	


exit(1)
# bordel
process = subprocess.run([bowtie2_home_directory+"/bowtie2-build", reference_kir_gens+"/KIR2DL1_gen.fasta", "kir2"])
process = subprocess.run([bowtie2_home_directory+"/bowtie2","-x", "kir2", "-1", ready+"/moje_art_data1.fq","-2",ready+"/moje_art_data2.fq", "-S", results+"/result2.sam"]) 
print(process)
process.stdout
print("konec")
#$BT2_HOME/bowtie2-build /home/kate/Dokumenty/FAV/Diplomka/existujicisw/bowtie2-2.4.1-linux-x86_64/example/reference/KIR2DL1_gen.fasta kir

# pak to prohnat cyklem
#$BT2_HOME/bowtie2 -x kir -1 /home/kate/Dokumenty/FAV/Diplomka/existujicisw/mojedata/moje_art_data1.fq -2 /home/kate/Dokumenty/FAV/Diplomka/existujicisw/mojedata/moje_art_data2.fq -S /home/kate/Dokumenty/FAV/Diplomka/existujicisw/mojedata/result.sam

# vedle se mužeš vygenerovat genom z artu a pak to nahaziš do bowtie a rovnou mu podsuneš referenční gen to znamená nějakej gen z KIR
# musim vzít referenční genomi a na nich zkusit porovnat jestli se to náhodou nějakýmu nehodí

# určení fitnes funkce
# dá se tam zakomponovat že v některcýh případech nějaké genomi spolu mají větší pravděpodobnost než jiné
# specializuješ se jen na KIR geny