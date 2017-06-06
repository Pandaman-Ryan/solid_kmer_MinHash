import os

os.system("./dsk -verbose 0 -file e_coli_genome.fasta -kmer-size 15 -abundance-min 10 -abundance-max 20 -out e_coli_genome.h5")
os.system("./dsk2ascii -verbose 0 -file e_coli_genome.h5 -out e_coli_genome.txt")
