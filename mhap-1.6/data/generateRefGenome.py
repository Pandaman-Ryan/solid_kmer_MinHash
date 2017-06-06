'''
Generate Reference Genome
@ An Zheng
This python code can read the genome sequence and cut a given-length segment from it.
The segment contains no "N"
'''

# module
import sys
import time
import os
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def main():
	genomePath = "../../../data/chromFa/hg19.fa"
	selectedChromesome = "chr21"
	segmentLength = 1000000
	savePath = "test_sample.fa"
	'''
	readGenome(genomePath, selectedChromesome, segmentLength, savePath)
	'''
	randomGenerator(segmentLength, savePath)

	print ("Sampling complete!")


def randomGenerator(segmentLength, savePath):
	genome = ""
	for index in range(segmentLength):
		nuc = random.choice(["A","T","C","G"])
		genome += nuc

	saveGenomeInFile(genome, savePath)

def readGenome(genomePath, selectedChromesome, segmentLength, savePath):
	genomeDict = {}

	for seq in SeqIO.parse(genomePath, "fasta"):
		chromID = seq.id
		genomeDict[chromID] = seq

	seqToSave = str((genomeDict[selectedChromesome]).seq)
	seqToSave_upper = seqToSave.upper()

	while True:
		startPos = random.randint(0, len(seqToSave_upper)-segmentLength)
		genome = seqToSave_upper[startPos: startPos+segmentLength]

		if "N" not in genome:
			break

	saveGenomeInFile(genome, savePath)

def saveGenomeInFile(genome, savePath):
	recordToSave = SeqRecord(Seq(genome), id = "chr", name = "chr", description = "chr")
	ofile = open(savePath, 'w')
	SeqIO.write(recordToSave, ofile, "fasta")
	ofile.close()

main()
