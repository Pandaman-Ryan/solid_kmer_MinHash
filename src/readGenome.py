'''
This file contains functions of reading genome from fa files
@ An Zheng

-----test version----
In this version, we only include chromosome 1 and 2.
-> read choromosome from file.i

'''
# module
import sys
import time
import os
from Bio import SeqIO

class genomeLoader_noMask(object):
	def __init__(self):
		print ("Start loading genome sequences...")

	# public functions
	def readGenome(self, genomePath):
		# read genomes
		genomeDict = {}
		startPosition = {}
		genomeLength = 0
		for seq in SeqIO.parse(genomePath, "fasta"):
			chromID = seq.id
			chromDes = seq.description

			seqToSave = str(seq.seq)
			seqToSave_upper = seqToSave.upper()
			genomeDict[chromID] = str(seqToSave_upper)

			chromStart = int((chromDes.split())[-1])
			startPosition[chromID] = chromStart

			genomeLength += len(seqToSave_upper)
		# END FOR
		#self.__saveGenomeInMultiFile(genomeDict, chromToLoad, genomeSaveDir, fileType)
		print ("number of reads: " + str(len(genomeDict)))
		return genomeDict, startPosition, genomeLength

	'''
	def __saveGenomeInMultiFile(self, genomeDict, chromToLoad, genomeSaveDir, fileType):
		genomeArr = []
		for chromName in chromToLoad:
			genomeArr.append(genomeDict[chromName])

		savePath = genomeSaveDir+"hg19_sample"+fileType
		ofile = open(savePath, 'w')
		SeqIO.write(genomeArr, ofile, "fasta")
		ofile.close()
	'''

###END#################################
