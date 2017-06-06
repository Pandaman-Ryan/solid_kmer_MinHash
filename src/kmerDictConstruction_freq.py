'''
This file contains a class which call DSK to find all the frequent Kmers in genomes
and then build a dictionary to save them.
@ An Zheng
'''
# module
import sys
import time
import os
import copy
import numpy
from Bio.Seq import Seq

# self-defined module
class dictConstructor(object):
	def __init__(self):
		print ("Start the construction of dictionary...")

	# public functions
	def constructKmerDict(self, kmerSize, genomePath, h5FilePath, kmerFilePath,\
				kmerFreq, smartFreq, genomeLength):
		# construct a trie of kmer
		kmerPath = self.__callDSKtofindKmer(kmerSize, kmerFreq, genomePath, h5FilePath, kmerFilePath)

		# read kmer
		kmerSet = self.__readKmerFromFile(kmerPath, smartFreq, genomeLength)
		kmerDict = self.__constructDict(kmerSet)
		print ("Trie construction complete!")
		return kmerDict

	# private functions
	def __callDSKtofindKmer(self, kmerSize, kmerFreq, genomePath, h5FilePath, kmerFilePath):
		# call DSK program to find kmer
		genomeFile = genomePath
		tempH5File = h5FilePath
		kmerFile = kmerFilePath
		
		print ("Using DSK to process genomes...")
                # /home/ryan/Code/replicatedSegment/Version6/DSK/dsk-2.0.7-Linux/bin/dsk -verbose 0 -file "input" -kmer-size 16 -abundance-min 1 -abundance-max 50 -out outputfile
		os.system("./DSK/dsk-2.0.7-Linux/bin/dsk -verbose 0 \
			-file "+ genomeFile + " -kmer-size "+str(kmerSize)+\
			" -abundance-min "+str(kmerFreq)+" -abundance-max 50 "+\
			"-out "+ tempH5File)
		os.system("./DSK/dsk-2.0.7-Linux/bin/dsk2ascii -verbose 0 \
			-file "+tempH5File+\
			" -out " + kmerFile)

		# return the path where kmer list are saved
		return kmerFile

	def __readKmerFromFile(self, kmerPath, smartFreq, genomeLength):
		# read kmer
		kmerSet = set()
		if not smartFreq:
			ifile = open(kmerPath)
			for line in ifile:
				kmer, freq = line.split()
				kmerSet.add(kmer)
			ifile.close()
		else:
			ifile = open(kmerPath)
			freqArr = []
			for line in ifile:
				kmer, freq = line.split()
				freqArr.append(int(freq))
			ifile.close()

			freqArr.sort(reverse=True)
			thresFreq = 2
			if genomeLength < len(freqArr):
				thresFreq = min(2,freqArr[genomeLength])
			print ("You are using the option \"--smartFreq\", the frequence threshold is: ",
				thresFreq)

			ifile = open(kmerPath)
			for line in ifile:
				kmer, freq = line.split()
				if int(freq) >= thresFreq:
					kmerSet.add(kmer)
			ifile.close()
		# END IF
		return kmerSet

	def __constructDict(self, kmerSet):
		# construct Trie
		print ("Start to construct Dict...")
		print ("There are ",len(kmerSet)," to add..")
		kmerCount = 0
		kmerDict = {}
		for index, kw in enumerate(kmerSet):
			kmerCount += 1
			if kmerCount%1000000 == 0:
				print ("Now processing ", kmerCount," -th kmer")
			r_kw = self.__computeReverseComplement(kw)
			kmerDict[kw] = index
			kmerDict[r_kw] = index

		print ("Dict constructed")
		return kmerDict

	def __computeReverseComplement(self, seq_toprocess):
		f_seq = Seq(seq_toprocess)
		r_seq = f_seq.reverse_complement()
		return str(r_seq)

###END#####################################
