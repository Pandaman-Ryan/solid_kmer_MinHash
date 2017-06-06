'''
This file contains a class which reads all Kmers in genomes,
convert them into integers,
and then use a dictionary to save them
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
	def constructKmerDict(self, kmerSize, genome):
		# construct a trie of kmer
		#kmerPath = self.__callDSKtofindKmer(kmerSize, genomePath, h5FilePath, kmerFilePath)

		# read kmer
		kmerSet = self.__readKmerFromFile(genome, kmerSize)
		kmerDict = self.__constructDict(kmerSet)
		print ("Trie construction complete!")
		return kmerDict

	# private functions
	def __readKmerFromFile(self, genome, kmerSize):
		# read kmer
		kmerSet = set()
		for readIndex in genome:
			read = genome[readIndex]
			readSize = len(read)
			for startPos in range(readSize-kmerSize):
				kmer = genome[readIndex][startPos:startPos+kmerSize]
				kmerSet.add(kmer)
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
