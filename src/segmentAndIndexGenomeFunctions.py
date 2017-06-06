'''
Version3: minHash powered

In this class, we implement several functions below:


'''

import math
import numpy
# self defined class
#from kmerTrieConstruction import trieConstructor
from .kmerDictConstruction_freq import dictConstructor as dictConstructor_freq
from .kmerDictConstruction_all import dictConstructor as dictConstructor_all

class si_tools(object):
	####################
	# public functions
	# private functions
	def indexGenome(self, genome, kmerSize, genomePath, h5FilePath, kmerFilePath,\
			freqKmer, kmerFreq, smartFreq, genomeLength):

		if freqKmer:
			dc = dictConstructor_freq()
			kmerDict = dc.constructKmerDict(kmerSize, genomePath, h5FilePath, kmerFilePath,\
					kmerFreq, smartFreq, genomeLength)
		else:
			dc = dictConstructor_all()
			kmerDict = dc.constructKmerDict(kmerSize, genome)

		print ("Start indexing the genome segments. This process may take a couple hours...")
		indexedReads = {}
		for readIndex in genome:
			read = genome[readIndex]
			readSize = len(read)
			indexedSeg = self.__indexGenomeWithKmer(read, kmerDict, kmerSize, readSize)
			indexedReads[readIndex] = indexedSeg

		return indexedReads

	# private functions
	def __indexGenomeWithKmer(self, genome, kmerDict, kmerSize, readSize):
		indexedGenome = []
		for startPos in range(readSize-kmerSize):
			kmer = genome[startPos:startPos+kmerSize]
			if kmer in kmerDict:
				indexedGenome.append(kmerDict[kmer])

		return indexedGenome		
#########END###################
