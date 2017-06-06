'''
Version3: minHash powered

In this class, we implement several functions below:


'''

import math
import numpy
# self defined class
from .jacTable import jaccardTableGenerator
from .readGenome import genomeLoader_noMask as genomeLoader
from .segmentAndIndexGenomeFunctions import si_tools
from .generateSketch import sketchGenerator_v1 as sketchGenerator
from .filterSegmentPairs import pairFilter

class dupDetector(object):
	# in init function, we set the parameters, load the genomes, and construct kmer-trie of the genome
	def __init__(self, genomePath, kmerSize):
		self.__kmerSize = kmerSize	
		# loading genome sequence
		gl = genomeLoader()
		self.genomes, self.startPositions, self.genomeLength = gl.readGenome(genomePath)
		print ("Duplication detector initialization complete!")

	####################
	# public functions
	def segmentAndIndexGenome(self, genomePath, h5FilePath, kmerFilePath, freqKmer, kmerFreq, smartFreq):
		############################
		# index the genome according to the info in trie
		# construct kmer-dict
		si = si_tools()
		genomes = self.genomes
		kmerSize = self.__kmerSize
		genomeLength = self.genomeLength
		indexedReads = si.indexGenome(genomes, kmerSize, genomePath, h5FilePath, kmerFilePath,\
					 freqKmer, kmerFreq, smartFreq, genomeLength)

		print ("Genome segments are indexed!")
		return indexedReads
		
	def createJaccardTable(self, indexedReads, threshold, specificity):
		jt = jaccardTableGenerator()
		jaccardTable, minHashArr = jt.generateTable(indexedReads, threshold, specificity)
		return jaccardTable, minHashArr
	
	def generatePossiblePair_loadCache(self, indexedReads, threshold, specificity, minHashArr):
		sk = sketchGenerator()
		overlappedReads, oldHashArr = sk.findPairWithMinHash(indexedReads, threshold, specificity, minHashArr)
		return overlappedReads
	
	def filterSegPair_loadCache(self, segmentPairs, indexedReads, threshold, jaccardTable):
		pf = pairFilter()
		segmentPairs_filtered = pf.filterPair(segmentPairs, indexedReads, threshold, jaccardTable)
		return segmentPairs_filtered

	def evaluate(self, readPairs, refPath, minOverlap, savePath, hThres, jThres):
		refFile = open(refPath)
		readPositions = {}
		#readSize = parameters.readSize
		for line in refFile:
			elems = line.split()
			name = elems[0]
			readStart = int(elems[5])
			readEnd = int(elems[6])
			if readEnd - readStart > minOverlap:#readSize*0.1:
				if name in readPositions:
					(readInRef_start_old, readInRef_end_old) \
							= readPositions[name]
					readInRef_start_new = int(elems[9])
					readInRef_end_new = int(elems[10])
					if (readInRef_end_new-readInRef_start_new) > \
						(readInRef_end_old-readInRef_start_old):
						readPositions[name] = (readInRef_start_new,\
									readInRef_end_new)
				else:
					readInRef = (int(elems[9]),int(elems[10]))
					readPositions[name] = readInRef

		countTruePairs = 0
		for key1 in readPositions:
			for key2 in readPositions:
				if key1 != key2 and int(key1[1:])< int(key2[1:]):
					ovlpRange = self.__checkOverlap(readPositions[key1],\
							 readPositions[key2])
					if ovlpRange > minOverlap:#readSize * 0.1:
						countTruePairs += 1
					

		countOvlp = 0
		for (key1, key2) in readPairs:
			if key1 in readPositions and key2 in readPositions:
				ovlpRange = self.__checkOverlap(readPositions[key1], readPositions[key2])
				if ovlpRange > minOverlap:#readSize * 0.1:
					countOvlp += 1
		print ("Find " + str (len(readPairs)) +" overlapped reads, in which there are "+\
			str(countOvlp) + " is correct. And there should be "+str(countTruePairs)+\
			" pairs in all.")

		if len(readPairs) == 0:
			precision = 0
			print("The precision is: 0")
		else:
			precision = countOvlp/float(len(readPairs))
			print ("The precision is: "+str(precision))

		if countTruePairs == 0:
			recall = 0
			print("The recall is: 0")
		else:
			recall = countOvlp/float(countTruePairs)
			print ("The recall is: " + str(recall))

		if savePath != None:
			saveFile = open(savePath, 'a')
			saveFile.write(str(hThres)+'\t'+str(jThres)+'\t'+str(precision)+'\t'+str(recall)+'\n')
			saveFile.close()

		return precision, recall

	# private functions
	def __checkOverlap(self, range1, range2):
		start1, end1 = range1
		start2, end2 = range2

		if end1<start2 or end2<start1:
			return 0
		else:
			ovlpRange = min(end1,end2) - max(start1, start2)
			return ovlpRange

	def __computeError(self, seq1, seq2):
		return 0, 0, 0	

	def __saveInfoIndex(self, dataGenome, path):
		# the format is [[[1,2,3],[4,5,6]], [[7,8,9],[10,11,12]], ...]
		# save format is 1|2|3 4|5|6
		#		 7|8|9 10|11|12
		ofile = open(path,'w')
		for dataChr in dataGenome:
			for dt in dataChr:
				saveStr = ""
				if len(dt) != 0:
					for index in dt:
						saveStr += str(index)+"|"
					saveStr = saveStr[:-1]
				ofile.write(saveStr+",")
			ofile.write("\n")

		ofile.close()

	def __readInfoIndex(self, path):
		dataGenome = []

		ifile = open(path)
		for line in ifile:
			dataChr = []
			# elemArrayList = ["101|102|103", "104","105"]
			elemArrayList = line.split(",")
			# there are an empty one at the end
			# (because of the drawback in the saveInfoIndex function)
			elemArrayList.pop()
			for elemArr in elemArrayList:
				data = []
				if len(elemArr) != 0:
					# elems = ["101","102","103"]
					elems = elemArr.split("|")
					for em in elems:
						data.append(int(em))
				dataChr.append(data)
			dataGenome.append(dataChr)

		self.indexedGenomeSeg = dataGenome
		return dataGenome

	def __saveSegmentPairs(self, segmentPairs, path):
		# pair: ((chromIndex1, pos1), (chromIndex2, pos2))
		ofile = open(path, 'w')
		for pair in segmentPairs:
			((chromIndex1, segIndex1), (chromIndex2, segIndex2)) = pair
			strToSave = str(chromIndex1)+","+str(segIndex1)+"\t"+\
					str(chromIndex2)+","+str(segIndex2)+"\n"
			ofile.write(strToSave)

		ofile.close()
	
		
#########END###################
