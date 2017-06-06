'''
This file contains a class which use minHash method to generate sketch
for each indexed segments.
@An Zheng

-------Test Version 1------------
-> use minHash and minHashLSH function in datasketch lib
* run out of cache when running it

-------Test Version 2------------
-> use minHash only and "manually" match the pairs

'''
import sys
import time
# version 1
from datasketch import MinHash, MinHashLSH


######Version 1#########################################
class sketchGenerator_v1(object):
	def __init__(self):
		print ("Start to generate sketches...")

	def findPairWithMinHash(self, indexedReads, threshold, specificity, minHashArr_prevComp):
		minHashArr, lsh = self.__buildLSH(indexedReads, threshold, specificity, minHashArr_prevComp)
		readPairs = self.__pairSegment(minHashArr, lsh)
		print ("Find all sketch pairs! We found "+str(len(readPairs))+" pairs in all.")
		return readPairs, minHashArr


	# private functions
	def __buildLSH(self, indexedReads, thres, specificity, minHashArr_prevComp):
		if minHashArr_prevComp == None:
			minHashArr = []
			for readIndex in indexedReads:
				if len(indexedReads[readIndex]) != 0:
					m = MinHash(num_perm = specificity)
					for kmerIndex in indexedReads[readIndex]:
						kmerIndex_str = str(kmerIndex)
						m.update(kmerIndex_str.encode("utf8"))
					minHashFeat = (readIndex,m)
				else:
					minHashFeat = (readIndex,"NULL")
				minHashArr.append(minHashFeat)
			print ("Now start to build the LSH. There are "+str(len(minHashArr))+ " segments in all...")
		else:
			minHashArr = minHashArr_prevComp

		lsh = MinHashLSH(threshold=thres, num_perm = specificity)
		for key, minH in minHashArr:
			if type(minH) != type("NULL"):
				lsh.insert(key, minH)
		print ("LSH has been constructed!")

		return minHashArr, lsh

	def __pairSegment(self, minHashArr, lsh):
		segPairs = []
		print ("Start to pair segments...")
		print ("There are "+str(len(minHashArr))+ " segments" )
		for key, minH in minHashArr:
			if type(minH) != type("NULL"):
				candidates = lsh.query(minH)
				for keyOfCandidate in candidates:
					segPairs.append((key, keyOfCandidate))

		return segPairs
###END##########################
