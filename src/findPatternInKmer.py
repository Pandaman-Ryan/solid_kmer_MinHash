
'''
This file contains functions of reading genome from file and find the kmer using DSK
@ An Zheng

-----test version----
In this version, we only include chromosome 1 and 2.
1: read choromosome from file.
2: find 20-mer in the chormosome using DSK
3: construct an aho corasick trie
-> 4: find rep segments in index array

'''
# module
import sys
import time
import os
import copy
import itertools

# self-defined module
import readGenome
import findAndIndexKmer
import parameters
import extendKmer
import evaluate

def readInfoIndex(path):
	dataGenome = []

	ifile = open(path)
	for line in ifile:
		dataChr = []
		# elemArrayList = ["101|102|103", "104","105"]
		elemArrayList = line.split()
		for elemArr in elemArrayList:
			data = []
			# elems = ["101","102","103"]
			elems = elemArr.split("|")
			for em in elems:
				data.append(int(em))
			dataChr.append(data)
		dataGenome.append(dataChr)

	return dataGenome

def findRepFromIndex(indexedGenome, indexKmerSize):
	#i record appearing positions for each kmer
	repIndexKmer = set()
	kmerDict = {}
	for chromIndex in range(len(indexedGenome)):
		indexedChrom = indexedGenome[chromIndex]
		for startPos in range(len(indexedChrom)-(indexKmerSize-1)):
			indexKmerList = enumerateIndexKmer(indexedChrom, startPos, indexKmerSize)
			# indexKmer = tuple(indexedChrom[startPos:startPos+indexKmerSize])
			for indexKmer in indexKmerList:
				if indexKmer in kmerDict:
					repIndexKmer.add(indexKmer)
					kmerDict[indexKmer].append((startPos,chromIndex))
				else:
					kmerDict[indexKmer] = [(startPos,chromIndex)]

	print ("find "+str(len(repIndexKmer))+" segments in all.")

	# extend from seed
	count = 0
	segmentPairSet = set()
	for kmer in repIndexKmer:
		count += 1
		if count% 100000 == 0:
			print(count)
			
		posList = kmerDict[kmer]
		lenOfSeed = len(kmer)
		for i,j in itertools.combinations(range(len(posList)), 2):
			startPos1, chromIndex1 = posList[i]
			startPos2, chromIndex2 = posList[j]
			# extend:
			# output: a tuple of a sequence of index which corresponds the index
			# segment = (startPos, endPos, chromsomeIndex)
			segment1, segment2 = extendKmer.extend(indexedGenome, lenOfSeed,\
						startPos1, chromIndex1, startPos2, chromIndex2)
			###prevent intersections... combine two intersections
			if ((segment1,segment2) not in segmentPairSet) and\
				((segment2,segment1) not in segmentPairSet):
				startPos1 = segment1[0]
				endPos1 = segment1[1]
				chormIndex2 = segment1[2]
				startPos2 = segment2[0]
				endPos2 = segment2[1]
				chormIndex2 = segment2[2]
				segmentPair = copy.deepcopy((segment1,segment2))
				segmentPairSet.add(segmentPair)

	print ("find "+str(len(segmentPairSet))+" segment pairs in all")

	return segmentPairSet		
			
def enumerateIndexKmer(indexedChrom, startPos, indexKmerSize):
	# indexKmer = tuple(indexedChrom[startPos:startPos+indexKmerSize])
	currentPos = 0
	kmerList = recursionEnum(indexedChrom, startPos, indexKmerSize, currentPos)

	return kmerList

def recursionEnum(indexedChrom, startPos, indexKmerSize, currentPos):
	kmerListRenew = []
	if currentPos == indexKmerSize-1:
		currentChoices = indexedChrom[startPos+currentPos]
		for index in currentChoices:
			newKmer = (index,)
			kmerListRenew.append(newKmer)
		return	kmerListRenew
	elif currentPos >= indexKmerSize:
		print ("Erro of currentPos @ func recursionEnum")
		sys.exit(1)
		return

	nextPos = currentPos + 1
	kmerList = recursionEnum(indexedChrom, startPos, indexKmerSize, nextPos)

	currentChoices = indexedChrom[startPos+currentPos]
	for kmer in kmerList:
		for index in currentChoices:
			newKmer = (index,)+kmer
			kmerListRenew.append(newKmer)
	
	return kmerListRenew


def convertIndexToSeq(indexSegmentPairSet, indexPositionPath):
	indexPosition = readInfoPos(indexPositionPath)
	
	genomeSegmentPairSet = []
	for pair in indexSegmentPairSet:
		startPos1 = pair[0][0]
		endPos1 = pair[0][1]
		chromIndex1 = pair[0][2]
		startPos2 = pair[1][0]
		endPos2 = pair[1][1]
		chromIndex2 = pair[1][2]

		#print (indexedGenome[chromIndex1][startPos1:endPos1])
		#print (indexedGenome[chromIndex2][startPos2:endPos2])
		startInGenome1 = indexPosition[chromIndex1][startPos1]
		endInGenome1 = indexPosition[chromIndex1][endPos1-1]
		startInGenome2 = indexPosition[chromIndex2][startPos2]
		endInGenome2 = indexPosition[chromIndex2][endPos2-1]


		segment1 = (startInGenome1, endInGenome1, chromIndex1)
		segment2 = (startInGenome2, endInGenome2, chromIndex2)
		genomeSegmentPairSet.append((segment1, segment2))

	return genomeSegmentPairSet
	
def readInfoPos(path):
	dataGenome = []

	ifile = open(path)
	for line in ifile:
		dataChr = []
		elemArrayList = line.split()
		for elem in elemArrayList:
			dataChr.append(int(elem))
		dataGenome.append(dataChr)

	return dataGenome
	
def savePairSet(genomeSegmentPairSet, segPairPath):
	ofile = open(segPairPath, 'w')
	for pair in genomeSegmentPairSet:
		startPos1 = pair[0][0]
		endPos1 = pair[0][1]
		chromIndex1 = pair[0][2]
		startPos2 = pair[1][0]
		endPos2 = pair[1][1]
		chromIndex2 = pair[1][2]
		strToWrite = str(startPos1)+" "+str(endPos1)+" "+str(chromIndex1)+" "+\
				str(startPos2)+" "+str(endPos2)+" "+str(chromIndex2)
		ofile.write(strToWrite+"\n")

	ofile.close()

	
def readPairSet(segPairPath):
	ifile = open(segPairPath)
	genomeSegmentPairSet = []
	for line in ifile:
		elems = line.split()
		startPos1 = int(elems[0])
		endPos1 = int(elems[1])
		chromIndex1 = int(elems[2])
		startPos2 = int(elems[3])
		endPos2 = int(elems[4])
		chromIndex2 = int(elems[5])
		genomeSegmentPairSet.append( ((startPos1, endPos1, chromIndex1),\
					(startPos2, endPos2, chromIndex2)) )

	return genomeSegmentPairSet

def sortPairs(pairSet):
	

	return pairSet

if __name__ == "__main__":
	startTime = time.time()

	# parameters
	numOfGenome = parameters.numOfGenome
	genomeFileDir = parameters.genomeFileDir
	minKmerSize = parameters.minKmerSize
	indexKmerSize = parameters.indexKmerSize
	indexedGenomePath = parameters.indexedGenomePath
	indexPositionPath = parameters.indexPositionPath
	fileType = parameters.fileType
	segPairPath = parameters.segPairPath
	trueSetPath = parameters.trueSetPath
	overlapThreshold = parameters.overlapThreshold

	'''
	# read genome
	genome = readGenome.readGenome(genomeFileDir, numOfGenome, fileType)

	# index the genome sequences using kmer trie
	trie = constructTrieFromKmerFile(numOfGenome, minKmerSize, genomeFileDir)
	indexedGenome, indexPosition = indexGenomeWithKmer(genome, trie, minKmerSize)
	'''

	'''
	# find replicated index in the indexedGenome
	indexedGenome = readInfoIndex(indexedGenomePath)
	indexSegmentPairSet = findRepFromIndex(indexedGenome, indexKmerSize)
	genomeSegmentPairSet = convertIndexToSeq(indexSegmentPairSet, indexPositionPath)
	savePairSet(genomeSegmentPairSet, segPairPath)
	'''
	
	# evaluating
	genomeSegmentPairSet = readPairSet(segPairPath)
	genomeSegmentPairSet_sorted = sortPairs(genomeSegmentPairSet)
	precision, recall = evaluate.evaluate(trueSetPath, genomeSegmentPairSet_sorted, overlapThreshold)
	print ("Precision rate is: "+str(precision))
	print ("Recall rate is: "+str(recall))
	
	# program done, the program statistics
	endTime = time.time()
	print ("The program works for "+str(endTime-startTime)+" in total.")
	print ("Done")


