'''
Vesion4: minHash powered
Do All file
@ An Zheng

-----test version----
Do all the steps below
In this version, we only include chromosome 1 and 2.
1: read choromosome from file.
2: find k-mer in the chormosome using DSK
3: construct an aho corasick trie
4: divide each genomes into N-length segments and index the frequent kmers
5: compute sketch for every segments
6: find pairs
7: evaluate pairs
8: extend the pairs found
'''
# module
import sys
import time
import os
import copy
import itertools

# self-defined module
from .dupDetectorFile import dupDetector
# import evaluate

class overlapDectector(object):
	def detectOverlaps(self, genomePath, kmerSize, h5FilePath, kmerFilePath,\
			specificity, minTimes, maxTimes,\
			minThreshold, maxThreshold, gap,\
			jac_threshold, hash_threshold,\
			refPath, minOverlap,\
			savePath, freqKmer, kmerFreq, smartFreq, testing):	
		startTime = time.time()
		# initialization: set parameters, load the genome data, construct Kmer-Trie
		dd = dupDetector(genomePath, kmerSize)

		# segment the genomes
		indexedReads = dd.segmentAndIndexGenome(genomePath, h5FilePath, kmerFilePath, freqKmer, kmerFreq, smartFreq)

		print ("The whole program has used: "+str(time.time()-startTime)+" seconds")
		if testing:
			thresholdBasis = minThreshold/maxTimes
			jaccardTable, minHashArr = dd.createJaccardTable(indexedReads, thresholdBasis, specificity)
			print ("The whole program has used: "+str(time.time()-startTime)+" seconds")

			for JHTimes in range(minTimes, maxTimes+1):
				ready_to_break = False
				for specificThreshold in [minThreshold+increase*gap for increase in range(int((maxThreshold-minThreshold)/gap+1))]:
					threshold = specificThreshold/JHTimes
					print (">>>>>>>>>>>>")
					print ("hThres: "+str(threshold)+"\t jThres: "+str(specificThreshold))
					# perform min-hash on the genomes
					overlappedReads = dd.generatePossiblePair_loadCache(indexedReads, threshold, specificity, minHashArr)			
					# filter the unlikely pairs using the kmer-index alignment
					readPairs_filtered = dd.filterSegPair_loadCache(overlappedReads, indexedReads, specificThreshold, jaccardTable)
					# calculate the precision and recall rate
					precision, recall = dd.evaluate(readPairs_filtered, refPath, minOverlap, savePath, threshold, specificThreshold)
					print ("The whole program has used: "+str(time.time()-startTime)+" seconds")
					
					# when precision and recall becomes zero, break in advance
					if ready_to_break and (precision == 0  and recall == 0):
						break
					elif precision == 0  and recall == 0:
						ready_to_break = True
					
			# END FOR
		else:
			threshold = hash_threshold
			specificThreshold = jac_threshold
			print (">>>>>>>>>>>>")
			print ("hThres: "+str(threshold)+"\t jThres: "+str(specificThreshold))
			minHashArr = None
			jaccardTable = None
			overlappedReads = dd.generatePossiblePair_loadCache(indexedReads, threshold, specificity, minHashArr)			
			readPairs_filtered = dd.filterSegPair_loadCache(overlappedReads, indexedReads, specificThreshold, jaccardTable)
			precision, recall = dd.evaluate(readPairs_filtered, refPath, minOverlap, savePath, threshold, specificThreshold)
			print ("The whole program has used: "+str(time.time()-startTime)+" seconds")
		# END IF
		
		# program done, the program statistics
		print (">>>>>>>>>>>>")
		endTime = time.time()
		print ("The program works for "+str(endTime-startTime)+" seconds in total.")
		print ("Done")

####END####################################
