'''
This file contain a class which deals with segment pair filtering 
@ An Zheng
'''

import random
import copy

class pairFilter(object):
	def __init__(self):
		print ("Start to filter the segment pairs...")

	def filterPair(self, segmentPairs, indexedReads, threshold, jarccardTable):
		segmentPairs_filtered = []
		print ("There are ",len(segmentPairs)," pairs to 2nd check...")
		count = 0
		for key1, key2 in segmentPairs:
			# order filter
			# only keep the pair which chrIndex1 <= chrIndex2
			# if chrIndex1 == chrIndex2 then compare segIndex
			count += 1
			if count%100000 == 0:
				print ("Now processing: ",count,"-th")

			index1 = int(key1[1:])
			index2 = int(key2[1:])
			if index1 < index2:
				if jarccardTable != None and key1 in jarccardTable and key2 in jarccardTable[key1]:
					jarccard = jarccardTable[key1][key2]
				else:
					seg1 = indexedReads[key1]
					seg2 = indexedReads[key2]
					jarccard = self.__checkJarccard(seg1, seg2)
				if jarccard >= threshold:
					segmentPairs_filtered.append((key1, key2))

		# segmentPairs_matched = self.__findMatchedPairs(segmentPairs_filtered)
		segmentPairs_matched = segmentPairs_filtered

		print ("There are ", len(segmentPairs_matched)," left.")
		return segmentPairs_matched

	def __checkOverlap(self, sample, listToCheck):
		sample = copy.deepcopy(sample)
		listToCheck = copy.deepcopy(listToCheck)

		sample.sort()
		listToCheck.sort()
		
		iter_sample = 0
		iter_list = 0
		while iter_sample < len(sample) and iter_list < len(listToCheck):
			if sample[iter_sample] == listToCheck[iter_list]:
				return True
				break
			elif sample[iter_sample] < listToCheck[iter_list]:
				iter_sample += 1
			else:
				iter_list += 1
		
		return False

	def __checkJarccard(self, sample, listToCheck):
		sample = copy.deepcopy(sample)
		listToCheck = copy.deepcopy(listToCheck)

		sample.sort()
		listToCheck.sort()
		
		iter_sample = 0
		iter_list = 0
		count_overlap = 0
		while iter_sample < len(sample) and iter_list < len(listToCheck):
			if sample[iter_sample] == listToCheck[iter_list]:
				count_overlap += 1
				iter_sample += 1
				iter_list += 1
			elif sample[iter_sample] < listToCheck[iter_list]:
				iter_sample += 1
			else:
				iter_list += 1

		return float(count_overlap)/float(min(len(sample), len(listToCheck)))
###END#######################################
