'''
Before this, you need:
export FMH=$(pwd)
Then run this program
'''

from optparse import OptionParser
import sys
import os
import subprocess
import random
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# from source
#from src_freq_kmer.doall import overlapDectector
#from src_all_kmer.doall import overlapDectector as overlapDectector_cmp
from src.doall import overlapDectector

def main():
	print ("====START=======================================")
	os.environ['FMH'] = os.getcwd()
	##########################
	# read options and args
	###########################
	usage = "usage: %prog [options] <output_file_path>"
	parser = OptionParser(usage)
	# index of this experiment
	parser.add_option('--experiment', dest='exp_index', default='0', help='The index of experiment [Default: %default]')

	# options for loading reference genome
	parser.add_option('--randSource', dest='randRef', default=False,  action='store_true', help='select a randomly generated sequence as the reference [Default: %default]')
	parser.add_option('--lenRef', dest='len_ref', default=1000000, type='int', help='Length of the reference genome [Default: %default]')
	parser.add_option('--ref', dest='ref_path', default=None, help="path of reference genome from file")

	# options for generating reads
	parser.add_option('--read', dest='read_path', default=None, help="path of reads from file")
	parser.add_option('--randRead', dest='randRead', default=False,  action='store_true', help='randomly generate reads [Default: %default]')
	parser.add_option('--nRead', dest='num_reads', default=500, type='int', help='Number of reads [Default: %default]')
	parser.add_option('--lenRead', dest='len_read', default=10000, type='int', help='Length of each read [Default: %default]')
	parser.add_option('--mut', dest='mut_rate', default=0.05, type='float', help='Mutation rate [Default: %default]')
	parser.add_option('--ins', dest='ins_rate', default=0.01, type='float', help='Insertion rate [Default: %default]')
	parser.add_option('--del', dest='del_rate', default=0.01, type='float', help='Deletion rate [Default: %default]')

	# options for true set
	parser.add_option('--truthPath', dest='true_set', default=None, help='true set path [Default: %default]')
	parser.add_option('--minReadLen', dest='minReadLen', default=200, type='int', help='the minimal read length setting used in Blasr [Default: %default]')
	parser.add_option('--nproc', dest='nproc', default=8, type='int', help='the number of proc used in Blasr [Default: %default]')
	parser.add_option('--bestn', dest='bestn', default=10, type='int', help='the bestn setting in Blasr [Default: %default]')

	# options for overlapping finding
	parser.add_option('--useFreqKmer', dest='freqKmer', default=False,  action='store_true', help='use freqent kmer to find the overlapped read pairs [Default: %default]')
	parser.add_option('--kmerFreq', dest='kmerFreq', default=2,  type='int', help='the minimal frequence of solid k-mers [Default: %default]')
	parser.add_option('--smartFreq', dest='smart', default=False,  action='store_true', help='smartly select kmer frequence [Default: %default]')
	parser.add_option('--kmerSize', dest='kmer_size', default=16, type='int', help='The size of frequent kmer (input of DSK) [Default: %default]')
	parser.add_option('--speci', dest='specificity', default=512, type='int', help='Specificity in minHash [Default: %default]')
	parser.add_option('--h5', dest='h5_path', default=None, help='The path which save h5 file generated by DSK [Default: %default]')
	parser.add_option('--kmer', dest='kmer_path', default=None, help='The path which save kmer file generated by DSK [Default: %default]')
	parser.add_option('--minOvlp', dest='min_ovlp', default=2000, type='float', help='the minimal overlaps of two reads which are considered as a qualified overlapped pair [Default: %default]')

	# is testing or regularly running
	parser.add_option('--testing', dest='testing', default=False,  action='store_true', help='Is testing or not [Default: %default]')

	# IF RUNNING NORMALLY
        # settings for filters' thresholds
	parser.add_option('--hThres', dest='hash_threshold', default=0.03, type='float', help='The Jaccard threshold for the 1nd filter [Default: %default]')
	parser.add_option('--jThres', dest='jac_threshold', default=0.1, type='float', help='The Jaccard threshold for the 2nd filter [Default: %default]')
	
	# IF TESTING
        # settings for filters' thresholds
	parser.add_option('--minjThres', dest='min_jac_threshold', default=0.01, type='float', help='The minimal Jaccard threshold for the 2nd filter [Default: %default]')
	parser.add_option('--jGap', dest='thres_gap', default=0.01, type='float', help='The gap of jThres increasing for setting different testing thresholds [Default: %default]')
	parser.add_option('--maxjThres', dest='max_jac_threshold', default=0.10, type='float', help='The maximal Jaccard threshold for the 2nd filter [Default: %default]')
	parser.add_option('--minTimes', dest='min_times', default = 5, type='int', help='the minimal times of jThres than hThres [Default: %default]')
	parser.add_option('--maxTimes', dest='max_times', default = 5, type='int', help='the maximal times of jThres than hThres [Default: %default]')

	# save file
	parser.add_option('--saveFreq', dest='save_path_freq', default=None, help="If you want to save the precision and recall, set this option with a path")
	parser.add_option('--saveAll', dest='save_path_all', default=None, help="If you want to save the precision and recall, set this option with a path")

	(options, args) = parser.parse_args()

	############################
	# generate reference genome
	############################
	random.seed()
	print (">>>>>>>>>>>>")
	print ("Loading the reference genome...")
	if options.randRef:
		savePath = ("%s/data/refGenome_rand_len"+str(options.len_ref)+"_experi"+options.exp_index+".fa") \
						%os.environ['FMH']
		genome = ""
		for index in range(options.len_ref):
			nuc = random.choice(["A","T","C","G"])
			genome += nuc

		recordToSave = SeqRecord(Seq(genome), id = "chr", name = "chr",\
					description = "chr")
		ofile = open(savePath, 'w')
		SeqIO.write(recordToSave, ofile, "fasta")
		ofile.close()

		refGenomePath = savePath
	else:
		if options.ref_path == None:
			print ("Specify the input genome, or use --randSource option")
			sys.exit(1)
		else:
			refGenomePath = options.ref_path
	# END IF
	print ("The input genome is saved in: "+refGenomePath)

	############################
	# generate reads
	############################
	print (">>>>>>>>>>>>")
	print ("Loading the reads...")
	if options.randRead:
		readsPath = ("%s/data/read_rand_n"+str(options.num_reads)\
				+"_len"+str(options.len_read)+"_"\
				+str(options.mut_rate)+"_"\
				+str(options.ins_rate)+"_"\
				+str(options.del_rate)+"_"\
				+"experi"+options.exp_index+".fa") %os.environ['FMH']
		cmd = "java -cp ./mhap-1.6/mhap-1.6.jar edu.umd.marbl.mhap.main.KmerStatSimulator %i %i %f %f %f %s > %s" %(options.num_reads, options.len_read, options.mut_rate, options.ins_rate, options.del_rate, refGenomePath, readsPath)
		subprocess.call(cmd, shell=True)
	else:
		print ("Warining: if you use the reads from the file, you need to specify the length of the reads!")
		if options.read_path == None:
			print ("Specify the input reads, or use --randRead option")
			sys.exit(1)
		else:
			readsPath = options.read_path
	# END IF
	print ("The input read is saved in: "+readsPath)

	################################
	# compute the true set
	################################
	print (">>>>>>>>>>>>")
	print ("Use Blasr to compute the true set...")
	current_path = os.environ['FMH']
	path = "%s/blasr_install/blasr/libcpp/alignment:%s/blasr_install/blasr/libcpp/hdf:%s/blasr_install/blasr/libcpp/pbdata:%s/blasr_install/hdf5/hdf5-1.8.16-linux-centos6-x86_64-gcc447-shared/lib/"%(current_path, current_path, current_path, current_path)
	os.environ["LD_LIBRARY_PATH"] += os.pathsep + path

	if options.true_set == None:
		true_set = '%s/data/ref_%s.m4'%(os.environ['FMH'],options.exp_index)
	else:
		true_set = options.true_set

	cmd = "./blasr_install/blasr/blasr %s %s --minReadLength %i --nproc %i --bestn %i -m 4 --out %s" %(readsPath, refGenomePath, options.minReadLen, options.nproc, options.bestn, true_set)
	subprocess.call(cmd, shell=True)
	print ("The reference result is saved in: "+true_set)
		

	#################################
	# set parameters for overlapping
	#################################

	print (">>>>>>>>>>>>")
	print ("Start to find overlapped reads...")
	if options.h5_path == None or options.kmer_path == None:
		h5_path = '%s/data/kmer_%s.h5'%(os.environ['FMH'], options.exp_index)
		kmer_path = '%s/data/kmer_%s.fa.txt'%(os.environ['FMH'], options.exp_index)
	else:
		h5_path = options.h5_path
		kmer_path = options.kmer_path

	if options.testing:
		print ("Using frequent kmer:")
		od = overlapDectector()
		freqKmer = True
		od.detectOverlaps(readsPath, options.kmer_size, h5_path, kmer_path, \
			options.specificity, options.min_times, options.max_times,\
			options.min_jac_threshold, options.max_jac_threshold, options.thres_gap,\
			options.jac_threshold, options.hash_threshold,\
			true_set, options.min_ovlp,\
			options.save_path_freq, freqKmer, options.kmerFreq, options.smart, options.testing)

		print ("Using all kmer:")
		od_control = overlapDectector()
		freqKmer = False
		od_control.detectOverlaps(readsPath, options.kmer_size, h5_path, kmer_path, \
			options.specificity, options.min_times, options.max_times,\
			options.min_jac_threshold, options.max_jac_threshold, options.thres_gap,\
			options.jac_threshold, options.hash_threshold,\
			true_set, options.min_ovlp,\
			options.save_path_all, freqKmer, options.kmerFreq, options.smart, options.testing)

	elif options.freqKmer:
		print ("Using frequent kmer:")
		od = overlapDectector()
		od.detectOverlaps(readsPath, options.kmer_size, h5_path, kmer_path, \
			options.specificity, options.min_times, options.max_times,\
			options.min_jac_threshold, options.max_jac_threshold, options.thres_gap,\
			options.jac_threshold, options.hash_threshold,\
			true_set, options.min_ovlp,\
			options.save_path_freq, options.freqKmer, options.kmerFreq, options.smart, options.testing)
	else:
		print ("Using all kmer:")
		od = overlapDectector()
		od.detectOverlaps(readsPath, options.kmer_size, h5_path, kmer_path, \
			options.specificity, options.min_times, options.max_times,\
			options.min_jac_threshold, options.max_jac_threshold, options.thres_gap,\
			options.jac_threshold, options.hash_threshold,\
			true_set, options.min_ovlp,\
			options.save_path_all, options.freqKmer, options.kmerFreq, options.smart, options.testing)
	print ("====END=======================================")
#######
# main
#######
if __name__ == "__main__":
	main()
##END###################################################
