import subprocess

'''
The details of each options is explained in doall.py
You can check it by using:
python3 doall.py --help
'''
'''
# example 1 (same as last Version):
# compare the method using frequent kmers and the method using all kmers
# all other settings are the same
cmd = "python3 doall.py --randSource --randRead"
subprocess.call(cmd, shell=True)
# Note: here the "--reuse" option enables you to use the exact same genome
# 	and reads used in the previous run
cmd = "python3 doall.py --freqKmer --reuse"
subprocess.call(cmd, shell=True)
'''

# example 2 (updated):
# run multiple times with different min-hash threshold (1st filter threshold)
# and the jaccard threshold (2nd filter threshold)
# Note1: you may want to save the final precision and recall. In order to do this,
# you need to use the option "--save" followed by the file you use to save the result
# the save format is:
# ###file start:
# hThreshold1, jThreshold1, precision1, recall1
# hThreshold2, jThreshold2, precision2, recall2
# ........

# Note2:
# here the jThreshold (the 2nd threshold) can be set using --minjThres and --maxjThres.
# These two options provide the lower and upper bounds of the jThreshold.
# The option jGap gives the increase gap. Namely, the program will try jThreshold:
# minjThres, minjThres + jGap, minjThres + 2*jGap, ...

# Note3:
# here the hThreshold is obtained according to jThreshold.
# hThreshold = jThreshold/ given_times
# And you can set the range of the times using --minTimes and --maxTimes
# e.g. --minTimes 2 --maxTimes 3 means the program will try both:
# hThreshold = jThreshold/2 and hThreshold = jThreshold/3
save_path = "data/precision_and_recall_all_kmer.txt"
cmd = "python3 doall.py --randSource --randRead --minjThres 0.01 --maxjThres 0.10 --jGap 0.02 --minTimes 2 --maxTimes 3 --save %s" % save_path
subprocess.call(cmd, shell=True)
# sample output in the save_path
# 0.005   0.01    0.8675496688741722      0.9978237214363439
# 0.015   0.03    0.9969861362266426      0.8998911860718172
# 0.025   0.05    0.9983818770226537      0.6713819368879217
# 0.034999999999999996    0.06999999999999999     1.0     0.4542981501632209
# 0.045   0.09    1.0     0.1441784548422198
# 0.0033333333333333335   0.01    0.8629107981220657      1.0 
# 0.01    0.03    0.9958992384299942      0.9249183895538629
# 0.016666666666666666    0.05    0.9984313725490196      0.6926006528835691
# 0.02333333333333333     0.06999999999999999     1.0     0.48258977149075083
# 0.03    0.09    1.0     0.2578890097932535

save_path = "data/precision_and_recall_freq_kmer.txt"
cmd = "python3 doall.py --reuse --freqKmer --minjThres 0.01 --maxjThres 0.10 --jGap 0.02 --minTimes 2 --maxTimes 3 --save %s" % save_path
subprocess.call(cmd, shell=True)


# example 3 (same as previous Version):
# you may want to use diffrent settings of genome and reads as input data
# let's take the "frequent kmer" mode as an example
cmd = "python3 doall.py --randSource --lenRef 100000 --randRead --nRead 500 --lenRead 1000 --minOvlp 200 --freqKmer"
subprocess.call(cmd, shell=True)

# Also, you can set different mutation rate, indel rate when generating reads
# Note: when you want to use diffrent genome or reads or both from the previous command, do NOT use "--reuse" option.
# Because, if you do, the setting won't change, and you will use the same data as previous one.
cmd = "python3 doall.py --randSource  --randRead --mut 0.1 --ins 0.02 --del 0.02 --freqKmer"
subprocess.call(cmd, shell=True)

# of course you can use your own genome, reads intead of random generated ones
# Note: to do so, do NOT use the "--randSource" and "--randRead" option
cmd = "python3 doall.py --ref data/refGenome_rand_len1000000.fa  --read data/read_rand_n500_len10000_0.05_0.01_0.01.fa --freqKmer"
subprocess.call(cmd, shell=True)
# and you also can only use your own genome and randomly generate reads from it
cmd = "python3 doall.py --ref data/refGenome_rand_len1000000.fa  --randRead --freqKmer"
subprocess.call(cmd, shell=True)


# example 4 (same as previous Version):
# As we discussed, the "--smartFreq" option enables you to smartly choose the kmer frequence range.
# the program will sort the kmer list and only use the first (genome)-length kmers for use
# the frequence is at least 2
cmd = "python3 doall.py --randSource --randRead --freqKmer --smartFreq"
subprocess.call(cmd, shell=True)

###END##########################################
