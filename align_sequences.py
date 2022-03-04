#!/usr/bin/env python

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError
import numpy as np

# a simple function to read the name and sequence from a file
# The file is expected to have just one contig/sequence. This function
# checks the assumption and complains if it is not the case.
def read_single_contig_fasta(filename):
    names = []
    sequences = []
    with open(filename, 'r') as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T"]:
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]

def smith_waterman(seq1, seq2, match, mismatch, gapopen, gapextend):
    max_score = 0
    alnseq1 = ""
    alnseq2 = ""
    
    M = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype=int).astype(float)
    M[range(1,len(seq1)+1),0] = np.ones(len(seq1), dtype=int).astype(float)*-1
    M[0,range(1,len(seq2)+1)] = np.ones(len(seq2), dtype=int).astype(float)*-1

    for i in range(1,len(seq1)+1):
        k = 0 #reset size of gap for each new row
        for j in range(1,len(seq2)+1):
            if seq1[i-1] == seq2[j-1]:
                score = match #match score
            else:
                score = mismatch #mismatch score
            gap = gapopen + k*gapextend #gap score
            M[i,j] = max(0, M[i-1,j]-gap, M[i-1,j-1]+score, M[i,j-1]-gap)
            if M[i,j] == (M[i-1,j]-gap) or M[i,j] == (M[i,j-1]-gap):
                k = k + 1 #we chose a gap, so size of gap increments by 1
            else:
                k = 0 #reset gap size if you have a match 
                
    #initialize aligned sequences to fill in the traceback algorithm
    alnseq1 = ["*" for r in range(np.max([len(seq1),len(seq2)]))] 
    alnseq2 = ["*" for r in range(np.max([len(seq1),len(seq2)]))]

    max_score = np.max(M) #save the max score from M as output
    
    idx_i = np.argmax(np.argmax(M,axis=1)) #starting points for the traceback algorithm
    idx_j = np.argmax(np.argmax(M,axis=0))
    i = idx_i
    j = idx_j
    #a loop to trace the optimal path back through M to return aligned sequences
    while i>0 or j>0:
      chosen_direction = max(M[i-1,j-1],M[i,j-1],M[i-1,j]) #choose direction that gives you largest score
      if chosen_direction == M[i-1,j-1]:
        if len(seq1)>len(seq2): #length of the aligned sequences = length of longest initial sequence, so we save assign characters accordingly
          alnseq1[i-1] = seq1[i-1]
          alnseq2[i-1] = seq2[j-1]
        else:
          alnseq1[j-1] = seq1[i-1]
          alnseq2[j-1] = seq2[j-1]
        i = i-1
        j = j-1
      if chosen_direction == M[i,j-1]:
        if len(seq1)>len(seq2):
          alnseq1[i-1] = ["-"] #insert a gap
          alnseq2[i-1] = seq2[j-1]
        else:
          alnseq1[j-1] = ["-"] #insert a gap
          alnseq2[j-1] = seq2[j-1]      
        j = j-1
      if chosen_direction == M[i-1,j]:
        if len(seq1)>len(seq2):
          alnseq1[i-1] = seq1[i-1]
          alnseq2[i-1] = ["-"] #insert a gap
        else:
          alnseq1[j-1] = seq1[i-1]
          alnseq2[j-1] = ["-"] #insert a gap
        i = i-1

    return max_score, alnseq1, alnseq2
    
def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and 
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(seq1, seq2, 
                                  match, mismatch, gapopen, gapextend)
    
    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:],
                     "hm:x:g:e:",
                     ["help", "match=", "mismatch=", "gapopen=", "gapextend="])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)
