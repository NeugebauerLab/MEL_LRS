'''
This script filters out unmapped reads and reads with soft-clipped alignments of polyA or polyT sequence (non-templated As and Ts)
Author: Kirsten Reimer
Usage: python filter_polyA.py file.sam > filtered_file.sam
Sep 2019
'''

import os
import sys
import re

file = sys.argv[1]
with open(file, 'r') as f:

    for line in f:
        line = line.strip('\n')
        col = line.split('\t')
        if col[0][0] == '@': # write header lines into output file
            print(line)
            continue
        cigar = col[5] # gets cigar string from SAM file

        if cigar.count('S') >= 1: # searches for S in cigar string indicating soft-clipped
            index = cigar.find('S') # finds position of S in cigar string
            if index <= 3: # if index is small then clipped bases are on the front of the read...should be polyT
                length_clipped = int(cigar.split("S")[0]) # get position of first S in cigar string
                clipped_bases = col[9][:length_clipped] # get sequence of clipped bases from start to first S
                if  clipped_bases.count('T') >= 4 and \
                    clipped_bases.count('T')/len(clipped_bases) >=0.9: # if minimum 4 T's and T content is >90% of soft-clipped bases
                        continue # skip lines that have polyT at beginning of read
                else:
                    print(line) # print lines that do not have polyT

            if (index + 1) == len(cigar): # if the S is at the end of the cigar string
                m = re.search('[0-9]{0,3}(?=S)' , cigar) # find number before last S in cigar string
                length_clipped = int(m.group(0)) #get length of clipped bases at end of read...should be polyA
                clipped_bases = col[9][:(-1 * (length_clipped + 1)):-1] # get sequence of clipped bases before last S
                if  clipped_bases.count('A') >= 4 and \
                    clipped_bases.count('A')/len(clipped_bases) >=0.9: # if minimum 4 A's and A content is >90% of soft-clipped bases

                        continue # do not write lines that have polyA at end of read
                else:
                    print(line) # write lines that do not have polyA at the end
        else:
            print(line)
