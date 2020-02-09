'''
This script fills replaces soft-clipped bases with mapped bases in the CIGAR string of long reads mapped to the GLOBE vector so that readthrough transcript coverage can be computed.
Input must be a SAM file.
'''
import os
import sys
import re

file_name = "/path/to/input.sam"
f = open(file_name)

out_good = open("/path/to/output_new_cigar.sam", "w")

for line in f:
    line = line.strip('\n')
    col = line.split('\t')
    if col[0][0] == '@': # write header lines into output file
        out_good.write(line + "\n")
        continue
    cigar = str(col[5]) # gets cigar string from SAM file
    #last100 = col[9][:100] # gets sequence of 100 leftmost bases in read
    map_start = col[3] # position of mapping start (end of soft-clipping)
    if "S" in cigar: # if soft-clipping present in read and...
        if re.findall('[A-Z]', cigar)[0] == 'S' and re.findall('[A-Z]', cigar)[1] == 'M': # ...if first cigar operator is 'S' and next operator is 'M'
            position_first_S = cigar.find('S') # finds position of left-most 'S' in cigar string
            length_first_S = int(cigar.split("S")[0]) # get digits before first 'S' in CIGAR string
            position_first_M = cigar.find('M') # find first M
            length_first_M = int(cigar.split("M")[0][position_first_S+1:]) # get length of first M
            new_M_length = length_first_S + length_first_M
            remainder = cigar[position_first_M+1:]
            new_cigar = str(new_M_length) + "M" + str(remainder)
            new_POS = int(map_start) - int(length_first_S)

            out_good.write(col[0] + "\t" + col[1] + "\t" + col[2] + "\t" + str(new_POS) + "\t" + col[4] + "\t" + str(new_cigar) + "\t" + col[6] + "\t" + col[7] + "\t" + col[8] + "\t" + col[9] + "\t" + col[10] + "\n")
        else:
            out_good.write(line + "\n")
    else:
        out_good.write(line + "\n")
