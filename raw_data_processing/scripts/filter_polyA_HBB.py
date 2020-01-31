'''
This script filters out reads with long soft-clipped alignments and reads with non-templated polyT sequences from HBB targeted long read sequencing
Note this script is specifically only for reads mapped to GLOBE vector!!!
Author: Kirsten Reimer
'''
import os
import sys
import re

# file = sys.argv[1]
# with open(file, 'r') as f:
file_name = "/path/to/input/file.sam"
f = open(file_name)

out_readthrough_bad = open("/path/to/output/readthrough_bad.sam", "w")
out_chimera_bad = open("/path/to/output/chimera_bad.sam", "w")
out_polyA_good = open("/path/to/output/polyA_good.sam", "w")
out_good = open("/path/to/output/polyA_filtered.sam", "w")

for line in f:
    line = line.strip('\n')
    col = line.split('\t')
    if col[0][0] == '@': # write header lines into output file
        out_good.write(line + "\n")
        out_polyA_good.write(line + "\n")
        out_chimera_bad.write(line + "\n")
        out_readthrough_bad.write(line + "\n")
        continue
    cigar = str(col[5]) # gets cigar string from SAM file
    last100 = col[9][:100] # gets sequence of 100 leftmost bases in read
    if "S" in cigar: # if soft-clipping present in read and...
        if re.findall('[A-Z]', cigar)[0] == 'S': # ...if read contains soft-clipping as the first cigar operator
            index = cigar.find('S') # finds position of left-most S in cigar string
            length_clipped = int(cigar.split("S")[0]) # get digits before first S in CIGAR string
            clipped_bases = col[9][:length_clipped] # get sequence of clipped bases from start to first S
            read_end = int(col[3]) - length_clipped # read end is start of mapped sequence plus clipped bases (subtract because all reads are on - strand)
            clip_start = int(col[3]) - 1 # start of mapped sequence minus 1
            if read_end < 2448 and last100.count("T") >= 70: # if soft clipping begins after the HBB insertion (mapping to downstream genomic region)
                out_readthrough_bad.write(line + "\n") #filter out readthrough reads with polyT tails
            elif clip_start in range(4426,4480) and clipped_bases.count('T')/len(clipped_bases) >=0.7: # if clipping starts within 50nt of canonical polyA site and clipped bases are greater than 70% T
                out_polyA_good.write(line + "\n") # filter out normal polyT reads
            elif clip_start in range(4426,4480) and clipped_bases.count('T') >=4: # if clipping starts within 50 nt of canonical polyA site and has more than 4T (accounts for even very small tails)
                out_polyA_good.write(line + "\n") # filter our normal polyT reads
            elif read_end > 2448 and length_clipped >= 20: # any other long mismatch within GLOBE sequnce
                out_chimera_bad.write(line + "\n")
            elif read_end < 2448 and clip_start > 2448: # any other long mismatch within GLOBE sequence
                out_chimera_bad.write(line + "\n")
            else:
                out_good.write(line + "\n")
        else:
            out_good.write(line + "\n")
    else:
        out_good.write(line + "\n")
