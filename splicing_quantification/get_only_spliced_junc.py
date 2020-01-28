'''
get_only_spliced_junc.py
Author: Kirsten Reimer
Usage: python3 get_only_spliced_junc.py input_junction_spanning_reads.bed > output_spliced_junction_spanning_reads.bed
This script filters only spliced reads from all intron-spanning reads.
Input file should be a BED12 file that has been intersected with a BED6 file of intron coordinates +/- 1 nt using bedtools intersect, e.g. "bedtools intersect -b introns_plus_1.bed -a data.bed -split -F 1 > junction_spanning_reads.bed"
'''

import os
import sys
import re

file = sys.argv[1]
with open(file, 'r') as bed_file:

  for line in bed_file:
      col = line.rstrip().split("\t")
      int_coord1 = col[1]
      int_coord2 = col[2]
      intron = str(int(int_coord1)+1) + str(int(int_coord2)-1) # create a string of the intron coordinates

      readStart = col[6]
      readEnd = col[7]
      blockStarts = col[11].split(",")
      blockSizes = col[10].split(",")

      n = int(col[9]) # get number of blocks to use as counter
      junclist = [] # create an empty list for coordinates of junctions in each read

      for i in range(1,n):
          juncStart = int(readStart) + int(blockSizes[i-1]) + int(blockStarts[i-1])
          juncEnd = int(readStart) + int(blockStarts[i])
          junc = str(juncStart) + str(juncEnd)
          junclist.append(junc) # add each junction to a list of junctions in the read

      if intron in junclist: # look for exact intron coordinates in list of junctions in the read
          print(line.rstrip() + "\n")
