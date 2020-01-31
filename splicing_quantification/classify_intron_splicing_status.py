'''
classify_intron_splicing_status.py
Author: Kirsten Reimer
This script filters only spliced reads from all intron-spanning reads.
Input file should be a BED12 file that has been intersected with a BED6 file of intron coordinates +/- 1 nt using bedtools intersect, e.g. "bedtools intersect -b introns_plus_1.bed -a data.bed -split -F 1 > junction_spanning_reads.bed"
'''

bed_filepath = "../input_junction_spanning_reads.bed"
bed_file = open(bed_filepath)
output = open("../junction_spanning_reads_splicing_status.bed", "w")

for line in bed_file:
    col = line.rstrip().split("\t")
    int_coord1 = col[1] 
    int_coord2 = col[2]
    intron = ((int(int_coord1)+1), (int(int_coord2)-1)) # create a pair of the intron coordinates

    readStart = col[6]
    readEnd = col[7]
    blockStarts = col[11].split(",")
    blockSizes = col[10].split(",")

    n = int(col[9]) # get number of blocks to use as counter
    junclist = [] # create an empty list for coordinates of junctions in each read
    
    # to make a search window around the intron annotation to check if a ready is really unspliced or spliced to a neighboring intron 
    false_count = 0 #false_count to be used for if statments below 
    intron_range0A = range(intron[0]-10,intron[0])
    intron_range0B = range(intron[0],intron[0]+10)
    intron_range1A = range(intron[1]-10,intron[1])
    intron_range1B = range(intron[1],intron[1]+10) 

    for i in range(1,n):
        juncStart = int(readStart) + int(blockSizes[i-1]) + int(blockStarts[i-1])
        juncEnd = int(readStart) + int(blockStarts[i])
        junc = (juncStart,juncEnd)
        junclist.append(junc) # add each junction to a list of junctions in the read 

            # if the junction annotation is within a window of the search window, false_count = 1. The long-read will not be counted as unspliced
        if junc[0] in intron_range0A:
            false_count = 1 
        elif junc[0] in intron_range0B: 
            false_count = 1
        elif junc[0] in intron_range1A: 
            false_count = 1
        elif junc[0] in intron_range1B: 
            false_count = 1
        elif junc[1] in intron_range0A:
            false_count = 1 
        elif junc[1] in intron_range0B: 
            false_count = 1
        elif junc[1] in intron_range1A: 
            false_count = 1
        elif junc[1] in intron_range1B: 
            false_count = 1
            

    if intron in junclist: # look for both exact intron coordinates in list of junctions in the read
    	output.write(line.rstrip() + "\t" + "spliced_exact_SJ" + "\n")

    elif intron[0] in (i[0] for i in junclist): # look for first intron coordinate in list of first junctions in the read
    	output.write(line.rstrip() + "\t" + "spliced_half_SJ" + "\n")

    elif intron[1] in (i[0] for i in junclist): # look for second intron coordinate in list of second junctions in the read
    	output.write(line.rstrip() + "\t" + "spliced_half_SJ" + "\n")

    elif false_count == 1: # if the junction is found within a window, dis-regard the long read 
            output.write(line.rstrip() + "\t" + "false_count" + "\n")

    else:
    	output.write(line.rstrip() + "\t" + "unspliced" + "\n")
