```
calculate_distance_to_splicing.py
This script calculates the block size of the last exon as well as coordinates of the upstream intron from a BED12 file
```

filepath = "../input.bed"
file = open(filepath)
output = open("../output_blocksize_introncoord.txt", "w")

for line in file:
    line_split = line.split('\t')
    exon_size_list = line_split[10] # gets list of block (exon) sizes from BED12 file
    intron_size_list = line_split[11] # gets list of block starts
    strand = line_split[5]
    chr = line_split[0]
    exon_size = exon_size_list.split(",")
    intron_size = intron_size_list.split(",")
    last_intron_size_plus = str(int(intron_size[-1]) - int(intron_size[-2]) - int(exon_size[-2]))
    last_intron_size_minus = str(int(intron_size[1]) - int(exon_size[0]))

# # print last exon size, as well as last intron size
#     if strand == '+':
#         output.write(exon_size[-1] + "\t" + last_intron_size_plus + "\n")
#     if strand == '-':
#         output.write(exon_size[0] + "\t" + last_intron_size_minus + "\n")

# print last exon size, as well as coordinates of immediately upstream intron
    if strand == '+':
        read_end = line_split[2]
        upstream_intron_start_plus = (int(read_end) - int(exon_size[-1]) - int(last_intron_size_plus))
        upstream_intron_end_plus = (int(read_end) - int(exon_size[-1]))
        output.write(exon_size[-1] + "\t" + chr + "\t" + str(upstream_intron_start_plus) + "\t" + str(upstream_intron_end_plus) + "\n")
    if strand == '-':
        read_end = line_split[1]
        upstream_intron_end_minus = (int(read_end) + int(exon_size[0]))
        upstream_intron_start_minus = (int(read_end) + int(exon_size[0]) + int(last_intron_size_minus))
        output.write(exon_size[0] + "\t" + chr + "\t" + str(upstream_intron_end_minus) + "\t" + str(upstream_intron_start_minus) + "\n")
