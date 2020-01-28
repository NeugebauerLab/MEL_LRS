```
classify_cigar.py
Author: Kirsten Reimer
This script sorts a SAM file into spliced and unspliced reads based on whether there is an "N" operator in the CIGAR string.
```

filepath = "../input.sam"
file = open(filepath)
spliced = open("../output_spliced.sam", "w")
unspliced = open("../output_unspliced.sam", "w")


for line in file:
    line_split = line.split('\t')
    if line_split[0][0] == '@':
        spliced.write(line)
        unspliced.write(line)
        continue
    cigar = line_split[5] # gets cigar string from SAM file

    if "N" in cigar:
        spliced.write(line)

    elif "N" not in cigar:
        unspliced.write(line)

spliced.close()
unspliced.close()
