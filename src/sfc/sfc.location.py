# scan raw bam file to find the location where covered with 
# multiple reads flagged with softclip, as I believe it might 
#contain the information of CNV.

import sys, pysam, re, os

bam = pysam.AlignmentFile(sys.argv[1], 'rb')
outfile = os.path.splitext(os.path.basename(sys.argv[1]))[0] + '.sfc.txt'
location = {}
# chr -> position -> freq

for read in bam:
    read_align = read.aligned_pairs

    if read.cigar[0][0] == 4:
        # Softclip on 5' end
        if read.reference_name not in location:
            location[read.reference_name] = {}

        if read.reference_start + 1 not in location[read.reference_name]:
            location[read.reference_name][read.reference_start +1] = 1
        else:
            location[read.reference_name][read.reference_start +1] += 1
        continue
    elif read.cigar[-1][0] == 4:
        # Softclip on 3' end
        if read.reference_name not in location:
            location[read.reference_name] = {}
            
        if read.reference_end + 1 not in location[read.reference_name]:
            location[read.reference_name][read.reference_end +1] = 1
        else:
            location[read.reference_name][read.reference_end +1] += 1
        continue
    if 'x' in read.reference_name or 'X' in read.reference_name:
        break

bam.close()
with open(outfile, 'w') as out:
    print('chr\tposition\tfreq', file = out)
    for chr_number in range(1,23):
        chr = 'chr' + str(chr_number)
        pos_set = sorted([i for i in location[chr].keys()])
        for pos in pos_set:
            print('%s\t%d\t%d' % (chr, pos, location[chr][pos]), file = out)
