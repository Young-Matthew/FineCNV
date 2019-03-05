"""
output reads with softclipped in bamFile
"""
import sys, pysam, os
if len(sys.argv) != 2:
    print('bam must be provided')
    exit()

def try_append(d, k, v):
    if k not in d:
        d[k] = [v]
    else:
        d[k].append(v)

bam = sys.argv[1]
outbam = os.path.splitext(os.path.basename(bam))[0] + '.sfc.bam'
inBam = pysam.AlignmentFile(bam, 'rb')
tmpBam = pysam.AlignmentFile('temp.bam', 'wb', template = inBam)
read_out = {}
# read_id -> read stored address in memory

last_chr = 'chr1'
for read in inBam:

    if read.reference_name != last_chr:
        last_chr = read.reference_name
        for read_id in read_out.keys():
            if any('S' in single.cigarstring for single in read_out[read_id]):
                for single in read_out[read_id]:
                    tmpBam.write(single)
        read_out = {}
    else:
        try_append(read_out, read.query_name, read)

tmpBam.close()
pysam.sort('-o', outbam, 'temp.bam')
os.remove('temp.bam')