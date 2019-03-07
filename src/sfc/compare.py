import sys
import os
import pyfaidx
import time
from datetime import datetime
import logging
from pysam import AlignmentFile
from nmismatch import nmismatch
from overlapanalysis import locationCombined, seqCombined

if len(sys.argv) < 2:
    print('python3 ' + sys.argv[0] + ' BAM ' + ' OUTDir ')
    exit()
t1 = datetime.now()

bam = AlignmentFile(sys.argv[1], 'rb')
lastChr = 'chr1'
position = 0
all_reads = {}
sfc_reads = {}
# position -> read_name -> read_address in memory

i = 0
for read in bam:    
    try:
        if read.cigarstring is None:
            continue
    except:
        continue    

    i += 1
    if i % 2 == 0:
        print(i // 2, 'molecule')

    if read.query_name in all_reads:
        res1 = list(locationCombined(all_reads[read.query_name], read))
        res2 = list(seqCombined(all_reads[read.query_name], read))
        r1 = all_reads[read.query_name]
        if res1 != res2:
            sys.stderr.write('Various compare result: %s\n%r\n%r\n' % (read.query_name, res1, res2))
            sys.stderr.write('%d: %s\n%d: %s\n\n' % (r1.reference_start, r1.query_sequence, read.reference_start, read.query_sequence))
        else:
            if res1[0] is None:
                print('Empty ', end = '\n')
            print('identical result %s' % read.query_name)

        del all_reads[read.query_name]
    
    else:
        all_reads[read.query_name] = read

t2 = datetime.now()
print((t2-t1).seconds, 'used')
