import sys
import os
import pyfaidx
import time
from datetime import datetime
import logging
from pysam import AlignmentFile
from nmismatch import nmismatch
from breakpoints import breakpoints

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

for read in bam:    
    try:
        if read.cigarstring is None:
            continue
    except:
        continue

    if read.reference_name != lastChr:

        for name, single_read in all_reads.items():
            if len(single_read) != 1:
                raise Exception('single read:', len(single_read))

            single_read = single_read[0]
            cigar = single_read.cigar
            if cigar[0][0] == 4:
                position = single_read.reference_end
            elif cigar[-1][0] == 4:
                position = single_read.reference_start
            else:
                continue

            if position in sfc_reads:
                sfc_reads[position][single_read.reference_name] = (single_read)
            else:
                sfc_reads[position] = {}

        all_reads = {}
        print(lastChr)
        breakpoints(sfc_reads)
        sfc_reads = {}
        lastChr = read.reference_name


    if read.query_name in all_reads:

        r1, r2 = all_reads[read.query_name][0], read
        if any((nmismatch(r1) >4, nmismatch(r2) >4)):
            del all_reads[read.query_name]
            continue

        if r2.is_reverse and  r1.is_reverse:
            del all_reads[read.query_name]
            continue

        if r2.is_reverse:
            forward, reverse = r1, r2
        else:
            forward, reverse = r2, r1

        if forward.cigar[-1][0] == 4:
            # Softclip on 3' end of forward sequence
            position = forward.reference_end
        elif reverse.cigar[0][0] == 4: 
            # Softclip on 5' end of reverse sequence
            position = reverse.reference_start 

            """
            if r1.cigar[-1][0] == 4:
                if r2.reference_end != r1.reference_start:
                    sys.stderr.write('PE: %s both has softclipped bases, location contradiction, r1:%d != r2:%d,' %
                        (read.query_name, r2.reference_end, r1.reference_start))
                    exit()
            """
        else:
            del all_reads[read.query_name]
            continue

        if position in sfc_reads:
            sfc_reads[position][read.reference_name] = (forward, reverse)
        else:
            sfc_reads[position] = {}

    else:
        all_reads[read.query_name] = [read]


bam.close()
t2 = datetime.now()
t_used = (t2 - t1).seconds
sys.stderr.write('%d seconds used' % (t_used))