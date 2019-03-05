import os, sys, re, pysam

usage = '''cut softclipped bases from the read and print according to it's position
eg: python3 sfc_short_bases.py bam chr1 10000  &
'''

if len(sys.argv) != 4:
    print(usage)
    exit()

outfile = ':'.join(sys.argv[2:4]) + '_' + os.path.splitext(os.path.basename(sys.argv[1]))[0] + '.txt'

with open(outfile, 'w') as out:
    start, end = int(sys.argv[3]) - 150 - 1, int(sys.argv[3]) + 150 - 1
    bam =  pysam.AlignmentFile(sys.argv[1], 'rb')
    extra_space, start_sfc_len = 100, 0
    sfc_read_number = 0

    for read in bam.fetch(sys.argv[2], start, end):
        if 'S' not in read.cigarstring:
            continue
        sfc_read_number += 1
        space_length = extra_space + (read.reference_start - start)
        if space_length < 0:
            space_length = 0

        if read.cigar[0][0] == 4:
            if read.cigar[0][1] > extra_space and read.cigar[0][1] > start_sfc_len:
                start_sfc_len =  read.cigar[0][1]
            space_length = extra_space + (read.reference_start - read.cigar[0][1] - start )

        for i in range(space_length):
            print(' ', end='', file = out)

        query_id = -1
        strand = '-' if read.is_reverse else '+'
        sfc_seq = ''

        for cigar in read.cigar:
            if cigar[0] == 0:
                query_id += cigar[1]
                for i in range(cigar[1]):
                    print(strand, end = '', file = out)

            elif cigar[0] == 1:
                query_id += cigar[1]

            elif cigar[0] == 2:
                for i in range(cigar[1]):

                    print(strand, end = '', file = out)

            elif cigar[0] == 4:
                sfc_seq = read.query_sequence[query_id + 1:query_id + 1 + cigar[1]]
                for i in range(cigar[1]):
                    query_id += 1
                    print(read.query_sequence[query_id], end ='', file = out)

            else:
                print('Unkown cigar', read.cigarstring, read.query_name)
                exit()
        print(file=out)
        #print(read.reference_start, read.cigar, read.query_name,'\n', read.query_sequence)

        print('>%s\n%s' %(read.query_name + '_' + str(sfc_read_number), sfc_seq))