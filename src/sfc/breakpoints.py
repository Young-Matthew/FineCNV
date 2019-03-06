from overlapanalysis import overlap

molecule = []
def breakpoints(reads):

    for position in list(reads.keys()):
        if len(reads[position].keys()) < 10:
            print('depth < 10: %d' % position)
            del reads[position]
            continue

        print(position)
        for name, read in reads[position].items():
            if len(read) == 1:
                continue

            r1, r2 = read
            """
            if r1.is_reverse and r1.reference_start < r2.reference_start:
                print('1:', r1.reference_name, r1.is_reverse, r1.reference_start, r2.reference_start)
                
            elif r2.is_reverse and r1.reference_start > r2.reference_start:
                print('2:', r1.reference_name, r1.is_reverse, r1.reference_start, r2.reference_start)
            """
            location, sequence = overlap(r1, r2)
            print(location, sequence)

        print()



