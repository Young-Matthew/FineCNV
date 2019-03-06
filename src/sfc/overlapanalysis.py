"Copied from https://github.com/OpenGene/fastp/blob/master/src/overlapanalysis.cpp"

def overlap(r1, r2=None):
    if r2 is None:
        return r1.reference_start - r1.query_alignment_start, r1.query_sequence

    if r1.reference_start == r2.reference_start:
        if r1.reference_end != r2.reference_end:
            raise Exception('Range should be same', r1.reference_name)
        return r1.reference_start - r1.query_alignment_start, r1.query_sequence

    else:
        r1_sfc_leng = r1.cigar[-1][1] if r1.cigar[-1] == 4 else 0

        if r1.reference_end + r1_sfc_leng <= r2.reference_start - r2.query_alignment_start:
        # overlap truly exist
            start, end = r2.reference_start - r2.query_alignment_start - (r1.reference_start - r1.query_alignment_start),
            r1.reference_end + r1_sfc_leng - (r1.reference_start - r1.query_alignment_start) + 1

            if mismatch(r1.query_sequence[start:], r2.query_sequence[:end]):
                return r1.reference_start - r1.query_alignment_start, r1.query_sequence + r2.query_sequence[end:]

            else:
                raise Exception('overlap issue:', r1.query_name)

        else:
            # Non overlap 
            empty = 'N' * (r2.reference_start - r2.query_alignment_start - (r1.reference_end + r1_sfc_leng ))
            return r1.reference_start - r1.query_alignment_start, r1.query_sequence + empty + r2.query_sequence

def mismatch(seq1, seq2):
    if len(seq1) != len(seq2):
        return False
    nmismatch = 0
    for i, base in enumerate(seq1):
        if base != seq2[i]:
            nmismatch += 1
        if nmismatch > 4:
            return False
    return True


