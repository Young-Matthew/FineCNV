import sys
"""
combine pair end reads with two method, one is based on the location of BWA mapping result, 
the other is based on overlap sequence
"""

def locationCombined(r1, r2=None):
    if r2 is None:
        return r1.reference_start - r1.query_alignment_start, r1.query_sequence

    r1_sfc_leng = r1.cigar[-1][1] if r1.cigar[-1] == 4 else 0
    r2_sfc_leng = r2.cigar[-1][1] if r2.cigar[-1] == 4 else 0

    r1_start_pad, r1_end_pad = r1.reference_start - r1.query_alignment_start, r1.reference_end + r1_sfc_leng
    r2_start_pad, r2_end_pad = r2.reference_start - r2.query_alignment_start, r2.reference_end + r2_sfc_leng

    r1_start_pad, r1_end_pad = r1.reference_start, r1.reference_end 
    r2_start_pad, r2_end_pad = r2.reference_start, r2.reference_end 

    if r1_end_pad >= r2_start_pad and  r2_end_pad >= r1_start_pad:
        offset = r2_start_pad - r1_start_pad

        if offset >=0:
            overlap_len = min(r1.query_alignment_length - offset, r2.query_alignment_length)
            # query_alignment_length: remove softclipped base both on start and end

            r1_query_start = r1.query_alignment_start + offset
            r1_query_end   = r1_query_start + overlap_len

            r2_query_start = r2.query_alignment_start
            r2_query_end   = r2_query_start + overlap_len

            combined_sequence, diff = baseCorrection(r1.query_sequence[r1_query_start:r1_query_end],
                r2.query_sequence[r2_query_start: r2_query_end], r1.query_qualities[r1_query_start:r1_query_end],
                r2.query_qualities[r2_query_start: r2_query_end], r1)
            
            if combined_sequence is None:
                sys.stderr.write('too many mismatches-1: %d\n' % diff)
                return None, offset, overlap_len, diff, 'F'
            return combined_sequence, offset, overlap_len, diff, 'F'
        
        else:
            offset      = abs(offset)
            overlap_len = min(r2.query_alignment_length - offset, r1.query_alignment_length)
            # query_alignment_length: remove softclipped base both on start and end

            r2_query_start = r2.query_alignment_start + offset
            r2_query_end   = r2_query_start + overlap_len

            r1_query_start = r1.query_alignment_start
            r1_query_end   = r1_query_start + overlap_len

            combined_sequence, diff = baseCorrection(r2.query_sequence[r2_query_start:r2_query_end],
                r1.query_sequence[r1_query_start:r1_query_end], r2.query_qualities[r2_query_start:r2_query_end],
                r2.query_qualities[r1_query_start:r1_query_end], r1)
            
            if combined_sequence is None:
                sys.stderr.write('too many mismatches-2: %d\n' % diff)
                return None, offset, overlap_len, diff, 'R'
            return combined_sequence, offset, overlap_len, diff, 'R'


    elif r2_end_pad >= r1_start_pad and r1_start_pad >= r2_start_pad:
        offset = r1_start_pad - r2_start_pad

        if offset >=0:
            offset      = abs(offset)
            overlap_len = min(r2.query_alignment_length - offset, r1.query_alignment_length)
            # query_alignment_length: remove softclipped base both on start and end

            r2_query_start = r2.query_alignment_start + offset
            r2_query_end   = r2_query_start + overlap_len

            r1_query_start = r1.query_alignment_start
            r1_query_end   = r1_query_start + overlap_len

            combined_sequence, diff = baseCorrection(r2.query_sequence[r2_query_start:r2_query_end],
                r1.query_sequence[r1_query_start:r1_query_end], r2.query_qualities[r2_query_start:r2_query_end],
                r2.query_qualities[r1_query_start:r1_query_end], r1)

            if combined_sequence is None:
                sys.stderr.write('too many mismatches-3: %d\n' % diff)
                return None, offset, overlap_len, diff, 'R'
            return combined_sequence, offset, overlap_len, diff, 'R'

        else:        
            offset = abs(offset)
            overlap_len = min(r1.query_alignment_length - offset, r2.query_alignment_length)
            # query_alignment_length: remove softclipped base both on start and end

            r1_query_start = r1.query_alignment_start + offset
            r1_query_end   = r1_query_start + overlap_len

            r2_query_start = r2.query_alignment_start
            r2_query_end   = r2_query_start + overlap_len

            combined_sequence, diff = baseCorrection(r1.query_sequence[r1_query_start:r1_query_end],
                r2.query_sequence[r2_query_start: r2_query_end], r1.query_qualities[r1_query_start:r1_query_end],
                r2.query_qualities[r2_query_start: r2_query_end], r1)
        
            if combined_sequence is None:
                sys.stderr.write('too many mismatches-4: %d\n' % diff)
                return None, offset, overlap_len, diff, 'F'
            return combined_sequence, offset, overlap_len, diff, 'F'

    else:
        sys.stderr.write('No overlap exist %s' % r1.query_name)
        return None, None, None, None

    if r1.reference_end + r1_sfc_leng  > r2.reference_end + r2_sfc_leng:
        sys.stderr.write('%s: R1 is %s\n' %(r1.query_name, 'forward' if r1.is_reverse else 'reverse'))
        sys.stderr.write('End position for forward read is large than it in reverse strand: %d > %d\n' %
            (r1.reference_end + r1_sfc_leng, r2.reference_end + r2_sfc_leng))
        sys.stderr.write('R1: %d + %d, R2: %d + %d\n\n'% (r1.reference_end, r1_sfc_leng, r2.reference_end, r2_sfc_leng))

    
#"Copied from https://github.com/OpenGene/fastp/blob/master/src/overlapanalysis.cpp"
def seqCombined(r1, r2):
    seq1, seq2 = r1.query_sequence, r2.query_sequence
    qual1, qual2 = r1.query_qualities, r2.query_qualities

    len1, len2 = len(seq1), len(seq2)
    complete_compare_require = 50
    overlap_len, offset, diff, overlapDiffLimit, overlapRequire = 0, 0, 0, 5, 30    

    # forward
    # a match of less than overlapRequire is considered as unconfident
    while offset < len1 - overlapRequire:
        # the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2)
        diff = 0 
        for i in range(overlap_len):
            if seq1[offset + i] != seq2[i]:
                diff += 1
                if diff >= overlapDiffLimit and i < complete_compare_require:
                    break

        if diff < overlapDiffLimit or (diff > overlapDiffLimit and i > complete_compare_require):
            return seq2[:overlap_len], offset, overlap_len, diff, 'F'

            """
            combined_sequence = seq1[:offset]
            for i in range(overlap_len):
                # correct various base in overlap area, Novaseq quality: 11, 25, 37, 40
                if seq1[offset + i] != seq2[i]:
                    if  qual1[offset + i] > qual2[i]:
                        combined_sequence += seq1[i]
                    else:
                        combined_sequence += seq2[i]
                else:
                    combined_sequence += seq1[i]

            combined_sequence += seq2[overlap_len:]

            return combined_sequence, offset, overlap_len, diff
            """
        offset += 1

    """
    reverse in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    check if distance can get smaller if offset goes negative
    this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    we go reversely
    """

    offset = 0
    while offset > -(len2-overlapRequire):
        # the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1,  len2- abs(offset))
        diff = 0

        for i in range(overlap_len):
            if seq1[i] != seq2[-offset + i]:
                diff += 1
                if diff >= overlapDiffLimit and i < complete_compare_require:
                    break

        if diff < overlapDiffLimit or (diff > overlapDiffLimit and i > complete_compare_require):
            return seq1[:overlap_len], offset, overlap_len, diff, 'R'

            if diff == 0:
                return seq1[:overlap_len], offset, overlap_len, diff
            else:
                combined_sequence = seq1[:overlap_len]
                for i in range(overlap_len):
                    if seq1[i] != seq2[-offset + i]:
                        if qual1[i] < qual2[-offset + i]:
                            combined_sequence[i] = seq2[-offset + i] 

                return combined_sequence, offset, overlap_len, diff
        offset += 1

    return None, None, None, None


def baseCorrection(seq1, seq2, qual1, qual2, r1):
    if len(seq1) != len(seq2):
        raise Exception(r1.query_name , ' :various read length: ', seq1, seq2, len(seq1), len(seq2))

    nmismatch = 0
    combined_sequence = seq1

    for i, base in enumerate(seq1):
        if base != seq2[i]:
            nmismatch += 1
            if qual1[i] < qual2[i]:
                combined_sequence = combined_sequence[:i] + seq2[i]

    if len(seq1) < nmismatch * 10:
        return None, nmismatch

    return combined_sequence, nmismatch


