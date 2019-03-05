

def nmismatch(read):
    try:
        md = read.get_tag('MD')

        deletion, snp = False, 0
        for x in md:
            if x == '^':
                deletion = True
            elif x in 'ATCG':
                snp += not deletion
            else:
                deletion = False

        return snp + sum(1 for op, length in read.cigartuples if op in (1, 2))

    except KeyError:
        pass

    try:
        nm = read.get_tag('NM')
        indel = sum(length-1 for op, length in read.cigartuples if op in (1, 2))

        return nm - indel

    except KeyError:
        pass

    except TypeError:
        pass

    return -1

