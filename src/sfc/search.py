
molecule = []
def search(reads):
    for position in list(reads.keys()):
        if len(reads[position].keys()) < 10:
            del reads[position]
            continue

        for name, read in reads[position].items():
            if len(read) == 1:
                continue


