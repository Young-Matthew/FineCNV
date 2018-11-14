from skgenome import tabio
from __future__ import print_function, absolute_import, division
from builtins import range
import logging
import pandas as pd
from merge import merge

# hopefully, most of funcs are copied from cnvkit.....................Good luck!

def do_target(bait_arr, annotate=None, do_short_names=False, do_split=False,
              avg_size=200/.75):
    """Transform bait intervals into targets more suitable for CNVkit."""
    tgt_arr = bait_arr.copy()
    # Drop zero-width regions
    tgt_arr = tgt_arr[tgt_arr.start != tgt_arr.end]
    tgt_arr = tgt_arr.subdivide(avg_size, 0)
    return tgt_arr

def subdivide(table, avg_size, min_size=0, verbose=False):
    return pd.DataFrame.from_records(
        _split_targets(table, avg_size, min_size, verbose),
        columns=table.columns)

def _split_targets(regions, avg_size, min_size, verbose):
    """Split large regions into smaller, consecutive regions.

    Output bin metadata and additional columns match the input dataframe.

    Parameters
    ----------
    avg_size : int
        Split regions into equal-sized subregions of about this size.
        Specifically, subregions are no larger than 150% of this size, no
        smaller than 75% this size, and the average will approach this size when
        subdividing a large region.
    min_size : int
        Drop any regions smaller than this size.
    verbose : bool
        Print a log message when subdividing a region.

    """
    for row in merge(regions).itertuples(index=False):
        span = row.end - row.start
        if span >= min_size:
            nbins = int(round(span / avg_size)) or 1
            if nbins == 1:
                yield row
            else:
                # Divide the region into equal-sized bins
                bin_size = span / nbins
                bin_start = row.start
                if verbose:
                    label = (row.gene if 'gene' in regions else
                             "%s:%d-%d" % (row.chromosome, row.start, row.end))
                    logging.info("Splitting: {:30} {:7} / {} = {:.2f}"
                                 .format(label, span, nbins, bin_size))
                for i in range(1, nbins):
                    bin_end = row.start + int(i * bin_size)
                    yield row._replace(start=bin_start, end=bin_end)
                    bin_start = bin_end
                yield row._replace(start=bin_start)
