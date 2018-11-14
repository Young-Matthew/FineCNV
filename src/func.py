from __future__ import print_function, absolute_import, division
from skgenome import tabio

from builtins import range
from builtins import zip
from past.builtins import basestring

import logging
import pandas as pd
import numpy as np
from merge import merge
import re

import math
import os.path
import time
from concurrent import futures
from Bio._py3k import StringIO
from cnary import CopyNumArray as CNA
import samutil
import core
import parallel
import pysam
import scatter, segmentation
import fix
from cmdutil import read_cna


#from params
"""Hard-coded parameters for CNVkit. These should not change between runs."""
# Filter thresholds used in constructing the reference (log2 scale)
MIN_REF_COVERAGE = -5.0
MAX_REF_SPREAD = 1.0
NULL_LOG2_COVERAGE = -20.0

# Fragment size for paired-end reads
INSERT_SIZE = 250

# Target/bin names that are not meaningful gene names
# (In some UCSF panels, "CGH" probes denote selected intergenic regions)
IGNORE_GENE_NAMES = ("-", ".", "CGH")
ANTITARGET_NAME = "Antitarget"
ANTITARGET_ALIASES = (ANTITARGET_NAME, "Background")


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


def do_antitarget(targets, access=None, avg_bin_size=150000,
                  min_bin_size=None):
    """Derive off-targt ("antitarget") bins from target regions."""
    if not min_bin_size:
        min_bin_size = 2 * int(avg_bin_size * (2 ** MIN_REF_COVERAGE))
        # min_bin_size = 9374
        # MIN_REF_COVERAGE = -5.0
        # MAX_REF_SPREAD = 1.0
        # NULL_LOG2_COVERAGE = -20.0
        # INSERT_SIZE = 250
    return get_antitargets(targets, access, avg_bin_size, min_bin_size)


def get_antitargets(targets, accessible, avg_bin_size, min_bin_size):
    """Generate antitarget intervals between/around target intervals.

    Procedure:

    - Invert target intervals
    - Subtract the inverted targets from accessible regions
    - For each of the resulting regions:

        - Shrink by a fixed margin on each end
        - If it's smaller than min_bin_size, skip
        - Divide into equal-size (region_size/avg_bin_size) portions
        - Emit the (chrom, start, end) coords of each portion
    """
    if accessible:
        # Chromosomes' accessible sequence regions are given -- use them
        accessible = drop_noncanonical_contigs(accessible, targets)
    else:
        # Chromosome accessible sequence regions not known -- use heuristics
        # (chromosome length is endpoint of last probe; skip initial
        # <magic number> of bases that are probably telomeric)
        TELOMERE_SIZE = 150000
        accessible = guess_chromosome_regions(targets, TELOMERE_SIZE)
    pad_size = 2 * INSERT_SIZE
    # 2 * 250
    bg_arr = (accessible.resize_ranges(-pad_size)
              .subtract(targets.resize_ranges(pad_size))
              .subdivide(avg_bin_size, min_bin_size))
    bg_arr['gene'] = ANTITARGET_NAME
    return bg_arr


def guess_chromosome_regions(targets, telomere_size):
    """Determine (minimum) chromosome lengths from target coordinates."""
    endpoints = [subarr.end.iat[-1] for _c, subarr in targets.by_chromosome()]
    whole_chroms = GA.from_columns({
        'chromosome': targets.chromosome.drop_duplicates(),
        'start': telomere_size,
        'end': endpoints})
    return whole_chroms


def drop_noncanonical_contigs(accessible, targets, verbose=True):
    """Drop contigs that are not targeted or canonical chromosomes.
    #canonical : standard
    Antitargets will be binned over chromosomes that:

    - Appear in the sequencing-accessible regions of the reference genome
      sequence, and
    - Contain at least one targeted region, or
    - Are named like a canonical chromosome (1-22,X,Y for human)

    This allows antitarget binning to pick up canonical chromosomes that do not
    contain any targets, as well as non-canonical or oddly named chromosomes
    that were targeted.
    """
    # TODO - generalize: (garr, by="name", verbose=True):

    access_chroms, target_chroms = compare_chrom_names(accessible, targets)
    # Filter out untargeted alternative contigs and mitochondria
    untgt_chroms = access_chroms - target_chroms
    # Autosomes typically have numeric names, allosomes are X and Y
    if any(is_canonical_contig_name(c) for c in target_chroms):
        chroms_to_skip = [c for c in untgt_chroms
                            if not is_canonical_contig_name(c)]
        # chrM is inside our access and target file,therefor won't be excluded
        # and there is no warning
    else:
        # Alternative contigs have longer names -- skip them
        max_tgt_chr_name_len = max(map(len, target_chroms))
        chroms_to_skip = [c for c in untgt_chroms
                            if len(c) > max_tgt_chr_name_len]
    if chroms_to_skip:
        logging.info("Skipping untargeted chromosomes %s",
                     ' '.join(sorted(chroms_to_skip)))
        skip_idx = accessible.chromosome.isin(chroms_to_skip)
        accessible = accessible[~skip_idx]
    return accessible

def compare_chrom_names(a_regions, b_regions):
    a_chroms = set(a_regions.chromosome.unique())
    b_chroms = set(b_regions.chromosome.unique())
    if a_chroms and a_chroms.isdisjoint(b_chroms):
        # isdisjoint return true when a and b contains no common member
        msg = "Chromosome names do not match between files"
        a_fname = a_regions.meta.get('filename')
        b_fname = b_regions.meta.get('filename')
        if a_fname and b_fname:
            msg += " {} and {}".format(a_fname, b_fname)
        msg += ": {} vs. {}".format(', '.join(map(repr, sorted(a_chroms)[:3])),
                                    ', '.join(map(repr, sorted(b_chroms)[:3])))
        raise ValueError(msg)
    return a_chroms, b_chroms

# CNVkit's original inclusion regex
re_canonical = re.compile(r"(chr)?(\d+|[XYxy])$")
# goleft indexcov's exclusion regex
re_noncanonical = re.compile("|".join((r"^chrEBV$",
                                       r"^NC|_random$",
                                       r"Un_",
                                       r"^HLA\-",
                                       r"_alt$",
                                       r"hap\d$",
                                       r"chrM",
                                       r"MT")))


def is_canonical_contig_name(name):
    #return bool(re_canonical.match(name))
    return not re_noncanonical.search(name)


def batch_write_coverage(bed_fname, bam_fname, out_fname, by_count, processes):
    """Run coverage on one sample, write to file."""
    cnarr = do_coverage(bed_fname, bam_fname, by_count, 0, processes)
    tabio.write(cnarr, out_fname)
    return out_fname

def do_coverage(bed_fname, bam_fname, by_count=False, min_mapq=0, processes=1):
    """Calculate coverage in the given regions from BAM read depths."""
    if not samutil.ensure_bam_sorted(bam_fname):
        raise RuntimeError("BAM file %s must be sorted by coordinates"
                            % bam_fname)
    samutil.ensure_bam_index(bam_fname)
    # ENH: count importers.TOO_MANY_NO_COVERAGE & warn
    cnarr = interval_coverages(bed_fname, bam_fname, by_count, min_mapq,
                               processes)
    return cnarr

def ensure_bam_sorted(bam_fname, by_name=False, span=50):
    """Test if the reads in a BAM file are sorted as expected.

    by_name=True: reads are expected to be sorted by query name. Consecutive
    read IDs are in alphabetical order, and read pairs appear together.

    by_name=False: reads are sorted by position. Consecutive reads have
    increasing position.
    """
    if by_name:
        # Compare read IDs
        def out_of_order(read, prev):
            return not (prev is None or
                        prev.qname <= read.qname)
    else:
        # Compare read locations
        def out_of_order(read, prev):
            return not (prev is None or
                        read.tid != prev.tid or
                        prev.pos <= read.pos)

    # ENH - repeat at 50%, ~99% through the BAM
    bam = pysam.Samfile(bam_fname, 'rb')
    last_read = None
    for read in islice(bam, span):
        if out_of_order(read, last_read):
            return False
        last_read = read
    bam.close()
    return True

def ensure_bam_index(bam_fname):
    """Ensure a BAM file is indexed, to enable fast traversal & lookup.

    For MySample.bam, samtools will look for an index in these files, in order:

    - MySample.bam.bai
    - MySample.bai
    """
    if os.path.isfile(bam_fname + '.bai'):
        # MySample.bam.bai
        bai_fname = bam_fname + '.bai'
    else:
        # MySample.bai
        bai_fname = bam_fname[:-1] + 'i'
    if not is_newer_than(bai_fname, bam_fname):
        logging.info("Indexing BAM file %s", bam_fname)
        pysam.index(bam_fname)
        bai_fname = bam_fname + '.bai'
    assert os.path.isfile(bai_fname), \
            "Failed to generate index " + bai_fname
    return bai_fname

def is_newer_than(target_fname, orig_fname):
    """Compare file modification times."""
    if not os.path.isfile(target_fname):
        return False
    return (os.stat(target_fname).st_mtime >= os.stat(orig_fname).st_mtime)


def interval_coverages(bed_fname, bam_fname, by_count, min_mapq, processes):
    """Calculate log2 coverages in the BAM file at each interval."""
    meta = {'sample_id': core.fbase(bam_fname)}
    start_time = time.time()

    # Skip processing if the BED file is empty
    with open(bed_fname) as bed_handle:
        for line in bed_handle:
            if line.strip():
                break
        else:
            logging.info("Skip processing %s with empty regions file %s",
                         os.path.basename(bam_fname), bed_fname)
            return CNA.from_rows([], meta_dict=meta)

    # Calculate average read depth in each bin
    if by_count:
        results = interval_coverages_count(bed_fname, bam_fname, min_mapq,
                                           processes)
        read_counts, cna_rows = zip(*results)
        read_counts = pd.Series(read_counts)
        cnarr = CNA.from_rows(list(cna_rows),
                              columns=CNA._required_columns + ('depth',),
                              meta_dict=meta)
    else:
        table = interval_coverages_pileup(bed_fname, bam_fname, min_mapq,
                                          processes)
        read_len = samutil.get_read_length(bam_fname)
        read_counts = table['basecount'] / read_len
        # read_len = 143
        table = table.drop('basecount', axis=1)
        cnarr = CNA(table, meta)

    # Log some stats
    tot_time = time.time() - start_time
    tot_reads = read_counts.sum()
    logging.info("Time: %.3f seconds (%d reads/sec, %s bins/sec)",
                 tot_time,
                 int(round(tot_reads / tot_time, 0)),
                 int(round(len(read_counts) / tot_time, 0)))
    logging.info("Summary: #bins=%d, #reads=%d, "
                 "mean=%.4f, min=%s, max=%s ",
                 len(read_counts),
                 tot_reads,
                 (tot_reads / len(read_counts)),
                 read_counts.min(),
                 read_counts.max())
    tot_mapped_reads = samutil.bam_total_reads(bam_fname)
    # number of all reads in bam file, 999973 in /haplox/users/wangzy/02.caizhenling/02.Faster.MrBam/00.data/1M.Caizhenling_css.cfDNA.dedup.bam
    if tot_mapped_reads:
        logging.info("Percent reads in regions: %.3f (of %d mapped)",
                     100. * tot_reads / tot_mapped_reads,
                     tot_mapped_reads)
    else:
        logging.info("(Couldn't calculate total number of mapped reads)")

    return cnarr


def interval_coverages_pileup(bed_fname, bam_fname, min_mapq, procs=1):
    """Calculate log2 coverages in the BAM file at each interval."""
    logging.info("Processing reads in %s", os.path.basename(bam_fname))
    if procs == 1:
        table = bedcov(bed_fname, bam_fname, min_mapq)
    else:
        chunks = []
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bed_chunk, bam_fname, min_mapq)
                         for bed_chunk in to_chunks(bed_fname))
            for bed_chunk_fname, table in pool.map(_bedcov, args_iter):
                chunks.append(table)
                rm(bed_chunk_fname)
        table = pd.concat(chunks, ignore_index=True)
    # Fill in CNA required columns
    if 'gene' in table:
        table['gene'] = table['gene'].fillna('-')
        # replace NaN in panda-produced file as -
    else:
        table['gene'] = '-'
    # User-supplied bins might be zero-width or reversed -- skip those
    spans = table.end - table.start
    # table's start position and end position
    # our target bed contains no zero-wideth or reversed line
    ok_idx = (spans > 0)
    table = table.assign(depth=0, log2=NULL_LOG2_COVERAGE)
    # NULL_LOG2_COVERAGE = -20.0
    # adding two columns to table 
    table.loc[ok_idx, 'depth'] = (table.loc[ok_idx, 'basecount']
                                  / spans[ok_idx])
    ok_idx = (table['depth'] > 0)
    table.loc[ok_idx, 'log2'] = np.log2(table.loc[ok_idx, 'depth'])
    return table

def bedcov(bed_fname, bam_fname, min_mapq):
    """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.

    i.e. mean pileup depth across each region.
    """
    # Count bases in each region; exclude low-MAPQ reads
    cmd = [bed_fname, bam_fname]
    if min_mapq and min_mapq > 0:
        cmd.extend(['-Q', bytes(min_mapq)])
    try:
        raw = pysam.bedcov(*cmd, split_lines=False)
        # only one line, waiting for tranferred to table format
        # adding number of bases cover in this target region: basecount in our pipeline
    except pysam.SamtoolsError as exc:
        raise ValueError("Failed processing %r coverages in %r regions. "
                         "PySAM error: %s" % (bam_fname, bed_fname, exc))
    if not raw:
        raise ValueError("BED file %r chromosome names don't match any in "
                         "BAM file %r" % (bed_fname, bam_fname))
    columns = detect_bedcov_columns(raw)
    table = pd.read_table(StringIO(raw), names=columns, usecols=columns)
    return table

def detect_bedcov_columns(text):
    """Determine which 'bedcov' output columns to keep.

    Format is the input BED plus a final appended column with the count of
    basepairs mapped within each row's region. The input BED might have 3
    columns (regions without names), 4 (named regions), or more (arbitrary
    columns after 'gene').
    """
    firstline = text[:text.index('\n')]
    tabcount = firstline.count('\t')
    if tabcount < 3:
        raise RuntimeError("Bad line from bedcov:\n%r" % firstline)
    if tabcount == 3:
        return ['chromosome', 'start', 'end', 'basecount']
    if tabcount == 4:
        return ['chromosome', 'start', 'end', 'gene', 'basecount']
        #This is the format of our result
    # Input BED has arbitrary columns after 'gene' -- ignore them
    fillers = ["_%d" % i for i in range(1, tabcount - 3)]
    return ['chromosome', 'start', 'end', 'gene'] + fillers + ['basecount']


def get_read_length(bam, span=1000):
    """Get (median) read length from first few reads in a BAM file.
    
    Illumina reads all have the same length; other sequencers might not.

    Parameters
    ----------
    bam : str or pysam.Samfile
        Filename or pysam-opened BAM file.
    n : int
        Number of reads used to calculate median read length.
    """
    was_open = False
    if isinstance(bam, basestring):
        bam = pysam.Samfile(bam, 'rb')
    else:
        was_open = True
    lengths = [read.query_length for read in islice(bam, span)
               if read.query_length > 0]
    if was_open:
        bam.seek(0)
    else:
        bam.close()
    return np.median(lengths)
    # 143


def bam_total_reads(bam_fname):
    """Count the total number of mapped reads in a BAM file.

    Uses the BAM index to do this quickly.
    """
    table = idxstats(bam_fname, drop_unmapped=True)
    return table.mapped.sum()

def idxstats(bam_fname, drop_unmapped=False):
    """Get chromosome names, lengths, and number of mapped/unmapped reads.

    Use the BAM index (.bai) to get the number of reads and size of each
    chromosome. Contigs with no mapped reads are skipped.
    """
    handle = StringIO(pysam.idxstats(bam_fname, split_lines=False))
    table = pd.read_table(handle, header=None,
                          names=['chromosome', 'length', 'mapped', 'unmapped'])
# table:
#    chromosome     length  mapped  unmapped
#0        chr1  249250621  999973         0
#1        chr2  243199373       0         0
#2        chr3  198022430       0         0

    if drop_unmapped:
        table = table[table.mapped != 0].drop('unmapped', axis=1)
    return table


def batch_run_sample(bam_fname, target_bed, antitarget_bed, ref_fname,
                     output_dir, male_reference, plot_scatter, plot_diagram,
                     rlibpath, by_count, skip_low, method, processes):
    """Run the pipeline on one BAM file."""
    # ENH - return probes, segments (cnarr, segarr)
    logging.info("Running the CNVkit pipeline on %s ...", bam_fname)
    sample_id = core.fbase(bam_fname)
    sample_pfx = os.path.join(output_dir, sample_id)

    raw_tgt = do_coverage(target_bed, bam_fname, by_count, 0,
                                   processes)
    tabio.write(raw_tgt, sample_pfx + '.targetcoverage.cnn')

    raw_anti = do_coverage(antitarget_bed, bam_fname, by_count, 0,
                                    processes)
    tabio.write(raw_anti, sample_pfx + '.antitargetcoverage.cnn')

    cnarr = fix.do_fix(raw_tgt, raw_anti, read_cna(ref_fname),
                       do_gc=True, do_edge=(method == "hybrid"), do_rmask=True)
    tabio.write(cnarr, sample_pfx + '.cnr')

    logging.info("Segmenting %s.cnr ...", sample_pfx)
    segments = segmentation.do_segmentation(cnarr, 'cbs',
                                            rlibpath=rlibpath,
                                            skip_low=skip_low,
                                            processes=processes,
                                            **({'threshold': 1e-6}
                                               if method == 'wgs'
                                               else {}))
    tabio.write(segments, sample_pfx + '.cns')

    if plot_scatter:
        scatter.do_scatter(cnarr, segments)
        pyplot.savefig(sample_pfx + '-scatter.pdf', format='pdf',
                       bbox_inches="tight")
        logging.info("Wrote %s-scatter.pdf", sample_pfx)

    if plot_diagram:
        is_xx = cnarr.guess_xx(male_reference)
        outfname = sample_pfx + '-diagram.pdf'
        diagram.create_diagram(cnarr.shift_xx(male_reference, is_xx),
                               segments.shift_xx(male_reference, is_xx),
                               0.5, 3, outfname)
        logging.info("Wrote %s", outfname)
