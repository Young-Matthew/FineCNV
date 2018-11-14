from __future__ import absolute_import, division, print_function
from past.builtins import basestring

import collections
import contextlib
import logging
import os
import re
import sys
from datetime import datetime

import pandas as pd
from Bio.File import as_handle
#from gary import GenomicArray as GA
from skgenome import tabio, GenomicArray as GA
import core
# git clone https://github.com/etal/cnvkit.git
# slight modify to corporate with myself

from func import do_target,do_antitarget,batch_write_coverage,batch_run_sample
import parallel
import reference

# x05 PLATFORM
normal_bams = "/data/haplox/users/wangzy/00.Tools/04.TempWork/03.cnvkit.Module/01.wenglihua/10M.css.gDNA.bam"
tumor_bam = "/data/haplox/users/wangzy/00.Tools/04.TempWork/03.cnvkit.Module/01.wenglihua/10M.css.cfDNA.bam"

#normal_bams = "/data/haplox/users/wangzy/03.Project/00.GreenLung/00.451/00.20180907/10.wenglihua/wenglihua_css.gDNA.dedup.bam"
#tumor_bam = "/data/haplox/users/wangzy/03.Project/00.GreenLung/00.451/00.20180907/10.wenglihua/wenglihua_css.cfDNA.dedup.bam"

target_bed = "/x03_haplox/users/chenyr/bin/wesplus/database/target_region//WESplus_gene.bed"
fasta = "/x03_haplox/users/chenyr/bin/wesplus/database/hg19/hg19_wesplus.sort.fa"
access_bed = "/x03_haplox/users/chenyr/bin/wesplus/database/target_region//access-5k-mappable.hg19.bed"

tgt_name_base, _tgt_ext = os.path.splitext(os.path.basename(target_bed))
new_target_fname = tgt_name_base + '.target.bed'

bait_arr = []
annotate = None
short_names = None
target_avg_size = None

bait_arr = tabio.read_auto(target_bed)
target_arr = do_target(bait_arr, annotate, short_names, True,**({'avg_size': target_avg_size} if target_avg_size else {}))
tabio.write(target_arr, new_target_fname, 'bed4')
target_bed = new_target_fname

antitarget_bed = tgt_name_base + '.antitarget.bed'
anti_kwargs = {}
if access_bed:
    anti_kwargs['access'] = tabio.read_auto(access_bed)
anti_arr = do_antitarget(target_arr, **anti_kwargs)
tabio.write(anti_arr, antitarget_bed, "bed4")

# hardest part build reference file
tgt_futures = []
anti_futures = []

nbam = normal_bams
sample_id = core.fbase(nbam)
output_dir = "/data/haplox/users/wangzy/00.Tools/04.TempWork/03.cnvkit.Module/out"
sample_pfx = os.path.join(output_dir, sample_id)
procs_per_cnn = 1
by_count = None

t1 = datetime.now()
print("start to build coverage on both target and antitarget area")
tgt_futures.append(batch_write_coverage(target_bed, nbam, sample_pfx + '.targetcoverage.cnn',by_count, procs_per_cnn))
anti_futures.append(batch_write_coverage(antitarget_bed, nbam,sample_pfx + '.antitargetcoverage.cnn',by_count, procs_per_cnn))
t2 = datetime.now()
t = (t2 - t1).seconds
print("it takes %d seconds to build coverage on target file" % (t))

#target_fnames = [tf.result() for tf in tgt_futures]
#antitarget_fnames = [af.result() for af in anti_futures]
target_fnames = tgt_futures
antitarget_fnames =  anti_futures

# Build reference from *.cnn
print("start to build reference.cnn")
ref_arr = reference.do_reference(target_fnames, antitarget_fnames, fasta, False, None, \
                                                        do_gc=True, \
                                                        do_edge=True,\
                                                        do_rmask=True)
t3 = datetime.now()
t = (t3 - t2).seconds
print("it takes %d seconds to build reference file" % (t))

output_reference = os.path.join(output_dir, "reference.cnn")
core.ensure_path(output_reference)
tabio.write(ref_arr, output_reference)

print("start to use tumor bam file to work to final step")
batch_run_sample(tumor_bam,target_bed, antitarget_bed, output_reference, \
                output_dir, False,False,False, \
                None,False,False,"hybrid",1)
