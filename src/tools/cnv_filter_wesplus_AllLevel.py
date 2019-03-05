#!/usr/bin/env python3

from os.path import splitext
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import logging
import re
import sys

class SingleMetavarHelpFormatter(RawDescriptionHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                parts.extend(action.option_strings)
                parts[-1] += ' %s' % args_string
            return ', '.join(parts)

def parse_args():
    description='Description:\nUsage:\n\
                -l  input cnvgene list file force\n\
                -c  input cnvkit cns file force\n\
                -o  out file force\n\
                -g  gff file\n'
                #-h  Help document\n'

    parser=ArgumentParser(formatter_class=SingleMetavarHelpFormatter,description=description)
    parser.add_argument('-l',)
    parser.add_argument('-c',)
    parser.add_argument('-o',)
    #parser.add_argument('-h')
    parser.add_argument('-g','--gff',help='gff-format file contains location information about exons')
    parser.add_argument('-b','--band',help='tab-delimited file contains cytoband info about homo sapiens')
    return  parser.parse_args()

description="Description:\nUsage:\n\
                -l  input cnvgene list file force\n\
                -c  input cnvkit cns file force\n\
                -o  out file force\n\
                -g  --gff gff file [/haplox/users/wangzy/Software/hg19.gff]\n\
                -b  --band cytoband file[/haplox/users/wangzy/Software/cytoBand.txt]\n"

def try_append(d,k,v):
    if k in d:
        d[k].append(v)
    else:
        d[k]=[v]

if __name__ == '__main__':

    o = parse_args()
    if o.l is None or o.c is None or o.o is None:
        print(description)
        exit(-1)

    if (o.gff is not None and o.band is None) or (o.gff is None and o.band is not None):
        print('--gff --band must coexist')
        exit(-1)

    gene=open(o.l,'r')
    #gene=open('final.cnv.genelist','r')
    gene_report={}
    for line in gene:
        line=line.rstrip()
        gene_report[line]=1
    gene.close()

    exon_info={}
    filter=False
    current_gene=''
    current_tranc=''
    current_chr=''
    exonID=1
    Nocnv=True
    header=True
    exonOutputFile=''
    bandChange=False
    bandheader=True

    if o.gff is not None:
        gff=open(o.gff,'r')
        for line in gff:

            line = line.rstrip().split('\t')
            if len(line)==1:
                continue

            if line[-1].startswith('Parent='):
                if filter:
                    continue
                else:
                    if line[2] != 'intron':
                        location = [int(line[3]), int(line[4])]
                        if line[2] == '3_UTR' or line[2] == '5_UTR':
                            id=line[2]
                        else:
                            id='exon' + str(exonID)
                            exonID+=1
                        exon_info[current_gene][current_chr][current_tranc][id] = location

            elif line[-1].startswith('ID='):
                match=re.search(r'ID=(.*); name=(.*);',line[-1])
                if match:
                    current_gene=match.group(2)
                    current_tranc=match.group(1)
                    current_chr=line[0]

                    if current_gene in gene_report:
                        filter=False
                        if current_gene not in exon_info:
                            exon_info[current_gene]={}
                        if current_chr not in exon_info[current_gene]:
                            exon_info[current_gene][current_chr]={}
                            # chrX and chrY might contains homozygous gene

                        exon_info[current_gene][current_chr][current_tranc]={}
                        exonID=1

                    else:
                        filter=True
                        continue
                else:
                    print('unnormal gene name')
                    print(line)
                    exit(-1)

            else:
                print('line format unnormal')
                print(line)
                exit(-1)

    cytoband={}
    if o.band is not None:
        band=open(o.band,'r')
        for line in band:
            info=line.rstrip().split()
            if info[0] not in cytoband:
                cytoband[info[0]]={}
            cytoband[info[0]][info[3]]=[int(info[1]),int(info[2])]

    cns=open(o.c,'r')
    basename, suffix = splitext(o.o)
    geneOutFile=basename +'.Gene'+suffix
    out=open(geneOutFile,'w')
    if o.gff:
        exonOutputFile=basename +'.Exon'+ suffix
        bandOutputFile = basename + '.Band' + suffix
        exonOut=open(exonOutputFile,'w')
        bandOut=open(bandOutputFile,'w')
    #cns=open('yangkai.tumor.dedupped.cns','r')
    #out=open('test.out.py','w')
    #print('%d genes stored ' % (len(exon_info.keys())))
    cnvkit={}
    # exon_info[current_gene][current_chr][current_tranc][id] = location
    for line in cns:
        info=line.rstrip().split('\t')
        genes=info[3].split(',')

        for gene in genes:
            exonName = ''
            exonLocation = ''
            if gene == 'AR' or gene == 'ARTX' or \
            (gene in gene_report and (float(info[4]) <= -0.9 or float(info[4]) >=0.5)):
                # Special genes on chrX and chrY

                for chr in exon_info[gene]:
                    exonName+=str(chr)+' '
                    for trans in exon_info[gene][chr]:
                        exonName+=str(trans)
                        for exon in exon_info[gene][chr][trans]:
                            if int(info[1]) <= exon_info[gene][chr][trans][exon][0] and \
                                int(info[2]) >= exon_info[gene][chr][trans][exon][1]:
                                exonName+=':' + str(exon)

                            elif not( int(info[2]) < exon_info[gene][chr][trans][exon][0] or \
                                int(info[1]) > exon_info[gene][chr][trans][exon][1]):
                                print("unusual exon location,location in gff file:\
                                      %d %d\nLocation in cnvkit: %d %d \
                                       %(exon_info[gene][chr][trans][exon][0],\
                                       exon_info[gene][chr][trans][exon][1],\
                                       info[1],info[2])")
                        exonName+=' '

                copynumber=str(round(2*2**float(info[4]),2))
                log2ratio=round(float(info[4]),2)
                type= 'loss' if log2ratio < 0 else 'gain'
                exonName='None' if len(exonName) == 0 else exonName
                if header:
                    print(*['Gene','Log2ratio','Copynumber','Type','Exon'],sep='\t',file=exonOut)
                    header=False
                print(*[gene,str(log2ratio),copynumber,type,exonName],sep='\t',file=exonOut)
                Nocnv=False

            #maintaing the same result from chenyr's script:cnv_filter_wesplus.pl
            ########################################################################
            if gene == 'AR' or gene == 'ARTX':
                cnvkit[gene]=info[4]
            if (gene in gene_report and (float(info[4]) <= -0.9 or float(info[4]) >= 0.5)):
                cnvkit[gene]=info[4] if (gene not in cnvkit or float(info[4]) > float(cnvkit[gene])) else cnvkit[gene]
            ##########################################################################

        if info[0] in cytoband and not -0.9<float(info[4])<0.5:
            for eachBand,pos in cytoband[info[0]].items():
                if int(info[1]) <= pos[0] and pos[1] <=int(info[2]):
                    copynumber = str(round(2 * 2 ** float(info[4]), 2))
                    type = 'loss' if float(info[4]) < 0 else 'gain'
                    if bandheader:
                        print(*['Chr','Band','Log2ratio', 'Copynumber', 'Type','StartInBand','EndInBand'],\
                              sep='\t', file=bandOut)
                    print(*[info[0],eachBand,str(round(float(info[4]), 2)), copynumber, type,\
                            str(pos[0]), str(pos[1]),], sep='\t', file=bandOut)
                    bandChange=True
                    bandheader=False
                    break
    if Nocnv:
        print('There is no cnv',file=exonOut)
    if not bandChange:
        print('There is no cnv on band level',file=bandOut)

    if len(cnvkit.keys()):
        print(*['Gene', 'Log2ratio', 'Copynumber', 'Type'], sep='\t', file=out)
        for gene,log2ratio in cnvkit.items():
            copynumber = str(round(2 * 2 ** float(log2ratio),2))
            type = 'loss' if float(log2ratio)< 0 else 'gain'
            print(*[gene, str(round(float(log2ratio),2)), copynumber, type], sep='\t', file=out)
    else:
        print('There is no cnv', file=out)

    cns.close()
    out.close()
    exonOut.close()
    bandOut.close()