#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Usage: python detect_chimera_splicing.py input.markdup.namesorted.bam
Input: name of the input bam file must be '/path/to/[sample_name].markdup.namesorted.bam'
Output: this script outputs '/path/to/[sample_name].sj.txt'
'''

import sys
import pysam
from intervaltree import Interval, IntervalTree


# files need to be specified by users
genome_fai = '/path/to/genome.fa.fai'
gene_model_gtf = '/path/to/gene_model.gtf'  # this should contain ONLY gene models, not TE models
te_model_gtf = '/path/to/TE_model.gtf'  # this should contain ONLY TE models, not gene models


# filenames
input_bam = sys.argv[1]
path_base = input_bam.rsplit('/', 1)[0]
sample_name = input_bam.split('/')[-1].replace('.markdup.namesorted.bam', '')
outfilename = '%s/%s.sj.txt' % (path_base, sample_name)


# parameters
MIN_BLOCK_LEN = 20
MIN_INTERSECT_LEN = 20




def strand_to_bool(s):
    if s == '+':
        return True
    return False

def get_gene_id(attr):
    for i in attr.split(';'):
        key, value = i.strip().split(' ')
        if key == 'gene_id':
            return value.replace('"', '')
    print('"gene_id" was not found in the input GTF file.', file = sys.stderr)
    exit(1)

def load_intervals():
    # make interval trees
    gene_intervals = {}
    te_intervals = {}
    with open(genome_fai) as infile:
        for line in infile:
            ls = line.split()
            chr = ls[0]
            tree1 = IntervalTree()
            gene_intervals[chr] = tree1
            tree2 = IntervalTree()
            te_intervals[chr] = tree2
    
    # load genes
    with open(gene_model_gtf) as infile:
        for line in infile:
            if line[0] == '#':
                continue
            ls = line.split('\t')
            if not ls[2] == 'exon':
                continue
            chr = ls[0]
            gene_id = get_gene_id(ls[-1])
            start = int(ls[3]) - 1
            end = int(ls[4])
            strand = strand_to_bool(ls[6])
            info = (gene_id, strand)
            gene_intervals[chr][start:end] = info

    # load TEs
    with open(te_model_gtf) as infile:
        for line in infile:
            if line[0] == '#':
                continue
            ls = line.split('\t')
            if not ls[2] == 'exon':
                continue
            chr = ls[0]
            gene_id = get_gene_id(ls[-1])
            start = int(ls[3]) - 1
            end = int(ls[4])
            strand = strand_to_bool(ls[6])
            info = (gene_id, strand)
            te_intervals[chr][start:end] = info
    return gene_intervals, te_intervals


# load intervals
gene_intervals, te_intervals = load_intervals()


# store read info
class Read:
    def __init__(self, read):
        self.chr = read.reference_name
        self.blocks = read.get_blocks()
        self.is_forward = read.is_forward
        self.seq = read.query_sequence


def process_reads(reads):
    if len(reads) == 0:
        return
    sjs = []
    for ri, read in enumerate(reads):
        # initial checks
        if len(read.blocks) <= 1:
            continue
        # check for gene and TE overlap
        chr = read.chr
        gi = gene_intervals[chr]
        ti = te_intervals[chr]
        block_overlaps = []
        for s, e in read.blocks:
            if e - s < MIN_BLOCK_LEN:
                block_overlaps.append(None)
                continue
            gos = gi[s:e]
            tos = ti[s:e]
            if len(gos) >= 1 and len(tos) == 0:
                gos = sorted(gos)
                overlaps = ['g', []]
                for go in gos:
                    gs, ge, ginfo = go
                    intersect_len = min(ge, e) - max(gs, s)
                    if intersect_len < MIN_INTERSECT_LEN:
                        continue
                    gene_id, strand = ginfo
                    if strand != read.is_forward:
                        continue
                    info = [gene_id, chr, s, e, ri]
                    overlaps[1].append(info)
                if len(overlaps[1]) >= 1:
                    block_overlaps.append(overlaps)
                else:
                    block_overlaps.append(None)
            elif len(tos) >= 1:
                tos = sorted(tos)
                overlaps = ['t', []]
                for to in tos:
                    ts, te, tinfo = to
                    intersect_len = min(te, e) - max(ts, s)
                    if intersect_len < MIN_INTERSECT_LEN:
                        continue
                    gene_id, strand = tinfo
                    info = [gene_id, chr, s, e, ri]
                    overlaps[1].append(info)
                if len(overlaps[1]) >= 1:
                    block_overlaps.append(overlaps)
                else:
                    block_overlaps.append(None)
            else:
                block_overlaps.append(None)
        # process overlaps
        for i in range(len(block_overlaps) - 1):
            l = block_overlaps[i]
            if l == None:
                continue
            elif l[0] == 'g':
                nex = block_overlaps[i+1]
                if nex == None:
                    continue
                elif nex[0] == 'g':
                    continue
                elif nex[0] == 't':
                    # gene, most right
                    gene_id, chr, s, e, ri = l[1][-1]
                    seq = reads[ri].seq
                    strand = '-'
                    if reads[ri].is_forward:
                        strand = '+'
                    # TE, most left
                    te_id, tchr, ts, te, ri = nex[1][0]
                    # junction
                    sj = [chr, e, ts, strand, 'gene:te', gene_id, te_id, seq]
                    sjs.append(sj)
            elif l[0] == 't':
                nex = block_overlaps[i+1]
                if nex == None:
                    continue
                elif nex[0] == 'g':
                    # TE, most right
                    te_id, tchr, ts, te, ri = l[1][-1]
                    # gene, most left
                    gene_id, chr, s, e, ri = nex[1][0]
                    seq = reads[ri].seq
                    strand = '-'
                    if reads[ri].is_forward:
                        strand = '+'
                    # junction
                    sj = [chr, te, s, strand, 'te:gene', te_id, gene_id, seq]
                    sjs.append(sj)
                elif nex[0] == 't':
                    continue
    # format
    if len(sjs) == 0:
        return
    lines = []
    for sj in sjs:
        tmp = '\t'.join([ str(i) for i in sj ])
        lines.append(tmp + '\n')
    return ''.join(lines)


# read bam
outfile = open(outfilename, 'w')
samfile = pysam.AlignmentFile(input_bam, 'rb')
read1s = []
read2s = []
current_rname = 'any'
n_processed = 0
for read in samfile:
    if read.is_unmapped:
        continue
    if read.is_duplicate:
        continue
    nh = read.get_tag('NH')
    if nh >= 2:
        continue
    rname = read.query_name
    if not rname == current_rname:
        ret1 = process_reads(read1s)
        ret2 = process_reads(read2s)
        for ret in (ret1, ret2):
            if ret != None:
                outfile.write(ret)
        read1s = []
        read2s = []
        current_rname = rname
    r = Read(read)
    if read.is_read1:
        read1s.append(r)
    else:
        read2s.append(r)
    n_processed += 1
    if n_processed % 100_000 == 0:
        print(n_processed, 'reads processed...')
ret1 = process_reads(read1s)
ret2 = process_reads(read2s)
for ret in (ret1, ret2):
    if ret != None:
        outfile.write(ret)
outfile.close()
samfile.close()
print(n_processed, 'reads processed.')
