#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Usage: python summarize_metadata.py
'''

import sys
import pandas as pd
from intervaltree import Interval, IntervalTree


# files need to be specified by users
genome_fai = '/path/to/genome.fa.fai'
gene_model_gtf = '/path/to/gene_model.gtf'  # this should contain ONLY gene models, not TE models
input_file = 'summary_splicing_junctions.tsv'
outfilename = 'summary_splicing_junctions_metadata.tsv'


# parameters
MIN_BLOCK_LEN = 20


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
    exon_intervals = {}
    with open(genome_fai) as infile:
        for line in infile:
            ls = line.split()
            chr = ls[0]
            tree1 = IntervalTree()
            gene_intervals[chr] = tree1
            tree2 = IntervalTree()
            exon_intervals[chr] = tree2
    
    # load genes
    gene_exon = {'gene', 'exon'}
    gene_names = {}
    with open(gene_model_gtf) as infile:
        for line in infile:
            if line[0] == '#':
                continue
            ls = line.split('\t')
            if not ls[2] in gene_exon:
                continue
            chr = ls[0]
            gene_id = get_gene_id(ls[-1])
            start = int(ls[3]) - 1
            end = int(ls[4])
            strand = ls[6]
            if ls[2] == 'gene':
                gene_intervals[chr][start:end] = gene_id
                gene_names[gene_id] = (start, end, strand)
            else:
                exon_intervals[chr][start:end] = gene_id
    return gene_intervals, exon_intervals, gene_names


# load intervals
gene_intervals, exon_intervals, gene_names = load_intervals()

# load positions
df = pd.read_table(input_file, index_col = 0)

data = []
for i in df.index:
    chr, pos, left, right, strand = i.split(':')
    s, e = pos.split('-')
    s = int(s)
    e = int(e)
    if left in gene_names:
        gene = left
        te = right
        gene_pos = s
        te_start = e
        te_end = e + MIN_BLOCK_LEN
    else:
        gene = right
        te = left
        gene_pos = e
        te_start = s - MIN_BLOCK_LEN
        te_end = s
    gene_overlap = gene_intervals[chr][te_start:te_end]
    exon_overlap = exon_intervals[chr][te_start:te_end]
    # judge isin exon
    isin_exon = False
    for i in exon_overlap:
        _, _, gene_name = i
        if gene_name == gene:
            isin_exon = True
            break
    # judge isin gene
    isin_gene = False
    for i in gene_overlap:
        _, _, gene_name = i
        if gene_name == gene:
            isin_gene = True
            break
    # judge isin intron
    isin_intron = False
    if isin_gene and not isin_exon:
        isin_intron = True
    # judge upstream or downstream
    is_up, is_down = False, False
    if not isin_gene:
        gene_start, gene_end, strand = gene_names[gene]
        if strand == '+':
            if te_end < gene_start:
                is_up = True
            else:
                is_down = True
        else:
            if te_start >= gene_end:
                is_up = True
            else:
                is_down = True
    # summarize
    if strand == '+':
        sjdonor = s
        sjaccep = e
        donor = left
        accep = right
    else:
        sjdonor = e
        sjaccep = s
        donor = right
        accep = left
    tmp = [
        chr,
        sjdonor,
        sjaccep,
        donor,
        accep,
        strand,
        isin_exon,
        isin_intron,
        isin_gene,
        is_up,
        is_down,
    ]
    data.append(tmp)

columns = [
    'chr',
    'sj_donor_pos0',
    'sj_acceptor_pos0',
    'donor_annot',
    'acceptor_annot',
    'strand',
    'TE_isin_exon',
    'TE_isin_intron',
    'TE_isin_gene',
    'TE_upstream',
    'TE_downstream',
]
df = pd.DataFrame(data, index = df.index, columns = columns)

df.to_csv(outfilename, sep = '\t')
print('Shape of metadata table = (%d, %d)' % df.shape)
