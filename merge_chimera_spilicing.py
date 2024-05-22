#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Usage: python merge_chimera_spilicing.py
'''

import sys,glob
import numpy as np
import pandas as pd
import collections


# Users need to specify the pass to '[sample_name].sj.txt'
# Here, you can use a wild card '*'
sj_file_path = '/path/to/*.sj.txt'

# filenames
outfilename = 'summary_splicing_junctions.tsv'



# file check
files = glob.glob(sj_file_path)
files = sorted(files)
if len(files) == 0:
    print('"[sample_name].sj.txt" was not found.', file = sys.stderr)
    exit(1)
else:
    print('n=%d files will be merged...' % len(files))

# main
df = pd.DataFrame()
for f in files:
    sample = f.split('/')[-1].split('.')[0]
    counter = collections.Counter()
    with open(f) as infile:
        for line in infile:
            ls = line.split()
            s = ls[5]
            e = ls[6]
            junction = '%s:%s:%s' % (s, e, ls[3])
            sj = '%s:%s-%s:%s' % (*ls[:3], junction)
            counter[sj] += 1
    sjs = sorted(counter.keys())
    vals = np.array([ counter[s] for s in sjs ]).reshape(-1, 1)
    tmp = pd.DataFrame(vals, index = sjs, columns = [sample])
    df = df.join(tmp, how = 'outer')

df = df.fillna(0)
df = df.astype(int)
df.index.name = 'junction'

df.to_csv(outfilename, sep = '\t')
print('Shape of merged count table = (%d, %d)' % df.shape)
