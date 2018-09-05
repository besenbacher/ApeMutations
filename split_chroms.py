#!/usr/bin/env python
from bx.seq.twobit import *
import sys

#TwoBit File:
t = TwoBitFile(open(sys.argv[1]))

#minimum number of consecutive N's before split
minN = int(sys.argv[2])
minRange = int(sys.argv[3])

for chrom in t:
    if 'Un' in chrom:
        continue
    clen = len(t[chrom])
    last = 0
    i = minRange
    while i<clen:
        if t[chrom][i-1:i].upper() == 'N':
            start = i
            end = i
            while start > 0 and t[chrom][start-1:start].upper() == 'N':
                start -= 1
            while end < clen and t[chrom][end-1:end].upper() == 'N':
                end += 1
            if end - start > minN:
                new_last =  start + ((end-start)/2)
                if end == clen:
                    new_last = clen
                print chrom, last, new_last
                last = new_last
                i = end + minRange - 1
        i += minN
    if last < clen:
        print chrom, last, clen
