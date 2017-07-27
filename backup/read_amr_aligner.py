#!/usr/bin/env python2.6
import sys
import os
import time
import random
import amr
import smatch
from smatch import *
import amr_graph
from amr_graph import *
def delete_unaligned(s, alignment):
    return

def get_concept(amr_hg, concept_index):
    return

if __name__ == "__main__":
    amr_file = sys.argv[1]
    f = open(amr_file, 'r')
    for i in xrange(4):
        f.readline()
    line = f.readline().strip()
    alignments = line.split()
    alignments = alignments[2:-5]
    alignments = [tuple(align.split('|')) for align in alignments]
    starts = []
    ends = []
    fragments = []
    for indexes in alignments:
        s_index = indexes[0]
        f_index = indexes[1]
        interval = s_index.split('-')
        starts.append(interval[0])
        ends.append(interval[1])
        frag = f_index.split('+')
        fragments.append(frag)

    f1 = AMRFragment(bitarray('10'), bitarray('0100'))
    f2 = AMRFragment(bitarray('10'), bitarray('0100'))
    print f1 == f2
    s = set()
    for i in xrange(3):
        s.add(f1)
        s.add(f2)
    print len(s)

    while True:
        curr_amr = get_amr_line(f)
        #curr_amr_hg = amr.AMR.parse_AMR_line(curr_amr, fragments)
        #curr_amr_hg.__str__()
        hg = AMRGraph(curr_amr)
        print str(hg)
        frag = hg.retrieve_fragment('0.0.1.0+0.0.1.0.0+0.0.1.0.0.0')
        print str(frag)
        #print str(hg)
        #hg.print_info()
        #curr_amr_hg.__str__()
        break
