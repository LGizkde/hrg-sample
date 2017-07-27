#!/usr/bin/python
import sys
import time
import pickle
import os
import cPickle
import alignment
import hypergraph
from fragment_hypergraph import FragmentHGNode, FragmentHGEdge
import smatch
from smatch import get_amr_line
import amr_graph
from amr_graph import *
import logger
import gflags
from data_divider import main
import HRGSample
from HRGSample import *
from rule import Rule
def main(argv):
    cycle_set = set()
    cycle_file = argv[1]
    with open(cycle_file, 'r') as f:
        for line in f:
            cycle_set.add(int(line.strip()))

    amr_file = argv[2]
    align_file = argv[3]
    sent_file = argv[4]
    cycle_amr = argv[5]
    cycle_align = argv[6]
    cycle_sent = argv[7]

    count = 0
    with open(amr_file, 'r') as f:
        with open(cycle_amr, 'w') as wf:
            amr_line = get_amr_line(f)
            while amr_line != '':
                if count in cycle_set:
                    wf.write('%s\n\n' % amr_line)
                amr_line = get_amr_line(f)
                count += 1
            wf.close()
            f.close()

    with open(align_file, 'r') as f:
        with open(cycle_align, 'w') as wf:
            alignments = f.readlines()
            for (i, line) in enumerate(alignments):
                if i in cycle_set:
                    wf.write(line)
            wf.close()
            f.close()

    with open(sent_file, 'r') as f:
        with open(cycle_sent, 'w') as wf:
            sents = f.readlines()
            for (i, line) in enumerate(sents):
                if i in cycle_set:
                    wf.write(line)
            wf.close()
            f.close()

if __name__ == '__main__':
    main(sys.argv)
