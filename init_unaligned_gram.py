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
from rule import Rule
from util.hgraph.hgraph import *
if __name__ == '__main__':
    f = open(sys.argv[1], 'r')
    voc_set = set()
    for i, line in enumerate(f):
        words = line.strip().split()
        for w in words:
            voc_set.add(w)

    unalign_f = open(sys.argv[2], 'w')
    for word in voc_set:
        unalign_f.write('[A1] ||| %s ||| (.) ||| 1.0 0.0 0.0\n' % word)
    unalign_f.close()
