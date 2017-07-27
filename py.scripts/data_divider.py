#!/usr/bin/python
import cPickle
import sys
import os
from smatch import get_amr_line
import amr_graph
from amr_graph import *
def main(argv):
    amr_file = argv[1]
    divide_dir = argv[2]
    f = open(amr_file, 'r')
    amr_line = get_amr_line(f)
    graphs = []
    #os.system('rm -rf %s' % divide_dir)
    #os.mkdir(divide_dir)
    curr_index = 0
    while amr_line and amr_line != '':
        #print amr_line
        amr_graph = AMRGraph(amr_line)
        #print >> sys.stderr, str(amr_graph)
        graphs.append(amr_graph)
        if len(graphs) % 5000 == 0:
            curr_dump_file = os.path.join(divide_dir, 'graph_%d' % curr_index)
            curr_f = open(curr_dump_file, 'wb')
            cPickle.dump(graphs, curr_f)
            curr_f.close()
            curr_index += 1
            graphs[:] = []
        amr_line = get_amr_line(f)

    curr_dump_file = os.path.join(divide_dir, 'graph_%d' % curr_index)
    curr_f = open(curr_dump_file, 'wb')
    cPickle.dump(graphs, curr_f)
    curr_f.close()

if __name__ == '__main__':
    main(sys.argv)


    #print 'Total number of amr pairs is %d' % len(graphs)
    #for i in xrange(0, len(graphs), 500):
    #    curr_file_index = i / 500
    #    curr_dump_file = os.path.join(divide_dir, 'graph_%d' % curr_file_index)
    #    curr_f = open(curr_dump_file, 'wb')
    #    end = i + 500
    #    if end > len(graphs):
    #        end = len(graphs)
    #    cPickle.dump(graphs[i:end], curr_f)
    #    curr_f.close()
    print 'Finished dumping all the graphs'

